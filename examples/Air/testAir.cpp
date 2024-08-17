/*--------------------------------*- C++ -*----------------------------------*\
|     ____ _                    ____  _           _  ___                      |
|    / ___| |__   ___ _ __ ___ |  _ \| | __ _ ___| |/ (_)_ __                 |
|   | |   | '_ \ / _ \ '_ ` _ \| |_) | |/ _` / __| ' /| | '_ \                |
|   | |___| | | |  __/ | | | | |  __/| | (_| \__ \ . \| | | | |               |
|    \____|_| |_|\___|_| |_| |_|_|   |_|\__,_|___/_|\_\_|_| |_|               |
|                                                                             |
|   A Freeware for Unified Gas-Plasma Kinetics Simulation                     |
|   Version:      1.0.0 (July 2024)                                           |
|   License:      GNU LESSER GENERAL PUBLIC LICENSE, Version 2.1              |
|   Author:       Xiao Shao                                                   |
|   Organization: King Abdullah University of Science and Technology (KAUST)  |
|   Contact:      xiao.shao@kaust.edu.sa                                      |
\*---------------------------------------------------------------------------*/

// Test case of air plasma kinetics using external profiles of electron density and electric field
// Reference: https://doi.org/10.1088/0022-3727/46/46/464010; https://doi.org/10.1016/j.combustflame.2022.111990

#include "cantera/core.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/ext/bolos/Logger.h"
#include "cantera/kinetics/Boltzmann.h"
#include "cantera/zerodim.h"
#include "cantera/numerics/Integrator.h"

#include "../../src/plasmaReactor.h"

using namespace Cantera;

/* ------------------------ PREPARE FEW USEFUL FUNCTIONS ------------------------ */
// Get number density of species [#/cm^3]
double getNumberDens(const shared_ptr<ThermoPhase>& gas, const size_t i ){
    return 1e-6 * gas->moleFraction(i) * Avogadro * gas->molarDensity();
}

// Function to read CSV data
void readCSV(std::string fileName, std::vector<std::pair<double, double>> &y_t) {
    std::ifstream file(fileName);
    // Check if the file is open
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + fileName);
    }
    std::string line, value;

    while (getline(file, line)) {
        std::stringstream ss(line);
        double t_val, y_val;

        getline(ss, value, ',');
        t_val = stod(value);

        getline(ss, value, ',');
        y_val = stod(value);

        y_t.emplace_back(t_val, y_val);
    }
}

// Quadratic interpolation function
double interpolate(const std::vector<std::pair<double, double>> &y_t, double queryPoint) {
    if (y_t.size() < 3) {
        throw std::runtime_error("Need at least 3 data points for quadratic interpolation.");
    }

    for (size_t i = 0; i < y_t.size() - 2; ++i) {
        if (queryPoint >= y_t[i].first && queryPoint <= y_t[i + 2].first) {
            // Three points for quadratic interpolation
            double x0 = y_t[i].first, y0 = y_t[i].second;
            double x1 = y_t[i + 1].first, y1 = y_t[i + 1].second;
            double x2 = y_t[i + 2].first, y2 = y_t[i + 2].second;

            // Coefficients of the quadratic polynomial ax^2 + bx + c
            double a, b, c;

            double denom = (x0 - x1) * (x0 - x2) * (x1 - x2);
            a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / denom;
            b = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / denom;
            c = (x1 * x2 * (x1 - x2) * y0 + x2 * x0 * (x2 - x0) * y1 + x0 * x1 * (x0 - x1) * y2) / denom;

            // Interpolated value
            return a * queryPoint * queryPoint + b * queryPoint + c;
        }
    }

    std::cerr << "Warning: Query point '" << queryPoint << "' is out of the data range! ";
    double boundaryValue;
    if (queryPoint < y_t[0].first) {
        boundaryValue = y_t[0].second;
    } else {
        boundaryValue = y_t[y_t.size() - 1].second;
    }
    std::cerr << "Returning " << boundaryValue << std::endl;

    return boundaryValue;
}

int main()
{
    CppBOLOS::currentLogLevel = CppBOLOS::LOG_INFO;

    // Set up time control and pulse number
    double runTime = 0.0;
    double dt = 0.1E-9;
    double endTime = 100E-9;
    bool thermalEffect = true;
    bool printReactionRates = false;

    /* ------------------ READ EXPERIMENTAL DATA ------------------ */
    std::vector<std::pair<double, double>> Vp_t_data;
    std::vector<std::pair<double, double>> Ne_t_data;
    // Read data from CSV
    readCSV("../Vp-t.csv", Vp_t_data);
    readCSV("../ne-t.csv", Ne_t_data);

    // Example usage of interpolate function
    double queryPoint = 5.0;  // Example query point
    try {
        double interpolatedVp = interpolate(Vp_t_data, queryPoint);
        double interpolatedNe = interpolate(Ne_t_data, queryPoint);
        std::cout << "Interpolated Vp(kV) value at t = " << queryPoint << " is: " << interpolatedVp << std::endl;
        std::cout << "Interpolated Ne(#/cm^3) value at t = " << queryPoint << " is: " << interpolatedNe << std::endl;
    } catch (const std::runtime_error &e) {
        std::cerr << e.what() << std::endl;
    }

    /* ------------------------------- SET UP BOLTZMANN SOLVER ------------------------------- */
    std::cout << "\n========  SETTING BOLTZMANN SOLVER ... ========\n" << std::endl;    // Species configured by CppBOLOS

    // Species configured by CppBOLOS
    std::map<std::string, double> BoltzmannSpecies = {
            // Be careful of possible difference of species names in LXCat and *yaml input,
            // especially for e/E, He/HE, Ar/AR. Always try to use unified names.
            {"N2", 0.774},
            {"O2", 0.186},
            {"O", 0.04},
            {"N2(A3)", 0}, {"N2(B3)", 0}, {"N2(a1)", 0}, {"N2(C3)", 0},
            {"O2(a1)", 0}, {"O(1D)", 0},
            {"N", 0}, {"N(2D)", 0},
            {"NO", 0},
            {"O3", 0},
//            {"CH4", 0},
//            {"Ar", 0},
    };

    // Read cross-section data
    std::string CS_data_file = "../../../data/LXCat/bolsigdb_air_NH3_H2.dat";
    std::stringstream ss = CppBOLOS::clean_file(CS_data_file);
    std::vector<CppBOLOS::Collision> collisions = CppBOLOS::parse(ss); // parse collision data

    // Set up grid. This can affect accuracy.
    BoltzmannRate::bsolver.set_grid("QuadraticGrid", 0, 60, 150);
    BoltzmannRate::bsolver.load_collisions(collisions);
    LOG_INFO("\nA total of " + std::to_string(BoltzmannRate::bsolver.number_of_targets()) +
             " targets have been loaded:\n" + BoltzmannRate::bsolver.targetNames());

    // Set T_gas, E/N, density and Initialization
    double Tgas = 1500; // gas teperature [K]
    double EN = 150; // reduced electirc field [Td]
    double nE = 1e13; // initial electron number density [#/cm^3]
    BoltzmannRate::bsolver.set_kT (Tgas);
    BoltzmannRate::bsolver.set_EN(EN);
    BoltzmannRate::bsolver.set_density(BoltzmannSpecies);
    BoltzmannRate::bsolver.init();

    // Solve EBE, serve as cache and preliminary convergence check
    BoltzmannRate::F0 = BoltzmannRate::bsolver.maxwell(4.0); // initial guess from Maxwell EEDF

    try{
        BoltzmannRate::F0 = BoltzmannRate::bsolver.converge(BoltzmannRate::F0, 200, 1e-5);
    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return -1;
    }

    double mean_energy = BoltzmannRate::bsolver.mean_energy(BoltzmannRate::F0);
    std::cout << "mean energy: " << mean_energy << "eV (" << BoltzmannRate::bsolver.get_Te() << "K)" << std::endl;

    /* ----------------------------- SET UP GAS PHASE AND REACTOR  -------------------------------- */
    std::cout << "\n========  SETTING GAS PHASE & REACTOR ... ========\n" << std::endl;

    auto sol = newSolution("../../../data/PAC_kinetics/gri30_plasma.yaml");
    std::cout <<"counting Boltzmann Processes: " << BoltzmannRate::NumProcess << std::endl;

    auto gas = sol->thermo();
    double P0 = OneAtm;

    double gas_density = 1e-6*P0/Tgas/CppBOLOS::KB; // [#/cm^3]

    // Create a mixture composition map
    Composition compMap = BoltzmannSpecies;
    compMap["Electron"] = nE/gas_density;
    // compMap["N2+"] = nE;

    gas->setMoleFractionsByName(compMap);
    gas->setState_TP(Tgas, P0);

    // Set the state
    // Initial electron number density does not matter; too small value may cause problem.
    int nsp = gas->nSpecies();
    std::cout<< "initial density: " << getNumberDens(gas, gas->speciesIndex("Electron")) << ", "<< gas->moleFraction("N2+") << " nSpecies: " << nsp << std::endl;

    // Create a gas-plasma reactor object
    ChemPlasReactor odes = ChemPlasReactor(sol, BoltzmannSpecies);
    odes.inertSpIndex = odes.findSpeciesIndex("N2");
    odes.nonThermal = !thermalEffect;
    odes.constPressure = false;

    /* ------------------------------ CREATE & INIT ODE INTEGRATOR --------------------------------- */
    unique_ptr<Integrator> integrator(newIntegrator("CVODE"));
    integrator->initialize(runTime, odes);
    // integrator->setTolerances(1e-5, 1e-10);

    // Access kinetic information
    auto kin = sol->kinetics();

    int irxns = kin->nReactions();
    vector<double> qf(irxns);
    vector<double> qr(irxns);
    vector<double> q(irxns);

    std::vector<std::string> species_names = gas->speciesNames(); // {"e", "N2(A)", "N2^+"};
    std::vector<int> index_list;

    // Create the csv output file
    std::ofstream outputFile("../outputAir.csv");
    outputFile << "Time(s), T_gas(K), N_gas(#/cm^3)";
    std::cout << "Writing information of: \n" << "Time, T_gas(K), N_gas(#/cm^3)";
    for (const auto& sp : species_names) {
        std::cout << ", " << sp;
        size_t index = gas->speciesIndex(sp);
        if (index < std::numeric_limits<size_t>::max()) {
            index_list.push_back(gas->speciesIndex(sp));
        } else {
            throw std::runtime_error("No valid species name found: " + sp);
        }
        outputFile << ", " << sp;
    }
    outputFile << std::endl;
    std::cout << std::endl;

    std::cout << "Time\t";
    for (const auto& name : species_names) {
        std::cout << name << "\t";
    }
    std::cout << std::endl;

    // Hold mole/mass fraction data
    std::vector<double> moleFraction_row (index_list.size());

    /* ------------------------------------- RUN THE SIMULATION ---------------------------------- */
    std::cout << "\n======== RUNNING SIMULATION ========\n" << std::endl;

    clock_t t0 = clock(); // save start time

    // Main time loop
    double disEnergy = 0.0;                 // initial discharge energy deposited
    while (runTime < endTime){

        // Relax dt after 25ns
        if (runTime > 25e-9) {
            dt = 5e-10;
        }

        std::cout << "\nrunTime [s]: "  << runTime << ", dt = " << dt << std::endl;

        // -------------------- Update E/N and nE from experimental data ----------------------
        try {
            EN = std::max(interpolate(Vp_t_data, runTime*1e9) * 1000 / (4e-3 * Avogadro*gas->molarDensity()) * 1e21, 0.1);
            nE = interpolate(Ne_t_data, runTime*1e9);
            std::cout << "Interpolated E/N = " <<
                      EN << ", n_E = " <<  nE << std::endl;
        } catch (const std::runtime_error &e) {
            std::cerr << e.what() << std::endl;
        }

        // Update E/N and mixture for the Boltzmann solver
        BoltzmannRate::bsolver.set_EN(EN);
        if (odes.updateBoltzmannMixture()) {
            BoltzmannRate::bsolver.set_density(BoltzmannSpecies);
        }

        std::cout << "Updating EEDF ..." << "\n";
        BoltzmannRate::bsolver.init();
        BoltzmannRate::updateBoltzmannSolver(200, 1e-5, 1E20/(EN*EN));
        // **NOTE**: 1E20/(EN*EN) is a good estimation of delta0 to speed up bsolver convergence,
        // especially for significantly varying E/N values.
        // 200: maximum iteration number; 1e-5: convergence tolerance
        // You may simply use default parameters by: BoltzmannRate::updateBoltzmannSolver()

        std::cout << "EN = " << BoltzmannRate::bsolver.get_EN() << ", Tgas = " << gas->temperature() << "[K]; Te = " << BoltzmannRate::bsolver.get_Te() << "[K]\n";

        // Update the integrator
        odes.imposeNe(nE);
        integrator->reinitialize(runTime, odes);

        /*----------------------  Advance time step --------------------*/
        runTime += dt;
        integrator->integrate(runTime);
        odes.updateState(integrator->solution());
        std::cout << "Internal integration steps count: " << odes.nSteps << std::endl;
        odes.nSteps = 0;
        /*--------------------------------------------------------------*/

        // Print out information and store data in csv
        std::cout << "Writing number density [#/cm^3] of " << gas->speciesName(index_list[0]) << ", "
                  << gas->speciesName(index_list[1]) << ", " << gas->speciesName(index_list[2]) << ", ... " << std::endl;
        moleFraction_row = {};
        outputFile << runTime << ", " << gas->temperature() << ", " << 1e-6*Avogadro*gas->molarDensity();
        for (const auto &i: index_list) {
            auto dens = getNumberDens(gas, i);
            moleFraction_row.push_back(dens); // Store number density [#/cm^3]
            std::cout << dens << " ";
            outputFile << ", " << dens;
        }
        outputFile << std::endl;
        std::cout << std::endl;

        // Check plasma power and energy deposited
        double power_temp = BoltzmannRate::bsolver.elec_power()
                            * Avogadro * gas->molarDensity()* ElectronCharge * getNumberDens(gas, odes.electronIndex);
        // unit: eV m^3/s * #/kmol * kmol/m^3 * coulomb(J/eV) * #/cm^3 = J/s/cm^3
        disEnergy = 1e-3 * odes.depositedPlasmaEnergy(); // mJ/cm^3
        std::cout << "Plasma Power [W/cm^3]: " << power_temp << " | Energy Deposition [mJ/cm^3]: " << disEnergy << "\n";

        odes.printTopRatesOfProgress();
        odes.printTopSpeciesProductionRates();

        if(printReactionRates){
            kin->getFwdRatesOfProgress(&qf[0]);
            kin->getRevRatesOfProgress(&qr[0]);
            kin->getNetRatesOfProgress(&q[0]);

            writelog("{:30s} {:>14s} {:>14s} {:>14s}  {:s}\n",
                     "Reaction", "Forward", "Reverse", "Net", "Unit");
            for (int i = 0; i < irxns; i++) {
                const auto& rxn = kin->reaction(i);
                writelog("{:30s} {:14.5g} {:14.5g} {:14.5g}  kmol/m3/s\n",
                         rxn->equation(), qf[i], qr[i], q[i]);
            }
        }

        std::cout << "-------------" << "\n";
    }

    clock_t t1 = clock(); // save end time
    double elapsed_time = static_cast<double>(t1 - t0) / CLOCKS_PER_SEC;
    std::cout << "\n======== ChemPlasKin finished. Execution time: " << elapsed_time << " (sec) ========\n" <<std::endl;

    return 0;

}
