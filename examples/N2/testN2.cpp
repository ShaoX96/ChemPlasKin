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

// Test case of N2 plasma kinetics using external profiles of electron density and electric field
// Reference: http://www.zdplaskin.laplace.univ-tlse.fr/external-profiles-of-electron-density-and-electric-field/index.html

#include "cantera/core.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/ext/bolos/Logger.h"
#include "cantera/kinetics/Boltzmann.h"
#include "cantera/zerodim.h"
#include "cantera/numerics/Integrator.h"

#include "plasmaReactor.h"

using namespace Cantera;

/* ------------------------ PREPARE FEW USEFUL FUNCTIONS ------------------------ */
// Get number density of species [#/cm^3]
double getNumberDens(const shared_ptr<ThermoPhase>& gas, const size_t i ){
    return 1e-6 * gas->moleFraction(i) * Avogadro * gas->molarDensity();
}

int main()
{
    CppBOLOS::currentLogLevel = CppBOLOS::LOG_INFO;

    // Set up time control and pulse number
    double runTime = 0.0;
    int nPulse = 10;
    bool thermalEffect = true;
    bool printReactionRates = true;

    /* ------------------------------- SET UP BOLTZMANN SOLVER ------------------------------- */
    std::cout << "\n========  SETTING BOLTZMANN SOLVER ... ========\n" << std::endl;    // Species configured by CppBOLOS

    std::cout << "Read external profiles of electron density and electric field" << std::endl;
    std::ifstream file("../data_in.dat");

    if (!file) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }

    std::string line;
    std::vector<std::vector<double>> data_in(3);
    int lineNumber = 0;

    // Skip the header line
    std::getline(file, line);
    lineNumber++;

    while (std::getline(file, line)) {
        lineNumber++;
        std::istringstream iss(line);
        double t, f, e;

        if (!(iss >> t >> f >> e)) {
            std::cerr << "Error reading line " << lineNumber << ": " << line << std::endl;
            continue;
        }

        data_in[0].push_back(t);
        data_in[1].push_back(f);
        data_in[2].push_back(e);
    }
    file.close();

    // Test reading by printing the data
    for (size_t i = 0; i < data_in[0].size(); ++i) {
        std::cout << "Time: " << data_in[0][i]
                  << " s, Field: " << data_in[1][i]
                  << " Td, Electrons: " << data_in[2][i]
                  << " cm^-3" << std::endl;
    }

    // Species configured by CppBOLOS
    std::map<std::string, double> BoltzmannSpecies = {
            {"N2", 1.0}
    };

    // Read cross-section data
    std::string filename = "../../../data/LXCat/bolsigdb_N2.dat";
    std::stringstream ss = CppBOLOS::clean_file(filename);
    std::vector<CppBOLOS::Collision> collisions = CppBOLOS::parse(ss);

    // Set up grid.
    BoltzmannRate::bsolver.set_grid("QuadraticGrid", 0, 30, 180);
    BoltzmannRate::bsolver.load_collisions(collisions);
    LOG_INFO("\nA total of " + std::to_string(BoltzmannRate::bsolver.number_of_targets()) +
             " targets have been loaded:\n" + BoltzmannRate::bsolver.targetNames());

    // Set T_gas, E/N, density and initialization
    double Tgas = 300; // gas teperature [K]
    double EN = 100; // reduced electirc field [Td]
    double nE = 1e6; // initial electron number density [#/cm^3]
    BoltzmannRate::bsolver.set_kT (Tgas);
    BoltzmannRate::bsolver.set_EN(EN);
    BoltzmannRate::bsolver.set_density(BoltzmannSpecies);
    BoltzmannRate::bsolver.init();

    // Solve EBE, serve as cache and preliminary convergence check
    BoltzmannRate::F0 = BoltzmannRate::bsolver.maxwell(2.0); // initial guess from Maxwell EEDF

    try{
        BoltzmannRate::F0 = BoltzmannRate::bsolver.converge(BoltzmannRate::F0, 200, 1e-5);
    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return -1;
    }

    /* ----------------------------- SET UP GAS PHASE AND REACTOR  -------------------------------- */
    std::cout << "\n========  SETTING GAS PHASE & REACTOR ... ========\n" << std::endl;

    auto sol = newSolution("../../../data/PAC_kinetics/plasmaN2.yaml");
    std::cout <<"counting Boltzmann Processes: " << BoltzmannRate::NumProcess << std::endl;

    auto gas = sol->thermo();
    double P0 = OneAtm;

    double gas_density = 1e-6*P0/Tgas/CppBOLOS::KB; // [#/cm^3]

    // Create a mixture composition map
    Composition compMap;
    compMap["N2"] = gas_density;
    compMap["e"] = nE;
    compMap["N2^+"] = nE;

    // Set the state
    gas->setMoleFractionsByName(compMap);
    gas->setState_TP(Tgas, P0);

    // Initial electron number density does not matter; too small value may cause problem.
    int nsp = gas->nSpecies();
    std::cout<< "initial density: " << getNumberDens(gas, gas->speciesIndex("e")) << ", "<< gas->moleFraction("N2^+") << " nSpecies: " << nsp << std::endl;

    // Create a gas-plasma reactor object
    ChemPlasReactor odes = ChemPlasReactor(sol, BoltzmannSpecies);
    odes.inertSpIndex = odes.findSpeciesIndex("N2");
    odes.nonThermal = !thermalEffect;

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
    std::ofstream outputFile("../outputN2.csv");
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

    for(int iPulse=0; iPulse<nPulse; iPulse++ ){

        for(int i=0; i<data_in[0].size()-1; i++){
            double dt = data_in[0][i+1] -  data_in[0][i];
            EN = 0.5*(data_in[1][i+1] +  data_in[1][i]);
            nE = 0.5*(data_in[2][i+1] +  data_in[2][i]);

            std::ostringstream oss; // Format output digit
            oss << "\nrunTime [s]: " << std::scientific << std::setprecision(9) << runTime+dt;
            std::cout << oss.str() << " | iPulse: " << iPulse << "\n";
            std::cout << "dt = " << dt << std::endl;

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

            std::cout << "E/N = " << BoltzmannRate::bsolver.get_EN() << "[Td]; Tgas = " << gas->temperature()
                      << "[K]; Te = " << BoltzmannRate::bsolver.get_Te() << "[K]\n";

            // Update the integrator
            odes.imposeNe(nE);
            integrator->reinitialize(runTime, odes);

            /*----------------------  Advance time step --------------------*/
            runTime += dt;
            integrator->integrate(runTime);
            odes.updateState(integrator->solution());
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
    }

    clock_t t1 = clock(); // save end time
    double elapsed_time = static_cast<double>(t1 - t0) / CLOCKS_PER_SEC;
    std::cout << "\n======== ChemPlasKin finished. Execution time: " << elapsed_time << " (sec) ========\n" <<std::endl;

    return 0;

}
