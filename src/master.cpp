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

#include "cantera/core.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/ext/bolos/Logger.h"
#include "cantera/kinetics/Boltzmann.h"
#include <iostream>
#include "cantera/zerodim.h"
#include "cantera/numerics/Integrator.h"
#include <cmath>

#include "plasmaReactor.h"
#include "readParameters.h"

using namespace Cantera;

/* ------------------------ PREPARE FEW USEFUL FUNCTIONS ------------------------ */
// Get number density of species [#/cm^3]
double getNumberDens(const shared_ptr<ThermoPhase>& gas, const size_t i ){
    return 1e-6 * gas->moleFraction(i) * Avogadro * gas->molarDensity();
}

// Logarithmically changed timestep
// K: grow/decay order | t_N: cover period of 0-t_N | Nsteps: assigned steps | n: step index
double logDynamicTimestep(const double& K, const double& t_N, const int& Nsteps, const int& n ) {
    return (t_N / K) * log10( (Nsteps + (n+1)*(pow(10, K) - 1)) /
                              (Nsteps + n*(pow(10, K) - 1)) );
}


int main(int argc, char *argv[]) {
    printHeader();
    CppBOLOS::currentLogLevel = CppBOLOS::LOG_INFO;      // Default log level
    std::string controlDictPath = "./controlDict";       // Default path
    std::string parameterPath = "./chemPlasProperties";  // Default path

    // Map string to LogLevel
    std::map<std::string, CppBOLOS::LogLevel> logLevels = {
            {"NONE", CppBOLOS::LOG_NONE},
            {"WARNING", CppBOLOS::LOG_WARNING},
            {"INFO", CppBOLOS::LOG_INFO},
            {"DEBUG", CppBOLOS::LOG_DEBUG}
    };

    // Parse command-line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-case" && i + 1 < argc) {
            std::string basePath = argv[++i];
            controlDictPath = basePath + "/controlDict";
            parameterPath = basePath + "/chemPlasProperties";
        } else if (arg == "-log" && i + 1 < argc) {
            std::string logLevelStr = argv[++i];
            auto it = logLevels.find(logLevelStr);
            if (it != logLevels.end()) {
                CppBOLOS::currentLogLevel = it->second;
            } else {
                std::cerr << "Unknown log level: " << logLevelStr << std::endl;
                return 1;
            }
        }
    }

    /* -------------------------------- SET UP PARAMETERS -------------------------------- */
    // Read time control
    double runTime = readParameter<double>(controlDictPath, "startTime");
    double t_end = readParameter<double>(controlDictPath, "endTime");
    const int nPulses = readParameter<int>(controlDictPath, "nPulses");

    double dt1 = readParameter<double>(controlDictPath, "dt1");
    double dt_max = readParameter<double>(controlDictPath, "dt_max");

    // Set T_gas, E/N, and other NRP parameters
    double Tgas = readParameter<double>(parameterPath, "Temperature");
    double P0 = readParameter<double>(parameterPath, "Pressure");
    double EN = readParameter<double>(parameterPath, "E/N(NRP)");
    double EN_DC = readParameter<double>(parameterPath, "E/N(DC)");
    bool stayEN_DC = readParameter<bool>(parameterPath, "stayDC");
    const double N_e_ini = readParameter<double>(parameterPath, "electronDens");
    double Ep = readParameter<double>(parameterPath, "pulseEnergy");
    const double f_NRP = readParameter<double>(parameterPath, "f_NRP");
    const double pulseWidth = readParameter<double>(parameterPath, "pulseWidth");
    const double pulseExten = readParameter<double>(controlDictPath, "pulseExten");
    const int NstepInPulse = readParameter<int>(controlDictPath, "NstepInPulse");
    const int NstepAfPulse = readParameter<int>(controlDictPath, "NstepAfPulse");
    const double dT_max = readParameter<double>(controlDictPath, "dT_max");
    const double TGAS_TOLERANCE = readParameter<double>(parameterPath, "TGAS_TOLERANCE");
    const double EN_TOLERANCE = readParameter<double>(parameterPath, "EN_TOLERANCE");

    // Default initialization
    double dt = dt_max;                     // initial time step
    double dt3, tau_NRP, nE_3, nE_1, nE_2;  // temporary variables used in the loop
    bool plasmaOn = true;                   // plasma switch flag
    bool dischargeOn = true;                // single nanosecond discharge switch flag
    int iPulse = 0;                         // pulse number index
    double disEnergy = 0.0;                 // initial discharge energy deposited
    const double pulsePeriod = 1.0/f_NRP;
    int iStepInPulse = 0, iStepAfPulse = 0; // step index in/after a pulse discharge
    double T_old = Tgas;
    double dTdt = 0.0, dTdt_max = 1.0;      // initial dT/dt and its maximum
    double time_max_dTdt = 0.0;             // ignition delay time (at maximum dT/dt)
    const double SMALL = 1.0e-16;
    bool fastExpansion = false;             // isentropic gas expansion flag
    t_end = std::max(pulsePeriod * nPulses, t_end); // overload t_end

    /* ------------------------------- SET UP BOLTZMANN SOLVER ------------------------------- */
    std::cout << "\n========  SETTING BOLTZMANN SOLVER ... ========\n" << std::endl;

    // Species configured by CppBOLOS
    auto bSpecies = readParameter<std::vector<std::string>>(parameterPath, "BoltzmannSpecies");
    auto mixture = readParameter<std::map<std::string, double>>(parameterPath, "Mixture");
    const auto FUELNAME = readParameter<string>(parameterPath, "fuelName");
    const double initFuelFraction = mixture[FUELNAME];

    std::map<std::string, double> BoltzmannSpecies;

    for (const auto& species : bSpecies) {
        // Check if species is in the mixture, if not assign 0.0
        double fraction = mixture.count(species) > 0 ? mixture[species] : 0.0;
        BoltzmannSpecies[species] = fraction;
        // Be careful of possible difference of species names in LXCat and YAML input,
        // especially for e/E, He/HE, Ar/AR. Always try to use unified names.
    }

    // Read cross-section data
    auto csDataFile = readParameter<string>(parameterPath, "csDataFile");
    std::stringstream ss = CppBOLOS::clean_file(csDataFile);
    std::vector<CppBOLOS::Collision> collisions = CppBOLOS::parse(ss);

    // Set up grid. This affects accuracy significantly.
    auto gridSize = readParameter<std::map<std::string, double>>(parameterPath, "gridSize");
    BoltzmannRate::bsolver.set_grid
            (
                    readParameter<string>(parameterPath, "gridType"),
                    gridSize["start"],
                    gridSize["end"],
                    (int)gridSize["points"]
            );
    BoltzmannRate::bsolver.load_collisions(collisions);

    LOG_INFO("\nA total of " + std::to_string(BoltzmannRate::bsolver.number_of_targets()) +
             " targets have been loaded:\n" + BoltzmannRate::bsolver.targetNames());

    BoltzmannRate::bsolver.set_kT (Tgas);
    BoltzmannRate::bsolver.set_EN(EN);
    BoltzmannRate::bsolver.set_density(BoltzmannSpecies);
    BoltzmannRate::bsolver.init();

    // Solve EBE, serve as cache and preliminary convergence check
    BoltzmannRate::F0 = BoltzmannRate::bsolver.maxwell(2.0); // initial guess from Maxwell EEDF
    BoltzmannRate::F0 = BoltzmannRate::bsolver.converge(BoltzmannRate::F0, 200, 1e-5);

    /* ----------------------------- SET UP GAS PHASE AND REACTOR  -------------------------------- */
    std::cout << "\n========  SETTING GAS PHASE & REACTOR ... ========\n" << std::endl;

    auto gas_file = readParameter<string>(parameterPath, "mechFile");
    auto sol = newSolution(gas_file, "", "mixture-averaged");
    auto gas_eq = newSolution(gas_file)->thermo();
    auto gas = sol->thermo();

    // Get equilibrate gas temperature.
    gas_eq->setState_TPX(Tgas, P0, mixture);
    gas_eq->equilibrate("HP"); // constant pressure
    auto T_eq = gas_eq->temperature();

    // Read plasma heating model
    std::string plasmaHeatModel = " ";
    try {
        plasmaHeatModel = readParameter<string>(parameterPath, "plasmaHeatModel");
    } catch (const std::runtime_error& e) {
        std::cerr << "Warning: " << e.what() << std::endl;
    }

    // Create a gas-plasma reactor object
    ChemPlasReactor odes = ChemPlasReactor(sol, BoltzmannSpecies, plasmaHeatModel);

    // User can assign electron mole fraction in 'Mixture' to override 'electronDens'
    double gas_numberDensity = 1e-6*P0/Tgas/CppBOLOS::KB;
    if (mixture.count(gas->speciesName(odes.electronIndex)) == 0) {
        mixture[gas->speciesName(odes.electronIndex)] = N_e_ini/gas_numberDensity;
    }
    gas->setState_TPX(Tgas, P0, mixture);
    odes.setConstPD(gas->pressure(), gas->density()); // Record pressure and density

    odes.nonThermal = readParameter<bool>(parameterPath, "nonThermal");
    odes.heatLoss = readParameter<bool>(parameterPath, "heatLoss");
    if (odes.heatLoss) {
        auto C0 = readParameter<double>(parameterPath, "C0");
        odes.T0 = Tgas;
        odes.C0 = C0;
        std::cout << "Start with constant-volume PAC reactor since heat loss model is on."<< std::endl;
        odes.constPressure = false;
    } else
    {
        if (readParameter<string>(parameterPath, "CP_or_CV") == "CP"){
            std::cout << "Start with constant-pressure PAC reactor"<< std::endl;
            odes.constPressure = true;
        } else {
            std::cout << "Start with constant-volume PAC reactor"<< std::endl;
            odes.constPressure = false;
        }
    }
    odes.inertSpIndex = odes.findSpeciesIndex(readParameter<string>(parameterPath, "inertSpecies"));

    /* ------------------------------ CREATE & INIT ODE INTEGRATOR --------------------------------- */
    //  - the default settings for CVodesIntegrator are used:
    //     solution method: BDF_Method
    //     problem type: DENSE + NOJAC
    //     relative tolerance: 1.0e-9
    //     absolute tolerance: 1.0e-15
    //     max step size: +inf
    unique_ptr<Integrator> integrator(newIntegrator("CVODE"));
    // initialize the integrator, specifying the start time and the RHS evaluator object.
    integrator->initialize(runTime, odes);
    auto Tolerances = readParameter<std::map<std::string, double>>(controlDictPath, "odes");
    integrator->setTolerances(Tolerances["reltol"], Tolerances["abstol"]);

    // 0D reactor setup in original Cantera
    IdealGasConstPressureReactor r;
    r.insert(sol); // 'insert' the gas into the reactor and environment.
    ReactorNet sim;
    sim.addReactor(r);
    sim.initialize();
    // sim.setTolerances(1e-9, 1e-18);
    // sim.setVerbose();

    // Access kinetic information
    auto kin = sol->kinetics();

    int irxns = kin->nReactions();
    vector<double> qf(irxns);
    vector<double> qr(irxns);
    vector<double> q(irxns);

    // Species to be stored and their indices
    std::vector<std::string> species_names;
    try {
        species_names = readParameter<std::vector<std::string>>(controlDictPath, "writeSpecies");
    } catch (const std::runtime_error& e) {
        std::cerr << "Warning: " << e.what() << std::endl;
        species_names = gas->speciesNames();
    }
    if (species_names.size() == 1 && species_names[0] == "all") {
        species_names = gas->speciesNames(); // Use all species
    }
    std::vector<int> index_list;

    // Create the csv output file
    std::ofstream outputFile(readParameter<string>(controlDictPath, "outputFile"));
    if (!outputFile.is_open()) {
        std::cerr << "Failed to open the output file. Check if the file path is valid." << std::endl;
        return 1; // Return with error code
    }
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

    // Hold mole/mass fraction data
    std::vector<double> moleFraction_row (index_list.size());

    /* ------------------------------------- RUN THE SIMULATION ---------------------------------- */
    std::cout << "\n======== RUNNING SIMULATION ========\n" << std::endl;

    clock_t t0 = clock(); // save start time

    // Main time loop
    const double dTbelowTeq = readParameter<double>(controlDictPath, "dTbelowTeq");
    while (runTime < t_end && gas->temperature() < (T_eq - dTbelowTeq)) {
        std::ostringstream oss; // Format output digit
        oss << "\nrunTime [s]: " << std::scientific << std::setprecision(9) << runTime;
        std::cout << oss.str() << " | iPulse: " << iPulse << "\n";

        double EN_temp = EN_DC;  // default E/N value

        if (iPulse < nPulses) { // TODO: there could be a more concise way to do dynamic time step
            dt3 = 0.;  // Initialize timestep at discharge quiescent zone
            if ((runTime + SMALL) < (iPulse * pulsePeriod + pulseExten)) { // Stages where dynamic timestep is applied
                if ((runTime + SMALL) < (iPulse * pulsePeriod + pulseWidth) && dischargeOn) {  // during the nanosecond pulse
                    if (iPulse <= 1) { // For the first two pulses, fix dt
                        dt = dt1;
                    } else { // From 3rd pulses, apply logarithmically reduced timestep
                        double K = std::max(std::round(log10(nE_2 / nE_1) * 100) / 100, 0.2);
                        std::cout << "Growth order K = " << K << "; tau_NRP = " << tau_NRP << std::endl;
                        dt = logDynamicTimestep(K, tau_NRP, NstepInPulse, iStepInPulse);
                    }
                    dt = std::max(dt, 1e-14);  // Make sure dt > 1e-14, avoid dt cascade
                    EN_temp = EN;
                    iStepInPulse += 1;
                } else { // After a nanosecond pulse, gradually relax timestep during pulseExten
                    if (iStepAfPulse == 0) { // First step after the discharge is off; record it
                        // sim.reinitialize();
                        integrator->reinitialize(runTime, odes);
                        tau_NRP = (runTime - iPulse * pulsePeriod); // Duration of the nanosecond discharge
                        std::cout << "NSD terminated at " << disEnergy << " [mJ/cm^3]" << ", tau_NRP = " << tau_NRP << "\n";
                        nE_2 = getNumberDens(gas, odes.electronIndex); // Store n_e at the end of a nanosecond pulse

                        // TODO: may need to update Boltzmann grid if E/N changes significantly

                        fastExpansion = true;
                    }
                    double K = (iPulse <= 1) ? -1 : std::min(round(log10(nE_3 / nE_2) * 100) / 100, -0.1); // fixed K for first two pulses
                    std::cout << "Decay order K = " << K << std::endl;
                    dt = (iStepAfPulse >= NstepAfPulse) ? // Make dt unchanged for iStepAfPulse >= NstepAfPulse to cut off dt growth
                         dt : logDynamicTimestep(K, pulseExten, NstepAfPulse, iStepAfPulse);

                    // Isentropic heat loss model. P, T, rho jump at 100 ns after the pulse
                    if (runTime > (iPulse * pulsePeriod + 1e-7) && odes.heatLoss && fastExpansion )
                    {
                        double T_new = gas->temperature() * pow( P0/gas->pressure(), 1 - gas->cv_mass()/gas->cp_mass() );
                        std::cout << "Temperature jump: " << gas->temperature() << " -> " << T_new << "\n";
                        gas->setState_TPY(T_new, P0, gas->massFractions());
                        odes.setConstPD(gas->pressure(), gas->density());
                        odes.constPressure = true;
                        fastExpansion = false; // only do once for a pulse
                    }

                    iStepAfPulse += 1;
                    if (iStepAfPulse == NstepAfPulse) { // At the end of pulseExten, record n_e
                        nE_3 = getNumberDens(gas, odes.electronIndex);
                    }
                }
            } else { // After pulseExten stage, using constant dt = dt_max
                dt3 = (iPulse + 1) * pulsePeriod - runTime;  // Avoid overstepping
                dt = std::min(dt3, dt_max);
            }
        } else { // NRP discharges finished; simply apply dt_max as tentative dt value
            if (not stayEN_DC) EN_temp = 0.0;
            dt = dt_max;
        }

        // Update dTdt_max and its time
        if (dTdt > dTdt_max & iStepAfPulse > 1) {
            dTdt_max = dTdt;
            time_max_dTdt = runTime;
        }

        // Timestep should also be reduced to capture steep temperature gradient
        double dt_ig = dt_max;  // Define dynamic step to capture ignition
        if (dTdt > 1000) dt_ig = dT_max / dTdt;
        dt_ig = std::min(dt_ig, dt_max);  // Make sure dt_ig don't exceed dt_max

        // In case the ignition happens at the end of a pulse period, make sure no overstepping
        if (iPulse < nPulses && 0 < dt3 && dt3 < dt_ig) dt_ig = dt3;

        // ======== Finally dt is determined =======
        dt = std::min(dt, dt_ig);

        if (dt < 1e-16) {
            std::cerr << "dt is too small !" << "\n";
            abort();
        }
        std::cout << "dt = " << dt << "; iStepInPulse = " << iStepInPulse << "; iStepAfPulse = " << iStepAfPulse
                  << std::endl;

        // Update Boltzmann solver as needed and advance time step.
        // This could be costly and should be optimized
        bool updateBoltzmannSolver = false;
        if (odes.updateBoltzmannMixture()) {
            BoltzmannRate::bsolver.set_density(BoltzmannSpecies);
            updateBoltzmannSolver = true;
        }
        if (std::abs(EN_temp - BoltzmannRate::bsolver.get_EN()) > EN_TOLERANCE) {
            BoltzmannRate::bsolver.set_EN(EN_temp);
            updateBoltzmannSolver = true;
        }
        if(std::abs(gas->temperature() - BoltzmannRate::bsolver.get_kT()) > TGAS_TOLERANCE) {
            BoltzmannRate::bsolver.set_kT(gas->temperature());
            updateBoltzmannSolver = true;
        }
        if (updateBoltzmannSolver) {
            BoltzmannRate::bsolver.init();
            std::cout << "Updating EEDF ..." << "\n";
            BoltzmannRate::updateBoltzmannSolver();
            integrator->reinitialize(runTime, odes);
            // BoltzmannRate::updateBoltzmannSolver(300, 1e-6); // update with new iteration settings
        }
        std::cout << "EN = " << BoltzmannRate::bsolver.get_EN() << "; Tgas = " << gas->temperature()
                  << "; Te = " << BoltzmannRate::bsolver.get_Te() << "\n";

        /*----------------------  Advance time step --------------------*/
        runTime += dt;
        // sim.advance(runTime);
        integrator->integrate(runTime);
        odes.updateState(integrator->solution());
        /*--------------------------------------------------------------*/

        // Print out information and store data in csv
        std::cout << "Writing number density [#/cm^3] of " << gas->speciesName(index_list[0]) << ", "
                  << gas->speciesName(index_list[1]) << ", " << gas->speciesName(index_list[1]) << "... " << std::endl;
        moleFraction_row = {};
        outputFile << runTime << ", " << gas->temperature() << ", " << 1e-6*Avogadro*gas->molarDensity();// r.temperature();
        for (const auto &i: index_list) {
            auto dens = getNumberDens(gas, i);
            moleFraction_row.push_back(dens); // Store number density [#/cm^3]
            std::cout << dens << " ";
            outputFile << ", " << dens;
        }
        outputFile << std::endl;
        std::cout << std::endl;
        std::cout << "Conversion ratio %: " << 100*(1 - gas->moleFraction(FUELNAME)/initFuelFraction) << std::endl;

        // Check plasma power and energy deposited
        double power_temp = BoltzmannRate::bsolver.elec_power()
                            * Avogadro * gas->molarDensity()* ElectronCharge * getNumberDens(gas, odes.electronIndex);
        // unit: eV m^3/s * #/kmol * kmol/m^3 * coulomb(J/eV) * #/cm^3 = J/s/cm^3
        disEnergy = 1e-3 * odes.depositedPlasmaEnergy(); // mJ/cm^3
        std::cout << "Power_input [W/cm^3]: " << power_temp << " | Energy deposited [mJ/cm^3]: " << disEnergy << "\n";
        // double P_elec = 1.5 * GasConstant * BoltzmannRate::bsolver.get_Te();

        // Check if terminate nanosecond discharge, cut off by designated power deposition Ep
        dischargeOn = (disEnergy < Ep);

        // Update dT/dt
        dTdt = (gas->temperature() - T_old) / dt;
        T_old = gas->temperature();

        // Trigger plasmaOn only ONCE for a new pulse
        if (runTime + SMALL >= (iPulse + 1) * pulsePeriod ) { // Entering next pulse period
            iPulse += 1;
            if (iPulse == nPulses) {
                std::cout << "\n >>>>>>>> NRP Discharges Finished <<<<<<<< \n" << std::endl;
                integrator->reinitialize(runTime, odes);
                plasmaOn = false;
            }
            if (plasmaOn) {
                std::cout << "\n**** Start Pulse Loop No. " << iPulse + 1 << " at runTime = " << runTime << " ****\n";
                if (odes.heatLoss) {
                    odes.setConstPD(gas->pressure(), gas->density());
                    odes.constPressure = false; // switch to constant volume
                }

                disEnergy = 0.0, dischargeOn = true;
                odes.resetDepositedPlasmaEnergy();
                iStepInPulse = 0, iStepAfPulse = 0;
                nE_1 = getNumberDens(gas, odes.electronIndex); // n_e at the beginning of a new pulse
            }
        }

        std::cout << "-------------" << "\n";
    }
    // End of the main time loop

    std::cout << "\nIgnition Delay Time [s]: " << time_max_dTdt << std::endl;

    clock_t t1 = clock(); // save end time
    double elapsed_time = static_cast<double>(t1 - t0) / CLOCKS_PER_SEC;
    std::cout << "\n======== ChemPlasKin finished. Execution time: " << elapsed_time << " (sec) ========\n" <<std::endl;


    return 0;

}