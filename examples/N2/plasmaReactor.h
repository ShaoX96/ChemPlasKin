// ChemPlasReactor class. Generate a unified gas-plasma ODE system

#ifndef CHEMPLASKIN_PLASMAREACTOR_H
#define CHEMPLASKIN_PLASMAREACTOR_H

#include "cantera/numerics/FuncEval.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/kinetics/Boltzmann.h"

namespace Cantera
{
class ChemPlasReactor : public FuncEval {
public:
    /*
     * Constructor
     * @param[in] sol Solution object specifying initial system state.
     */
    ChemPlasReactor
    (
            shared_ptr<Solution> sol,
            std::map<std::string, double>& BoltzmannSpecies,
            std::string plasmaHeatModel = "default"
    )
    : m_BoltzmannSpecies(BoltzmannSpecies)
    {
        /* ---------------------- INITIALIZE MEMBER VARS ---------------------- */
        // Pointer to the system's ThermoPhase object, kinetics manager and Transport object.
        // Updated by the solver during simulation to provide iteration-specific information
        m_gas = sol->thermo();
        m_kinetics = sol->kinetics();

        // Record pressure and density at initial state, which can be reset using setConstPD()
        m_pressure = m_gas->pressure(); // For constant pressure case (CP)
        m_density = m_gas->density(); // For constant volume case (CV)

        // Number of chemical species in the system.
        m_nSpecies = m_gas->nSpecies();
        std::cout <<"Creating ChemPlasReactor: " << m_nSpecies  << " species and "
        << m_kinetics->nReactions() << " reactions." << std::endl;

        // Resize the vector<double> storage containers for
        // species partial molar enthalpies/internal energy
        // and net production rates. Updated and used by the solver per iteration.
        m_hbar.resize(m_nSpecies);
        m_ubar.resize(m_nSpecies);
        m_wdot.resize(m_nSpecies);
        m_q_net.resize(m_kinetics->nReactions());
        m_deltaH.resize(m_kinetics->nReactions());

        // Number of equations in the ODE system.
        m_nEqs = m_nSpecies + 4;

        electronIndex = findSpeciesIndex({"E", "electron", "Electron", "e"});

        // Specify plasma heating model.
        // Detailed: calculating from each reaction steps.
        // Flitti & Pancheshnyi 2009 model: global modelling
        detailPlasmaHeatModel = (plasmaHeatModel != "Flitti");

        // Register energy transfer data for plasma reactions
        setPlasmaEnergyTransfer(sol->reaction_section);

        // Create a map for all species
        for (size_t i = 0; i < m_gas->nSpecies(); i++) {
            allSpeciesIndices[m_gas->speciesName(i)] = i;
        }

        // Register species configured by CppBOLOS
        for (const auto& [sp, _] : m_BoltzmannSpecies) {
            bSpeciesIndices[sp] = findSpeciesIndex(sp);
            allSpeciesIndices.erase(sp);  // Remove the species from the all-species map
        }
    }

    /**
     * Evaluate the ODE right-hand-side function, ydot = f(t,y).
     *   - overridden from FuncEval, called by the integrator during simulation.
     * @param[in] t time.
     * @param[in] y solution vector, length neq()
     * @param[out] ydot rate of change of solution vector, length neq()
     * @param[in] p sensitivity parameter vector, length nparams()
     *   - note: sensitivity analysis isn't implemented in this example
     */
    void eval(double t, double* y, double* ydot, double* p) override {
        // Time derivative of the solution vector [T, he, Ep, e_vib(N2), Y_1, Y_2, ... Y_K],
        // *ydot*, are defined for clear and convenient access to these vectors:
        nSteps ++;
//        std::cout << std::scientific << std::setprecision(10);
//        std::cout << "eval t: " << t << std::endl;
        double temperature = y[0];
        // std::cout << "evl::temperature: " << temperature << std::endl;
        m_he = y[1];
        m_Ep = y[2];
        eVibN2 = y[3];
        double *massFracs = &y[4];

        // Apply external electron number density profile
        if (imposedNe) {
            massFracs[electronIndex] = m_N_e*1e6/Avogadro*m_gas->molecularWeight(electronIndex) / m_density;
        }
        // std::cout << "Debug: ne=" << massFracs[electronIndex] << std::endl;

        double *dTdt = &ydot[0];
        double *HE_dot = &ydot[1];
        double *Ep_dot = &ydot[2];
        double *eVibN2_dot = &ydot[3];
        double *dYdt = &ydot[4];

        /* ------------------------- UPDATE GAS STATE ------------------------- */
        // The state of the ThermoPhase is updated to reflect the current solution
        // vector, which was calculated by the integrator.
        m_gas->setMassFractions(massFracs);
        if (constPressure) {
            m_gas->setState_TP(temperature, m_pressure); // Constant pressure
        } else {
            m_gas->setState_TD(temperature, m_density); // Constant volume
        }

        /* ----------------------- GET REQ'D PROPERTIES ----------------------- */
        double rho = m_gas->density();
        double cp = m_gas->cp_mass();
        double cv = m_gas->cv_mass();
        m_gas->getPartialMolarEnthalpies(&m_hbar[0]); // J/kmol
        m_gas->getPartialMolarIntEnergies(&m_ubar[0]); // J/kmol
        m_kinetics->getNetProductionRates(&m_wdot[0]); // kmol/m^3/s or kmol/m^2/s
        m_kinetics->getDeltaEnthalpy(&m_deltaH[0]); // J/kmol

        /* -------------------------- ENERGY EQUATION ------------------------- */
        // the rate of change of the system temperature is found using the energy
        // equation for a closed-system constant pressure ideal gas:
        //     m*cp*dT/dt = - sum[h(k) * dm(k)/dt] + plasmaHeatingSource + heatLoss

        if (nonThermal) {
            *dTdt = 0.0;
            *eVibN2_dot = 0.0;
            *HE_dot = 0.0;
        } else {
            // unit: [J/m^3/s]
            double hdot_vol = 0;
            double udot_vol = 0;
            for (size_t k = 0; k < m_nSpecies; k++) {
                hdot_vol += m_hbar[k] * m_wdot[k];
                udot_vol += m_ubar[k] * m_wdot[k];
            }

            if (detailPlasmaHeatModel) {
                // Calculate deposited plasma energy using sum{\epsilon_th^j * k_j}
                double eth_sum = 0.0;
                double ext_sum = 0.0; // Electron energy deposited
                double evibN2_sum = 0.0;
                double FGH_sum = 0.0;

                m_kinetics->getNetRatesOfProgress(&m_q_net[0]);
                for (const auto &[i, e_th]: e_th_map) {
                    // Unit: J/m^3/s (= eV * kmol/m^3/s * #/kmol * J/eV)
                    eth_sum += e_th * m_q_net[i] * Avogadro * ElectronCharge;
                    FGH_sum += (e_th * Avogadro * ElectronCharge - m_deltaH[i]) * m_q_net[i];
                }
                for (const auto &[i, e_ext]: e_ext_map) {
                    ext_sum += e_ext * m_q_net[i] * Avogadro * ElectronCharge;
                }
                for (const auto &[i, e_vib]: evib_N2_map) {
                    evibN2_sum += e_vib * m_q_net[i] * Avogadro * ElectronCharge;
                }

                *eVibN2_dot = 0;
                *Ep_dot = ext_sum;

                // Currently we only consider e_vib of N2.
                if (constPressure) {
                    // *dTdt = - hdot_vol / (rho * cp);
                    *HE_dot = -hdot_vol + eth_sum ;
                    *dTdt = (-hdot_vol + eth_sum) / (rho * cp);
                } else {
                    // *dTdt = - udot_vol / (rho * cv);
                    *HE_dot = -udot_vol + eth_sum ;
                    *dTdt = (-udot_vol + eth_sum) / (rho * cv);
                }
            }
            else // Flitti A, Pancheshnyi S 2009 model
            {
                // Plasma energy (carried by electrons) deposition rate
                // Unit: eV/s * electron discharge[C] * electron number density[#/m^3]
                // = (eV m^3/s * #/kmol * kmol/m^3) * e * (X_e * #/kmol * kmol/m^3) = J/s/m^3
                // Equivalent to: P_ext = e*N_e*E*U_e
                double P_ext = BoltzmannRate::bsolver.elec_power() * Avogadro * m_gas->molarDensity()* ElectronCharge
                               * m_gas->moleFraction(electronIndex) * Avogadro * m_gas->molarDensity();
                double P_elec = 1.5 * GasConstant * BoltzmannRate::bsolver.get_Te() * m_wdot[electronIndex];

                *eVibN2_dot = 0.0;
                *Ep_dot = P_ext;

                if (constPressure) {
                    *HE_dot = -hdot_vol + P_ext - P_elec;
                    *dTdt = (-hdot_vol + P_ext - P_elec) / (rho * cp);
                } else {
                    *HE_dot = -udot_vol + P_ext - P_elec;
                    *dTdt = (-udot_vol + P_ext - P_elec) / (rho * cv);
                }
            }
        }

        /* --------------------- SPECIES CONSERVATION EQS --------------------- */
        // the rate of change of each species' mass fraction is found using the closed-system
        // species conservation equation, applied once for each species:
        //     m*dY(k)/dt = dm(k)/dt
        // or equivalently:
        //     dY(k)/dt = dw(k)/dt * MW(k) / rho
        for (size_t k = 0; k < m_nSpecies; k++) {
            dYdt[k] = m_wdot[k] * m_gas->molecularWeight(k) / rho; // dYdt[k] is equivalent to *(dYdt + k)
            if(imposedNe) {
                dYdt[electronIndex] = 0.0;
            }
        }
    }

    // Judgment function on density change of CppBOLOS-configured species
    bool updateBoltzmannMixture() {
        const double RELATIVE_TOLERANCE = 0.1;
        const double ABSOLUTE_TOLERANCE = 0.01;
        const double WARNING_THRESHOLD = 1e-2;
        // TODO: adjust these values, check sensitivity
        bool change = false;

        for (auto& [sp, dens] : m_BoltzmannSpecies) {
            double moleFrac = m_gas->moleFraction(bSpeciesIndices.at(sp));
            double diff = abs(dens - moleFrac);
            if (diff > std::max(ABSOLUTE_TOLERANCE, dens * RELATIVE_TOLERANCE)) {
                change = true;
            }
        }

        // Check any species not configured by CppBOLOS exceeding WARNING_THRESHOLD
        for (const auto& [speciesName, index] : allSpeciesIndices) {
            double moleFrac = m_gas->moleFraction(index);
            if (moleFrac > WARNING_THRESHOLD) {
                std::cout << "[Warning]: Species '" << speciesName << "' not configured by CppBOLOS with mole fraction "
                          << moleFrac << " > " << WARNING_THRESHOLD << ".\n";
            }
        }

        if (change) {
            for (auto& [sp, dens] : m_BoltzmannSpecies) {
                dens = m_gas->moleFraction(bSpeciesIndices.at(sp)); // Update the density in BoltzmannSpecies map
            }
            std::cout << "CppBOLOS Configured Mixture Updated." << "\n";
        }
        return change;
    }

    /**
     * Number of equations in the ODE system.
     *   - overridden from FuncEval, called by the integrator during initialization.
     */
    size_t neq() const override {
        return m_nEqs;
    }

    /**
     * Provide the current values of the state vector, *y*.
     *   - overridden from FuncEval, called by the integrator during initialization.
     * @param[out] y solution vector, length neq()
     */
    void getState(double* y) override {
        // the solution vector *y* is [T, he, Ep, e_vib(N2), Y_1, Y_2, ... Y_K],
        // T: system temperature;
        // he: enthalpy or internal energy
        // Ep: plasma energy;
        // e_vib(N2): vibrational energy stored in N2
        // Y_k: mass fraction of species k.
        y[0] = m_gas->temperature();
        y[1] = m_he;
        y[2] = m_Ep;
        y[3] = eVibN2;
        m_gas->getMassFractions(&y[4]);
        double sumY = 0.0;
        for (size_t i = 4; i < m_nSpecies + 4; i++) {
            if (y[i] < 0.0) {
                std::cout << "Warning: Negative mass fraction for species " << m_gas->speciesName(i-4)
                << " (" << m_gas->massFraction(i-4) << ") " << "corrected to 0." << std::endl;
                y[i] = 0.0;
            }
            if (i != inertSpIndex + 4) {
                sumY += y[i];
            }
        }

        // Adjust the mass fraction of the idle species to ensure sum(Yi) = 1
        y[inertSpIndex + 4] = 1.0 - sumY;
    }

    void updateState(double* y)
    {
        checkFinite("y", y, m_nEqs);

        m_gas->setMassFractions_NoNorm(y+4);
        if (constPressure) {
            m_gas->setState_TP(y[0], m_pressure);
        } else {
            m_gas->setState_TD(y[0], m_density);
        }

        m_he = y[1];
        m_Ep = y[2];
        eVibN2 = y[3];
    }

    void setConstPD (double pressure, double density) {
        m_pressure = pressure;
        m_density = density;
    }

    /* ------------------- Functions to check production rates -------------------*/
    // Maximum production rate
    size_t w_dot_max_index() {
        auto max_it = std::max_element(m_wdot.begin(), m_wdot.end());
        return std::distance(m_wdot.begin(), max_it);
    }

    double w_dot_max() {
        return m_wdot[w_dot_max_index()];
    }

    std::string w_dot_max_name() {
        return m_gas->speciesName(w_dot_max_index());
    }

    // Electron production rate [kmol/s/m^3]
    double w_dot_e() {
        return m_wdot[electronIndex];
    }

    void printTopSpeciesProductionRates() {
        m_kinetics->getNetProductionRates(&m_wdot[0]); // kmol/m^3/s or kmol/m^2/s
        // Vector to store (index, rate of progress) pairs
        vector<std::pair<size_t, double>> rate_values;
        for (size_t i = 0; i < m_wdot.size(); ++i) {
            // Store the absolute value of the rate of progress
            rate_values.emplace_back(i, std::abs(m_wdot[i]));
        }

        // Sort by the rate of progress value, descending order
        std::sort(rate_values.begin(), rate_values.end(),
                  [](const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) {
                      return a.second > b.second; // Sort in descending order
                  });

        // Print top 5 species
        for (size_t j = 0; j < std::min(5, static_cast<int>(rate_values.size())); ++j) {
            auto index = rate_values[j].first;
            writelog("Sp({:1d})  {:30s} {:14.5g}   kmol/m3/s\n",
                     index, m_gas->speciesName(index), m_wdot[index]);
        }

    }

    void printTopRatesOfProgress() {
        m_kinetics->getNetRatesOfProgress(&m_q_net[0]);

        // Vector to store (index, rate of progress) pairs
        vector<std::pair<size_t, double>> rate_values;

        for (size_t i = 0; i < m_q_net.size(); ++i) {
            // Store the absolute value of the rate of progress
            rate_values.emplace_back(i, std::abs(m_q_net[i]));
        }

        // Sort by the rate of progress value, descending order
        std::sort(rate_values.begin(), rate_values.end(),
                  [](const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) {
                      return a.second > b.second; // Sort in descending order
                  });

        // Print top 5 reactions
        for (size_t j = 0; j < std::min(5, static_cast<int>(rate_values.size())); ++j) {
            auto index = rate_values[j].first;
            writelog("R({:1d})  {:30s} {:14.5g}   kmol/m3/s\n",
                     index, m_kinetics->reaction(index)->equation(), m_q_net[index]);
        }
    }

    /* ------------------- Functions to print fast gas heating info -------------------*/

    void printMaxFGH(){
        double maxFGH = 0.0;
        size_t  maxFGH_index;
        vector<double> deltaH(m_kinetics->nReactions());
        m_kinetics->getNetRatesOfProgress(&m_q_net[0]);
        m_kinetics->getDeltaEnthalpy(&deltaH[0]);
        for (const auto& [i, e_th] : e_th_map) {
            const auto& rxn = m_kinetics->reaction(i);
            // Unit: J/m^3/s (= eV * #/kmol * J/eV - J/kmol) * kmol/m^3/s
            double FGH = ( e_th_map[i] * Avogadro * ElectronCharge - deltaH[i] ) * m_q_net[i];
            if (FGH > maxFGH) {
                maxFGH = FGH;
                maxFGH_index = i;
            }

        }
        writelog("R({:1d})  {:30s} {:14.5g}   kmol/m3/s\n",
                 maxFGH_index, m_kinetics->reaction(maxFGH_index)->equation(), m_q_net[maxFGH_index]);
    }

    void printTopFGH(){
        m_kinetics->getNetRatesOfProgress(&m_q_net[0]);
        m_kinetics->getDeltaEnthalpy(&m_deltaH[0]);

        // Vector to store (index, FGH value) pairs
        vector<std::pair<size_t, double>> FGH_values;

        for (const auto& [i, e_th] : e_th_map) {
            // Calculate FGH
            double FGH = (e_th * Avogadro * ElectronCharge - m_deltaH[i]) * m_q_net[i];
//            double FGH = (e_th * Avogadro * ElectronCharge) * q_net[i];
            FGH = m_q_net[i];
            FGH_values.emplace_back(i, FGH);
        }

        // Sort by FGH value, descending order
        std::sort(FGH_values.begin(), FGH_values.end(),
                  [](const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) {
                      return a.second > b.second; // Sort in descending order
                  });

        // Print top 5 reactions
        for (size_t j = 0; j < std::min(10, static_cast<int>(FGH_values.size())); ++j) {
            auto maxFGH_index = FGH_values[j].first;
            writelog("R({:1d})  {:30s} {:14.5g}  kmol/m3/s {:14.5g}  J/m^3/s\n",
                     maxFGH_index, m_kinetics->reaction(maxFGH_index)->equation(), m_q_net[maxFGH_index], FGH_values[j].second);
        }
    }

    // Use external electron number density [#/cm^3] profile
    void imposeNe(double ne) {
        m_N_e = ne;
        imposedNe = true;
    }

    double depositedPlasmaEnergy() const {
        return m_Ep;
    }

    void resetDepositedPlasmaEnergy() {
        m_Ep = 0.0;
    }

    bool constPressure = true;

    bool nonThermal = false;

    bool heatLoss = false;

    double T0, C0; // Ambient temperature

    size_t electronIndex = -1;
    size_t inertSpIndex = -1;
    size_t nSteps = 0;

    size_t findSpeciesIndex(const std::vector<std::string>& legalNames) const {
        for (const auto& name : legalNames) {
            size_t index = m_gas->speciesIndex(name);
            if (index < std::numeric_limits<size_t>::max()) {
                return index;
            }
        }
        std::cout << "ChemPlasReactor::findSpeciesIndex ERROR: Registering species index of possible names: [";
        for (const auto& name : legalNames) {
            std::cout << name << " ";
        }
        throw std::runtime_error("]\nNo valid species name found! Consider add as dummy species.");
    }

    size_t findSpeciesIndex(const std::string name) {
        size_t index = m_gas->speciesIndex(name);
        if (index < std::numeric_limits<size_t>::max()) {
            return index;
        }
        std::cout << "ChemPlasReactor::findSpeciesIndex ERROR: Register species index of: " << name << std::endl;
        throw std::runtime_error("No valid species name found! Consider add as dummy species.");
    }

private:
    // Private member variables, to be used internally.
    shared_ptr<ThermoPhase> m_gas;
    shared_ptr<Kinetics> m_kinetics;
    vector<double> m_hbar;
    vector<double> m_ubar; // internal energy
    vector<double> m_wdot; // net production rate of species, kmol/m^3/s, Length: nSpecies
    vector<double> m_q_net; // net production rate of individual reactions, kmol/m^3/s, Length: nReactions()
    vector<double> m_deltaH;
    double m_pressure;
    size_t m_nSpecies;
    size_t m_nEqs;
    double m_density; // fixed density for constant volume case
    double m_Ep = 0.0; // Deposited energy. (J/m^3)
    double m_he = 0.0;

    std::map<std::string, double>& m_BoltzmannSpecies;
    std::map<std::string, size_t> allSpeciesIndices;
    std::map<std::string, size_t> bSpeciesIndices; // indices of species configured by CppBOLOS

    double m_N_e = 0.0; // Number density of e-
    bool imposedNe = false;

    bool detailPlasmaHeatModel = true;

    // Maps to store reaction number and corresponding energy_transfer values
    std::map<size_t, double> e_th_map; // energy transfer in the reaction steps, can be negative
    std::map<size_t, double> e_ext_map; // Externally deposited plasma energy Boltzmann reactions, must positive
    // Plasma energy stored in vibration states of species: N2, O2, NO, NH3, etc
    std::map<size_t, double> evib_N2_map;

    double eVibN2 = 0.0;

    // Register plasma energy transfer data
    void setPlasmaEnergyTransfer(const vector<AnyMap>& reactions) {
        for(size_t i = 0; i < reactions.size(); i++) {
            auto& reaction = reactions[i];
            if(reaction.hasKey("energy_transfer")) {
//            const YAML::Node& energy_transfer = reactions[i]["energy_transfer"];
                auto& energy_transfer = reaction["energy_transfer"].as<AnyMap>();

                // Check for e_th
                if(energy_transfer.hasKey("e_th")) {
                    double e_th = parseEnergyString(energy_transfer["e_th"].as<std::string>());
                    e_th_map[i] = e_th;
                    if (reaction.hasKey("type") && reaction["type"] == "Boltzmann") {
                        e_ext_map[i] = e_th;
                    }
                }

                // Check for evib_N2
                if(energy_transfer.hasKey("evib_N2")) {
                    double evib_N2 = parseEnergyString(energy_transfer["evib_N2"].as<std::string>());
                    evib_N2_map[i] = evib_N2;
                    if (reaction.hasKey("type") && reaction["type"] == "Boltzmann") {
                        e_ext_map[i] = evib_N2;
                    }
                }
            }
        }
    }

    // Utility function to parse custom energy format
    static double parseEnergyString(const std::string& energyString) {
        // Find the position of "_eV"
        size_t pos = energyString.find("_eV");
        if(pos != std::string::npos) {
            // Convert the part before "_eV" to double
            return std::stod(energyString.substr(0, pos));
        } else {
            // Throw an exception or return a sentinel value, based on your preference
            throw std::runtime_error("Invalid energy format: " + energyString);
        }
    }

};

}

#endif //CHEMPLASKIN_PLASMAREACTOR_H
