//
// Created by Xiao Shao on 2023/9/5.
//
#include "cantera/kinetics/Boltzmann.h"
#include "cantera/thermo/ThermoPhase.h"

namespace Cantera
{

CppBOLOS::BoltzmannSolver BoltzmannRate::bsolver;
Eigen::VectorXd BoltzmannRate::F0;
double BoltzmannData::reducedEfield = 0.0;
double BoltzmannData::electronTemp = 1.0;
size_t BoltzmannRate::NumProcess = 0;
bool BoltzmannRate::updateOnce = true;

bool BoltzmannRate::useCachedSolutions = false;
std::vector<std::string> BoltzmannRate::ProcessesList;
std::vector<double> BoltzmannRate::cachedSolutions;

// Update only once for each writing time step.
bool BoltzmannData::update(const ThermoPhase& phase, const Kinetics& kin)
{
    if(BoltzmannRate::updateOnce) {
        BoltzmannRate::updateOnce = false;
//        std::cout << "Update yes!" << std::endl;
        return true;
    } else {
//        std::cout << "Update No!" << std::endl;
        return false;
    }

}


BoltzmannRate::BoltzmannRate(const AnyMap& node, const UnitStack& rate_units)
{
    // Ensure the Boltzmann solver has been set up
    if (F0.size() == 0){
        std::cerr << "Error: BoltzmannRate::bsolver must be initialized before creating BoltzmannRate objects." << std::endl;
    }

    if (node.hasKey("process")) {
        process = node["process"].asString();
    } else{
        throw InputFileError("BoltzmannRate", m_input, "\"process\" is missing!");
    }

    LOG_INFO("Initial BoltzmannRate of " << process <<": " << bsolver.rate(F0, process) );
    ProcessesList.push_back(process);
    processIndex = NumProcess;
    NumProcess += 1;
}

// ... Implement other methods

}