//
// Created by Xiao Shao on 2023/9/5.
//

#ifndef CT_BOLTZMANN_H
#define CT_BOLTZMANN_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/Units.h"
#include "cantera/kinetics/ReactionData.h"
#include "ReactionRate.h"
#include "MultiRate.h"

#include "cantera/ext/bolos/solver.h"  // include CppBOLOS's solver
#include "cantera/ext/bolos/Logger.h"

namespace Cantera {

class AnyValue;
class AnyMap;

struct BoltzmannData : public ReactionData
{
    BoltzmannData() = default;

    bool update(const ThermoPhase& phase, const Kinetics& kin) override;
    using ReactionData::update;

    static double electronTemp; //!< electron temperature
    static double reducedEfield; //!< reduced electric field [Td]

};

class BoltzmannRate : public ReactionRate
{
public:
//
    BoltzmannRate(const AnyMap& node, const UnitStack& rate_units={});

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return make_unique<MultiRate<BoltzmannRate, BoltzmannData>>();
    }

    double evalFromStruct(const BoltzmannData& shared_data) const {
        LOG_DEBUG( "Calculating evalFromStruct of BoltzmannRate of " << process << ": " << Avogadro*bsolver.rate(F0, process) );
        if (useCachedSolutions) {
            return cachedSolutions[processIndex];
        }
        return Avogadro * bsolver.rate(F0, process);
//        std::cout << "Debug line 1. BoltzmannRate " << "process" << process << std::endl;
        // For each writing time step, only recalculate reaction rate of each process ONCE.
        // Avoid recalculating in each integration dt (CVODE)
        /*
        if(countEval < NumProcess){
            countEval += 1;
            rateValue = bsolver.rate(f1,process); // cache the result
//            std::cout << "Debug line 2. Updating BoltzmannRate. rateValue: " << rateValue << std::endl;
            return Avogadro*rateValue; // convert unit to m^3/(kmol*s)
        } else{
//            std::cout << "Debug line 2. NoUpdating BoltzmannRate. rateValue: " << rateValue << std::endl;
            return Avogadro*rateValue;
        } */
    }

    const string type() const override {
        return "Boltzmann";
    }

    static CppBOLOS::BoltzmannSolver bsolver;

    static void updateBoltzmannSolver(const int maxItr=200, const double rtol=1e-5, const double delta0=1E14){
        // bsolver.init();
        F0 = bsolver.converge(F0, maxItr, rtol, delta0); // converge() will update Te
//        BoltzmannData::electronTemp = bsolver.get_Te();
        BoltzmannData::electronTemp = std::max(bsolver.get_Te(), bsolver.get_kT()); // enforce T_e >= T_gas
        BoltzmannData::reducedEfield = bsolver.get_EN();
        updateOnce = true;
    }

    static size_t NumProcess;
    static bool updateOnce;

    static Eigen::VectorXd F0; // EEDF

    static bool useCachedSolutions;
    static std::vector<std::string> ProcessesList;
    static std::vector<double> cachedSolutions;

    static void cacheBoltzmannSolutions(const int maxItr=200, const double rtol=1e-5, const double delta0=1E14)
    {
        updateBoltzmannSolver(maxItr, rtol, delta0);

        // Create solution vector [k1, k2, ..., k_N, Te, E/N]
        cachedSolutions.resize(NumProcess + 2);
        size_t i = 0;
        while (i < ProcessesList.size()) {
            cachedSolutions[i] = Avogadro * bsolver.rate(F0, ProcessesList[i]);
            i ++;
        }
        cachedSolutions[i] = BoltzmannData::electronTemp;
        cachedSolutions[i+1] = BoltzmannData::reducedEfield;
    }

private:
    std::string process;
    size_t processIndex;

};

}
#endif //CT_BOLTZMANN_H

