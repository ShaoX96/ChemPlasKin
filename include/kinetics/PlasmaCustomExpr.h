//
// Created by Xiao Shao on 2023/9/17.
//

#ifndef PLASMACUSTOMEXPR_H
#define PLASMACUSTOMEXPR_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/Units.h"
#include "cantera/kinetics/ReactionData.h"
#include "ReactionRate.h"
#include "MultiRate.h"
#include "Boltzmann.h"
#include "cantera/ext/muparser/muParser.h"

namespace Cantera {

class AnyValue;
class AnyMap;

struct PlasmaCustomData : public BoltzmannData
{
    bool update(const ThermoPhase& phase, const Kinetics& kin) override;
    using ReactionData::update;

    // Electron temperature used for update judgement
    double temp_Te = 1.0;
};

//! Custom expression of plasma reaction rate. Flexible but may not be efficient
class PlasmaCustomExpr : public ReactionRate
{
public:

    PlasmaCustomExpr(){}

    explicit PlasmaCustomExpr(const AnyMap& node, const UnitStack& rate_units={});

    // In case the constructor is copied, the variables shall be redefined
    PlasmaCustomExpr(const PlasmaCustomExpr& other)
            : parser(other.parser), Tgas(other.Tgas), Te(other.Te), EN(other.EN), m_A(other.m_A) // and so on for other members
    {
        // std::cout << "Copy constructor called!" << std::endl;
        // Now that the members have been copied, re-register the variables with the parser.
        parser.DefineVar("Tgas", &Tgas);
        parser.DefineVar("Te", &Te);
        parser.DefineVar("EN", &EN);
        // Do the same for any other necessary variables.

//        if (rateExpr.find("Tgas") == std::string::npos) {
//            updateTgas = false;
//        }
//        if (rateExpr.find("Te") == std::string::npos) {
//            updateTe = false;
//        }

    }

    void setRateParameters(const AnyValue& rate,
                           const UnitSystem& units,
                           const UnitStack& rate_units);

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return make_unique<MultiRate<PlasmaCustomExpr, PlasmaCustomData>>();
    }

    double evalFromStruct(const PlasmaCustomData& shared_data) {
        Tgas = shared_data.temperature;
        Te = shared_data.electronTemp;
        EN = shared_data.reducedEfield;
        return m_A*parser.Eval();
    }

    const string type() const override {
        return "PlasmaCustomExpr";
    }

private:
    mu::Parser parser;
    std::string rateExpr;
    // Tgas, Te, and EN are defined only for mu::Parser evaluation
    double Tgas = 300.0;
    double Te = 1.0;
    double EN = 0.0;
    double m_A = NAN; //!< Pre-exponential factor
    string m_A_str = "A"; //!< The string for the pre-exponential factor

//    bool updateTgas = true;
//    bool updateTe = true;


// We give up this map method to avoid re-registration of parser varibles.
// std::map<std::string, double> variables = {{"Tgas", 300.0}, {"Te", 1.0}};
};

}

#endif //PLASMACUSTOMEXPR_H
