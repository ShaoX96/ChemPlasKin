//
// Created by Xiao Shao on 2023/9/17.
//

#include "cantera/kinetics/PlasmaCustomExpr.h"
#include "cantera/thermo/ThermoPhase.h"

namespace Cantera {

PlasmaCustomExpr::PlasmaCustomExpr(const AnyMap &node, const UnitStack &rate_units) {
    ReactionRate::setParameters(node, rate_units);
    if (node.hasKey("rateExpr")) {
        setRateParameters(node["rateExpr"], node.units(), rate_units);
    } else{
        throw InputFileError("PlasmaCustomExpr", m_input, "\"rateExpr\" is missing!");
    }

    // Check for presence of each variable and dynamically register it
    /*
    for (const auto& var : variables) {
        if (customExpr.find(var.first) != std::string::npos) {
            parser.DefineVar(var.first, const_cast<double*>(&var.second));
        }
    } */

    // Dynamic registration adds complexity. Use fixed variable registration
    parser.DefineVar("Tgas", &Tgas);
    parser.DefineVar("Te", &Te);
    parser.DefineVar("EN", &EN);

    try {
        // Parse the expression
        LOG_DEBUG( "mu::Parser reading custom rateExpr: " << rateExpr );
        parser.SetExpr(rateExpr);
        LOG_DEBUG("Testing parser.Eval() ..." << parser.Eval() << ", passed.");
    } catch (const mu::ParserError& e) {
        std::cerr << "rateExpr '" << rateExpr << "'" << std::endl;
        std::cerr << "mu::ParserError: " << e.GetMsg() << std::endl;

        if (e.GetCode() == mu::ecUNASSIGNABLE_TOKEN) {
            std::cerr << "Unknown variable detected. Please check your expression." << std::endl;
        }
        // std::cerr.flush(); // Ensure the error stream is flushed immediately
        // Force code termination immediately if there is mu::Parser encounter an error.
        std::exit(EXIT_FAILURE);
    }
}

void PlasmaCustomExpr::setRateParameters(
        const AnyValue& rate, const UnitSystem& units, const UnitStack& rate_units)
{
    if (rate.is<AnyMap>()) {
        auto& rate_map = rate.as<AnyMap>();
        m_A = units.convertRateCoeff(rate_map[m_A_str], conversionUnits());
        if (rate_map.hasKey("Expr")){
            rateExpr = rate_map["Expr"].asString();
        }
    } else {
        throw InputFileError("PlasmaCustomExpr", m_input, "Incorrect \"rateExpr\" format!");
    }
}

bool PlasmaCustomData::update(const ThermoPhase& phase, const Kinetics& kin)
{
    double T = phase.temperature();
    // If gas or electron temeprature does not change, no update
    // We may reduce update frequency by setting a tolerence of change
    if (std::abs(T - temperature) < 0.001 && temp_Te == electronTemp) {
        return false;
    }
    update(T);
    temp_Te = electronTemp;

    return true;
}

}