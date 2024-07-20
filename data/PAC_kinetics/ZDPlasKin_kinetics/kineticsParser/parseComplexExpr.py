#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 24 17:25:55 2023

@author: shaox

Parse ZDPlasKin reactions with complex expressions
"""
import re
import sympy as sp

def parse_variables(fortran_code):
    translated_variables = {}

    # Split the Fortran code into lines
    lines = fortran_code.strip().split('\n')
    
    # Define the regex pattern for matching variable definitions
    pattern_param = re.compile(r"\$\s*([\w\s,]+)::\s*([\w\s=.,\d*+-dD]+)(?:\s*!.*|$)")
    pattern = re.compile(r"\$\s+([a-zA-Z_0-9]+)\s*=\s*([^!]+)(?:\s*!.*|$)")

    # Strip out definitions with initilization
    for line in lines:
        match = pattern_param.match(line)
        if match:
            initializations = match.group(2)
            print(initializations, '\n')
            init_expressions = ['$ ' + expr.strip() for expr in initializations.split(',')]
            lines = init_expressions + lines
    
    # Extract variable and definitions                
    for line in lines:       
        match = pattern.match(line)
        if match:
            # Extract the variable name and expression
            var_name, expression = match.groups()
            
            # Translate the expression to Python syntax
            python_expression = translate_expression(expression)
            
            # Store the translated expression
            translated_variables[var_name] = python_expression
    
    return translated_variables

def translate_expression(expression):
    # Replace Fortran-specific syntax with Python syntax
    python_expression = expression.replace('d', 'e').strip()
    # python_expression = python_expression.replace('**', '^')
    # ... any other translation rules
    return python_expression


def substitute_and_simplify(expression, variables):
    # First, make the expression a sympy object
    expr_sympy = sp.sympify(expression)
    
    # Now substitute each predefined variable in
    for var, expr in variables.items():
        expr_sympy = expr_sympy.subs(var, sp.sympify(expr))
    
    # Now attempt to simplify the expression
    expr_simplified = sp.simplify(expr_sympy)
    
    # Further attempt to rationalize and simplify radicals
    # expr_simplified = sp.nsimplify(expr_simplified)
    # expr_simplified = sp.radsimp(expr_simplified)
    
    return expr_simplified


def preprocess_variables(parsed_variables):
    changes_made = True  # Flag to check if any changes were made during a pass
    
    while changes_made:
        changes_made = False
        
        # Make a copy of the dictionary to compare with after the pass
        original_vars = parsed_variables.copy()
        
        for var, expr in parsed_variables.items():
            simplified_expr = substitute_and_simplify(expr, parsed_variables)
            parsed_variables[var] = str(simplified_expr)  # Update with simplified expression
        
        # Check if any variable's expression has changed during this pass
        for var in parsed_variables:
            if parsed_variables[var] != original_vars[var]:
                changes_made = True
                break

    return parsed_variables



if __name__ == "__main__":

    
    # Example usage:
    fortran_code = """
$ double precision :: QvibH2, kVT10_H2O2, kVT10_H2H2O, kVT10_H2H2, kVT10_H2OH, kVT10_H2H, kVT10_H2O, kVT10_H2HE ! cm3.s-1
$ double precision :: QvibO2, kVT10_O2O2, kVT10_O2H2O, kVT10_O2H2, kVT10_O2OH, kVT10_O2H, kVT10_O2O, kVT10_O2HE ! cm3.s-1
$ double precision :: kVV10_H2H2, kVV10_O2O2 ! cm3.s-1
$ double precision :: M ! cm-3
$ double precision :: ktot, k0, kR1, kR2, kR3 ! cm3.s-1
$ double precision, parameter :: energy_vibH2= 0.516d0*11605.0d0, energy_vibO2=0.190d0*11605.0d0 ! K
$ QvibH2 = exp(-energy_vibH2/Tgas)
$ QvibO2 = exp(-energy_vibO2/Tgas)

# Rate constant is calculated at 1 atm.
$ M = 101325.0d0/1.38064852d-17/Tgas

$ kVT10_H2O2 = 0.08*exp(-19.94d0-144.9d0/Tgas**(1.0d0/3.0d0))
$ kVT10_H2H2O = 0.23*exp(-19.94d0-144.9d0/Tgas**(1.0d0/3.0d0)) 
$ kVT10_H2H2 = exp(-19.94d0-144.9d0/Tgas**(1.0d0/3.0d0))
$ kVT10_H2OH = 0.23*exp(-19.94d0-144.9d0/Tgas**(1.0d0/3.0d0)) 
$ kVT10_H2H = 1.73d-14/Tgas**0.5
$ kVT10_H2O = 2.2d-12
$ kVT10_H2HE = 1.31d-13*Tgas*exp(-95.2d0/Tgas**(1.0d0/3.0d0))/(1.0d0-QvibH2)

$ kVT10_O2O2 = exp(-18.32d0-157.0d0/Tgas**(1.0d0/3.0d0)) 
$ kVT10_O2H2O = 100*exp(-18.32d0-157.0d0/Tgas**(1.0d0/3.0d0)) 
$ kVT10_O2H2 = exp(-22.04d0-91.5d0/Tgas**(1.0d0/3.0d0))
$ kVT10_O2OH = 100*exp(-18.32d0-157.0d0/Tgas**(1.0d0/3.0d0)) 
$ kVT10_O2H = 1.7*exp(-22.04d0-91.5d0/Tgas**(1.0d0/3.0d0))
$ kVT10_O2O = Tgas*exp(-28.87d0-111.0d0/Tgas**(1.0d0/3.0d0)) 
$ kVT10_O2HE = 4.54d-14*Tgas*exp(-60.85d0/Tgas**(1.0d0/3.0d0))/(1.0d0-QvibO2)

$ kVV10_H2H2 = 4.23d-15*(300.0d0/Tgas)**(1.0d0/3.0d0)
$ kVV10_O2O2 = 9.0d-15*(Tgas/300)**1.5

$ ktot = 3.0d-11*Tgas**0.17*exp(-2726.0d0/Tgas)
$ k0 = 2.2d-37*Tgas**2.34*exp(-1350.0d0/Tgas)
$ kR1 = 6.8d-16*Tgas**1.51*exp(-3040.0d0/Tgas)
$ kR2 = ktot*(k0*M/(ktot+k0*M))
$ kR3 = ktot-kR1-kR2
    """
    
    parsed_variables = parse_variables(fortran_code)
    parsed_variables = preprocess_variables(parsed_variables)

    print(parsed_variables, "\n")
    
    reaction_expression = " QvibO2 *kVT10_H2H2O*2.06*exp(14.92e0*(1.882e0/Tgas)**0.5) "
    simplified_expression = substitute_and_simplify(reaction_expression, parsed_variables)
    print(simplified_expression)

