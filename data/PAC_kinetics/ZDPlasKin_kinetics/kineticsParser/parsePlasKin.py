"""
>>>>>>>>>>>><<< Part of ChemPlasKin Project <<<<<<<<<<<<<<<<<<<

A python parser translating plasma kinetics input in ZDPlasKin format to YAML format. 
Utilized to build a unified chemistry-plasma kinetic mechanism input.

Inline definition of variables will also be parsed, substituted and simplified.

Written by Xiao Shao @KAUST, 2023
"""
import re
import yaml
from collections import Counter
import sympy as sp
import argparse
import logging


def read_reactions(filename=''):
    """
    Read REACTION section lines with expanded reactions based on group syntax. 
    
    Returns
    -------
    variables: list of lines containing inline definition of variables started with "$"
    reactions: list lines containing reaction steps expanded based on group syntax
    """
    with open(filename, 'r') as file:
        lines = file.readlines()

    in_species_section = False
    in_reactions_section = False
    species = []
    variables = []
    reactions = []
    
    def expand_reaction(line, groups):
        """Helper function to expand a reaction line based on group syntax."""
        # Extract all group symbols from the line
        group_symbols = [symbol for symbol in groups.keys() if symbol in line]
    
        # Total number of substitutions for each symbol
        n_substitutions = len(groups[group_symbols[0]])
    
        # Create reactions for each combination of substitutions
        expanded_lines = []
        for i in range(n_substitutions):
            substituted_line = line
            for symbol in group_symbols:
                # Replace all occurrences of the symbol
                while symbol in substituted_line:
                    substituted_line = substituted_line.replace(symbol, groups[symbol][i], 1)
            expanded_lines.append(substituted_line)
            
        logging.info("Expanding '{}' :\n{}".format(line, '\n'.join(expanded_lines)))
    
        return expanded_lines

    for i in range(len(lines)):
        line = lines[i].strip()
        
        if line == "SPECIES":
            in_species_section = True     
        elif line == "REACTIONS":
            in_reactions_section = True
        elif line == "END":
            in_species_section = False
            in_reactions_section = False
        elif in_species_section and line and not line.startswith('#') and not line.startswith('!'):
            species.extend(line.split())
        elif in_reactions_section and line and not line.startswith('#') and not line.startswith('!'):
            if line.startswith('$'):
                variables.append(line)
            elif "=>" in line and "@" in line:
                # Gather group substitutions
                groups = {}
                j = i + 1  # Start from the next line
                group_pattern = re.compile(r'^\s*@\S+\s*=.+')  # in format of " @A = O  O2 ..."
                while group_pattern.match(lines[j]):
                    group_symbol, substitutions = lines[j].strip().split("=")
                    groups[group_symbol.strip()] = [s.strip() for s in substitutions.split()]
                    j += 1
                # Now expand the reaction
                expanded_lines = expand_reaction(line, groups)
                reactions.extend(expanded_lines)
                i = j  # Jump to line after the last group definition
            elif "=>" in line and "!" in line:
                reactions.append(line)

    return species, variables, reactions



def parse_variables(lines):
    """
    Parse all variables defined as inline Fortran code
    
    Returns
    -------
    A dictionary of variables: (name: definition)
    """
    translated_variables = {}
    
    # Define the regex pattern for matching variable definitions
    pattern_param = re.compile(r"\$\s*([\w\s,]+)::\s*([\w\s=.,\d*+-dD]+)(?:\s*!.*|$)")
    pattern = re.compile(r"\$\s+([a-zA-Z_0-9]+)\s*=\s*([^!]+)(?:\s*!.*|$)")

    # Strip out definitions after '::', which may hold initial values
    for line in lines:
        match = pattern_param.match(line)
        if match:
            initializations = match.group(2)
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


def preprocess_variables(parsed_variables):
    """
    Ensure all the variable expressions in their simplest forms 
    for more effective substitution into subsequent expressions
    """
    changes_made = True  # Flag to check if any changes were made during a pass
    
    while changes_made:
        changes_made = False
        
        # Make a copy of the dictionary to compare with after the pass
        original_vars = parsed_variables.copy()
        
        for var, expr in parsed_variables.items():
            simplified_expr, _ = substitute_and_simplify(expr, parsed_variables)
            parsed_variables[var] = str(simplified_expr)  # Update with simplified expression
        
        # Check if any variable's expression has changed during this pass
        for var in parsed_variables:
            if parsed_variables[var] != original_vars[var]:
                changes_made = True
                break

    return parsed_variables


def translate_expression(expression):
    """
    Replace Fortran-specific syntax with Python syntax
    """
    python_expression = expression.replace('d', 'e').strip()
    # python_expression = python_expression.replace('**', '^')
    # ... any other translation rules
    return python_expression


def substitute_and_simplify(expression, variables):
    """
    Generate simplified expression by substitution
    Returns
    -------
    expr_simplified : string

    """
    # First, make the expression a sympy object
    expr_sympy = sp.sympify(expression)
    
    # Store the initial expression for comparison later
    initial_expr = expr_sympy
    
    # Now substitute each predefined variable in
    for var, expr in variables.items():
        expr_sympy = expr_sympy.subs(var, sp.sympify(expr))
    
    # Now attempt to simplify the expression
    expr_simplified = sp.simplify(expr_sympy)
    
    # Further attempt to rationalize and simplify radicals
    # expr_simplified = sp.nsimplify(expr_simplified)
    # expr_simplified = sp.radsimp(expr_simplified)
    
    changed = (expr_simplified != initial_expr)
    
    return expr_simplified, changed
    

def parse_reactions(variables, reactions):
    """
    Parse reactions lines. Each reaction is registered as a dictionary with
    equation, type, and rate expression.
    Parameters
    ----------
    variables: list
        dictionary of parsed variables
    reactions : list
        lines of reaction steps.
    Returns
    -------
    parsed_reactions : list
        a list parsed reactions with reaction_data dictionary as elements.

    """
    parsed_reactions = []

    for reaction in reactions:
        reaction_data = {}
        
        # if '!' in reaction:
        equation, comment = re.split(r'\s*!\s*', reaction, maxsplit=1)
        
        # Re-format reaction equations
        equation = equation.replace('=>', ' => ')
        equation = re.sub(r'(?<!\^)\+', ' + ', equation) # insert space around '+', except '^+'
        equation = re.sub(r'\s+', ' ', equation).strip() # replace multiple spaces with one
        if 'ANY_NEUTRAL' in equation:
            equation = equation.replace('ANY_NEUTRAL', 'M')
            reaction_data['efficiencies'] = 0.0
            
        equation = re.sub(r'(?<=\s)(\d)(?=[A-Za-z])', r'\1 ', equation) # Separate digital starting species
        
        reaction_data['equation'] = equation
        
        # Determine the number of reactant species
        reactants, _ = equation.split('=>')
        reactant_count = len(reactants.split(' + '))

        # Determine the unit based on the reactant count
        if reactant_count == 1:
            unit = '1/s'
        elif reactant_count == 2:
            unit = 'cm^3/molec/s'
        elif reactant_count == 3:
            unit = 'cm^6/molec^2/s'
        else:
            raise ValueError(f'Unexpected number of reactant species: {reactant_count}')

        matchBOLSIG = re.search(r'BOLSIG\+?\s*(.*)', comment)
        if matchBOLSIG:
            reaction_data['type'] = 'Boltzmann'
            reaction_data['process'] = matchBOLSIG.group(1).strip()
        # if 'BOLSIG' in comment:
        #     reaction_data['type'] = 'Boltzmann'
        #     reaction_data['process'] = comment.split('BOLSIG')[1].strip()
        # elif re.match(r'\d+(\.\d+[deDE][+-]\d+)', comment):
        else:
            comment_old = re.sub('d', 'e', comment.strip())
            comment, changed = substitute_and_simplify(comment_old, variables)
            if changed:
                logging.info(f"Rephrasing '{comment_old}' to: {comment}")            
            if '*' in str(comment) or '/' in str(comment):
                reaction_data['type'] = 'PlasmaCustomExpr'
                expr = str(comment).replace('**', '^')
                reaction_data['rateExpr'] = { 'A': f'{1.0} {unit}', 'Expr': expr }
            else:
                try:
                    rate_constant = float(comment)  # can simplify output format
                except:
                    logging.info(f"\n !!!!!!!! Error in pharsing expression: '{comment}' !!!!!!!!")
                    raise SystemExit("Exiting the script due to parsing reaction rate failue.")
                reaction_data['rate-constant'] = {'A': f'{rate_constant} {unit}', 'b': 0.0, 'Ea': 0.0}

        parsed_reactions.append(reaction_data)
            
    # After parsing, count the occurrences of each equation
    equation_counts = Counter(reaction_data['equation'] for reaction_data in parsed_reactions)
    
    # Now mark duplicates in parsed_reactions
    for reaction_data in parsed_reactions:
        equation = reaction_data['equation']
        if equation_counts[equation] > 1:
            reaction_data['duplicate'] = "true"
            
    logging.info("\nReactions sucessfully parsed! \n")

    return parsed_reactions


# def write_yaml(parsed_reactions, filename='kinetics.yaml'):
#     with open(filename, 'w') as file:
#         yaml.dump({'reactions': parsed_reactions}, file, default_flow_style=None)
        
def write_yaml(charge_coefficient, parsed_reactions, markIndex = False, filename=''):
    """
    Write parsed reactions into a *yaml file. Avoid auto writing using:
        with open(filename, 'w') as file:
            yaml.dump({'reactions': parsed_reactions}, file, default_flow_style=None)
    to have a full control of output format.
    
    Parameters
    ----------
    charge_coefficient : string
        DESCRIPTION. Hard-coded charges efficiencies (all 0.0)
    parsed_reactions : list of dictionaries
        DESCRIPTION.
    filename : TYPE, optional
        DESCRIPTION. Output file. The default is 'kinetics.yaml'.
    Returns
    -------
    None.

    """
    with open(filename, 'w') as file:
        
        # file.write('reactions:\n')
        file.write('# /------------------ Plasma Kinetics Reactions ---------------------\n')
        logging.info(f"Writing reactions to {filename}:")
        i = 0
        for reaction_data in parsed_reactions:
            i += 1
            logging.info(f"R({i}): {reaction_data['equation']}")
            space = ' ' * (40 - len(reaction_data['equation']))
            file.write(f"- equation: {reaction_data['equation']}{space}{'# R(' + str(i) + ')' if markIndex else ''}\n")
            
            if 'duplicate' in reaction_data:
                file.write('  duplicate: ' + reaction_data['duplicate'] + '\n')
            
            if 'type' in reaction_data:
                file.write('  type: ' + reaction_data['type'] + '\n')
            
            if 'process' in reaction_data:
                file.write('  process: ' + reaction_data['process'] + '\n')
            
            if 'rate-constant' in reaction_data:
                rate_constant = reaction_data['rate-constant']
                file.write(f'  rate-constant: {{A: {rate_constant["A"]}, b: {rate_constant["b"]}, Ea: {rate_constant["Ea"]}}}\n')
            
            if 'rateExpr' in reaction_data:
                rateExpr = reaction_data['rateExpr']
                file.write(f'  rateExpr: {{A: {rateExpr["A"]}, Expr: {rateExpr["Expr"]}}}\n')
                
            if 'efficiencies' in reaction_data:
                file.write('  efficiencies: ' + charge_coefficient + '\n')
            
        file.write('# -------------------------- END ---------------------------------/\n')
                

def main(input_file, mark_index, output_file, quiet):
    if quiet:
        logging.basicConfig(level=logging.ERROR, 
                            format='%(message)s',
                            handlers=[logging.FileHandler("parseKinetics.log"), logging.StreamHandler()])
    else:
        logging.basicConfig(level=logging.INFO, 
                            format='%(message)s',
                            handlers=[logging.FileHandler("parseKinetics.log"), logging.StreamHandler()])

    logging.info("\n>>>>>>>> Reading REACTION section ... <<<<<<<<")
    species, variables, reactions = read_reactions(input_file)
    
    charged_species = [s for s in species if s == 'e' or s == 'E' or '^-' in s or '^+' in s]
    charge_coefficient = "{" + ", ".join([f"{s}: 0." for s in charged_species]) + "}"
    
    parsed_variables = parse_variables(variables)
    parsed_variables = preprocess_variables(parsed_variables)
    
    logging.info("\nVariable list defined inline (preprocessed):")
    for key, value in parsed_variables.items():
        logging.info(f"{key}: {value}")
        
    logging.info("\n>>>>>>>> Processing reactions ... <<<<<<<<")

    parsed_reactions = parse_reactions(parsed_variables, reactions)
    
    write_yaml(charge_coefficient, parsed_reactions, markIndex=mark_index, filename=output_file)
    
    logging.info(f"\nWrite successfully with {len(parsed_reactions)} reactions and {len(species)} species:\n{species}")


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Parse and convert plasma kinetics in ZDPlasKin format.',
                                 usage='''python parsePlasKin.py --input <inputfile> [--output <outputfile>] [--markIndex <True/False>] [--quiet <True/False>]
                                 \n\nExample: python parsePlasKin.py --input "plasmaH2O2.inp" --output "parsedPlasKin.yaml" --markIndex True --quiet False''')

    parser.add_argument('--input', required=True, help='Input file name')
    parser.add_argument('--markIndex', type=bool, default=True, help='Mark reaction index in output file (default: True)')
    parser.add_argument('--output', default='parsedPlasKin.yaml', help='Output YAML file name (default: parsedPlasKin.yaml)')
    parser.add_argument('--quiet', type=bool, default=False, help='Set logging to quiet mode (default: False)')
    
    args = parser.parse_args()
        
    main(args.input, args.markIndex, args.output, args.quiet)
    
