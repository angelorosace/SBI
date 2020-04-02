from argument_parser import parse_arguments
from printer import set_is_verbose
from InputParser import InputParser
from structure_builder import recursively_add_chains_to_structure
from Bio.PDB import PDBIO

if __name__ == '__main__':
    inputs = parse_arguments()
    if inputs['verbose']:
        set_is_verbose()
    i = InputParser(inputs["chain_interactions"])
    if i.root:
        model, rmsd = recursively_add_chains_to_structure(i.structure_dict[i.root], i.interactions, inputs["chain_limit"], inputs["stoichiometry"])
    else:
        model, rmsd = recursively_add_chains_to_structure(i.structure_dict.items()[0], i.interactions, inputs["chain_limit"], inputs["stoichiometry"])
    io = PDBIO()
    io.set_structure(model)
    io.save(inputs["output_file"])
