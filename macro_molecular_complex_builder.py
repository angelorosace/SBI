from argument_parser import parse_arguments
import InputParser
from sim_check import get_novel_sequences
from structure_builder import build_model

if __name__ == '__main__':
    # inputs = parse_arguments()
    # inputs_parser = InputParser.InputParser(inputs)
    inputs_parser = InputParser.InputParser(["draftAB.pdb", "draftAC.pdb", "draftAD.pdb", "draftBC.pdb", "draftBD.pdb", "draftCD.pdb"])
    # get novel/unique chains
    novel_sequences = get_novel_sequences(inputs_parser.seq_dict)
    # get non similar chains from the initial input
    useful_chains = [inputs_parser.chain_dict[chain_id] for chain_id in inputs_parser.chain_dict if chain_id in novel_sequences]
    # build complex
    build_model(useful_chains)
