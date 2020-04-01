from Bio.PDB import *


class InputParser:

    parser = PDBParser(QUIET=1)
    structure_dict = {}
    interactions = {}

    def __init__(self, interaction_list):
        self.generate_structures(interaction_list)

    def generate_structures(self, lst):
        all_pdbs = [self.parser.get_structure(filename, filename) for filename in lst]
        for pdb in all_pdbs:
            self.structure_dict[pdb.get_id()] = pdb
        self.populate_interactions()
        return all_pdbs

    def populate_interactions(self):
        for file_name, interaction_structure in self.structure_dict.items():
            interaction_ids = [chain.get_id() for chain in interaction_structure.get_chains()]
            self.add_interaction(interaction_structure, interaction_ids[0])
            self.add_interaction(interaction_structure, interaction_ids[1])

    def add_interaction(self, interaction, interaction_id):
        if interaction_id not in self.interactions:
            self.interactions[interaction_id] = [interaction]
        else:
            if interaction not in self.interactions[interaction_id]:
                self.interactions[interaction_id].append(interaction)
