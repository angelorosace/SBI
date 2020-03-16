from Bio.PDB import *

class InputParser:

    parser = PDBParser(QUIET=1)
    atom_dict = {}
    seq_dict = {}
    residue_map = {
        "ALA": "A",
        "CYS": "C",
        "ASP": "D",
        "GLU": "E",
        "PHE": "F",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LYS": "K",
        "LEU": "L",
        "MET": "M",
        "ASN": "N",
        "PRO": "P",
        "GLN": "Q",
        "ARG": "R",
        "SER": "S",
        "THR": "T",
        "VAL": "V",
        "TRP": "W",
        "TYR": "Y"
    }

    def __init__(self, interaction_list):
        self.populate_atom_dict(self.generate_structures(interaction_list))
        self.populate_seq_dict()

    def generate_structures(self, lst):
        allpdbs = [self.parser.get_structure(filename, filename) for filename in lst]
        return allpdbs

    def populate_atom_dict(self, lst):
        for pdb in lst:
            for chain in pdb.get_chains():
                chain_id = chain.get_id()
                atoms = []
                for atom in chain.get_atoms():
                    atoms.append(atom)
                if chain_id not in self.atom_dict.keys():
                    self.atom_dict[chain_id] = atoms

    def populate_seq_dict(self):
        for chain in self.atom_dict:
            self.get_seq_data(chain, self.atom_dict[chain])

    def get_seq_data(self, chain, atoms):
        seq = ""
        for atom in atoms:
            res = atom.get_parent().get_resname()
            map_res = self.residue_map[res]
            if res in self.residue_map.keys():
                if seq == "":
                    seq += map_res
                else:
                    if not seq[-1] == map_res:
                        seq += map_res
        self.seq_dict[chain] = seq

