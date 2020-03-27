from Bio.PDB import *


class InputParser:

    parser = PDBParser(QUIET=1)
    atom_dict = {}
    seq_dict = {}
    chain_dict = {}
    dna_bases = ["DC", "DG", "DA", "DT"]
    dna_map = {
        "DC": "D",
        "DG": "G",
        "DA": "A",
        "DT": "T"
    }
    rna_bases = ["G", "U", "C", "A"]
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
        self.populate_dicts(self.generate_structures(interaction_list))
        self.populate_seq_dict()

    def generate_structures(self, lst):
        allpdbs = [self.parser.get_structure(filename, filename) for filename in lst]
        return allpdbs

    def populate_dicts(self, lst):
        for pdb in lst:
            for chain in pdb.get_chains():
                is_dna = True
                is_rna = True
                chain_id = chain.get_id()
                atoms = []
                for atom in chain.get_atoms():
                    res_name = atom.get_parent().get_resname().lstrip()
                    if res_name not in self.dna_bases:
                        is_dna = False
                    if res_name not in self.rna_bases:
                        is_rna = False
                    atoms.append(atom)
                if is_dna:
                    self.populate_sub_dict("DNA", self.atom_dict, chain_id, atoms)
                    self.populate_sub_dict("DNA", self.chain_dict, chain_id, chain)
                elif is_rna:
                    self.populate_sub_dict("RNA", self.atom_dict, chain_id, atoms)
                    self.populate_sub_dict("RNA", self.chain_dict, chain_id, chain)
                elif chain_id not in self.atom_dict.keys():
                    self.chain_dict[chain_id] = chain
                    self.atom_dict[chain_id] = atoms

    @staticmethod
    def populate_sub_dict(seq_type, dictionary, chain_id, elem):
        if seq_type not in dictionary:
            dictionary[seq_type] = {}
            dictionary[seq_type][chain_id] = elem
        elif chain_id not in dictionary[seq_type].keys():
            dictionary[seq_type][chain_id] = elem

    def populate_seq_dict(self):
        for chain in self.atom_dict:
            if chain != "DNA" and chain != "RNA":
                seq = self.get_seq_data(self.atom_dict[chain])
                self.seq_dict[chain] = seq
            else:
                self.populate_sub_seq_dict(chain)

    def populate_sub_seq_dict(self, chain):
        for xna_chain in self.atom_dict[chain]:
            seq = self.get_seq_data(self.atom_dict[chain][xna_chain], chain)
            if chain not in self.seq_dict:
                self.seq_dict[chain] = {}
                self.seq_dict[chain][xna_chain] = seq
            else:
                self.seq_dict[chain][xna_chain] = seq

    def get_seq_data(self, atoms, mapping=None):
        seq = ""
        tmp = -1
        map = self.residue_map
        if mapping:
            if mapping == "DNA":
                map = self.dna_map
            else:
                map = self.rna_bases
        for atom in atoms:
            res = atom.get_parent()
            res_name = res.get_resname().lstrip()
            if type(map) == dict:
                map_res = map[res_name]
            else:
                map_res = res_name
            res_num = res.get_id()[1]
            if res_name in map:
                if seq == "":
                    seq += map_res
                    tmp = res_num
                else:
                    if tmp != res_num:
                        seq += map_res
                        tmp = res_num
        return seq
