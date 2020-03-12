class InputParser:

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
        self.populate_atom_dict(interaction_list)
        self.populate_seq_dict()

    def populate_atom_dict(self, lst):
        for interaction in lst:
            file = open(interaction, "r")
            values = [line.split() for line in file]
            chain_ids = list(set([value[4] for value in values]))
            chain_ids.sort()
            if chain_ids[0] not in self.atom_dict.keys():
                first_chain = [atom for atom in values if atom[4] == chain_ids[0]]
                self.atom_dict[chain_ids[0]] = first_chain
            if chain_ids[1] not in self.atom_dict.keys():
                second_chain = [atom for atom in values if atom[4] == chain_ids[1]]
                self.atom_dict[chain_ids[1]] = second_chain
            file.close()

    def populate_seq_dict(self):
        for chain in self.atom_dict:
            self.get_seq_data(chain, self.atom_dict[chain])

    def get_seq_data(self, chain, atoms):
        seq = ""
        for residue in atoms:
            if residue[3] in self.residue_map.keys():
                if seq == "":
                    seq += self.residue_map[residue[3]]
                else:
                    if not seq[-1] == self.residue_map[residue[3]]:
                        seq += self.residue_map[residue[3]]
        self.seq_dict[chain] = seq
