from Bio.PDB import *
from Bio import pairwise2

class Filter:

    parser = PDBParser()
    pdb_seq_dict = {}

    def __init__(self, pdb_list):
        self.filter_pdb(pdb_list)

    def filter_pdb(self, pdb_list):
        for i, pdb in enumerate(pdb_list):
            chains = self.split_chains(pdb, i)
            residues = self.get_sequence(chains[1])
            self.pdb_seq_dict[chains[0]] = residues

    def split_chains(self, pdb, label):
        structure = self.parser.get_structure(label, pdb)
        modules = [module for module in structure]
        chains = [chain for chain in modules]
        return label, chains

    def get_sequence(self, chains):
        residues = [residue for residue in chains]
        return residues

    def check_similarity(self, seq1, seq2):
        alignment = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)
        score = alignment[2]
        length = alignment[4]
        similarity = score/length
        if similarity >= 0.95:
            print("discarded")

