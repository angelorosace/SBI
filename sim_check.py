from Bio import pairwise2

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
    "TYR": "Y",
    # DNA
    " DC": "C",
    " DG": "G",
    " DA": "A",
    " DT": "T",
    # RNA
    "G": "G",
    "U": "U",
    "C": "C",
    "A": "A"
}


def chain_similarity(chain1, chain2):
    """
    Takes two protein chain sequences and returns true if their sequences have > 95% similarity and false if not
    """
    # Get the sequences of the chains passed
    seq1, seq2 = get_sequence(chain1), get_sequence(chain2)
    if seq1 and seq2:
        # Align the aa sequences to see if there is similiarity
        alignment = pairwise2.align.globalxx(seq1, seq2)
        score = alignment[0][2]
        length = max(len(seq1), len(seq2))
        sim_score = score / length
        if sim_score > 0.95:
            return structure_similarity(chain1, chain2)
        else:
            return False


def structure_similarity(chain1, chain2):
    # get overlapping atoms between chains and calculate the percentage of similarity of the chains
    # based on their structure
    overlapping_atoms = get_overlapping_atoms(chain1, chain2)
    length = min(len([atom for atom in chain1.get_atoms()]), len([atom for atom in chain2.get_atoms()]))
    percent_similarity = overlapping_atoms / length * 100
    return percent_similarity >= 80


def get_overlapping_atoms(chain1, chain2):
    atoms1 = [atom for atom in chain1.get_atoms()]
    atoms2 = [atom for atom in chain2.get_atoms()]
    length = min(len(atoms1), len(atoms2))
    atom_coordinates1 = [atom.get_coord() for atom in atoms1[:length]]
    atom_coordinates2 = [atom.get_coord() for atom in atoms2[:length]]
    atom_overlapping = 0
    i = 0
    # check for each pair of atoms in the two chains if they overlap
    while i < length:
        j = 0
        coords_overlapping = 0
        # check if all the coordinates for all the dimension correspond between the two atoms
        while j < 3:
            if atom_coordinates1[i][j] == atom_coordinates2[i][j]:
                coords_overlapping += 1
            j += 1
        if coords_overlapping == 3:
            atom_overlapping += 1
        i += 1
    return atom_overlapping


def get_sequence(chain):
    seq = ''
    for res in chain.get_residues():
        res_name = res.get_resname()
        if res_name in residue_map:
            map_res = residue_map[res_name]
            seq += map_res
    return seq
