from Bio.PDB import Structure, Model, Chain, NeighborSearch, PDBParser, Superimposer
from Bio import pairwise2

MAX_CHAINS_IN_STRUCTURE = 100

def testing():
    structure_id = 'complex0'
    model_id = 'model0'

    structure = Structure.Structure(structure_id)
    model = Model.Model(model_id)
    structure.add(model)

    parser = PDBParser()
    pdb_str = parser.get_structure('6mgh', '6gmh.pdb')

    a = None
    for c in pdb_str.get_chains():
        if c.id != 'A':
            model.add(c)
        else:
            a = c

    return structure, a

    # searcher = NeighborSearch(list(structure.get_atoms()))

    # for atom in a.get_atoms():
    #     clashes = searcher.search(atom.get_coord(), 2, 'A')

    #     if len(clashes):
    #         print(atom.id, clashes, a.id)

def has_clashes_with_structure(structure, atoms, clash_distance=2, minimum_atoms_for_clash=10):
    searcher = NeighborSearch(list(structure.get_atoms()))

    for atom in atoms:
        clashes = searcher.search(atom.get_coord(), clash_distance, 'A')
        if len(clashes) >= minimum_atoms_for_clash:
            return True

    return False

# def get_sequence(residues):
#     residue_mapping = {
#         "ALA": "A",
#         "CYS": "C",
#         "ASP": "D",
#         "GLU": "E",
#         "PHE": "F",
#         "GLY": "G",
#         "HIS": "H",
#         "ILE": "I",
#         "LYS": "K",
#         "LEU": "L",
#         "MET": "M",
#         "ASN": "N",
#         "PRO": "P",
#         "GLN": "Q",
#         "ARG": "R",
#         "SER": "S",
#         "THR": "T",
#         "VAL": "V",
#         "TRP": "W",
#         "TYR": "Y"
#     }

#     sequence = ''
#     for residue in residues:
#         residue_name = residue.get_resname()
#         if residue_name in residue_mapping:
#             sequence += residue_mapping[residue_name]

#     return sequence




# TODO
#   Figure out structure of stoichiometry
#   Loop through every current chain in structure to see if adding this chain is valid
def is_compatible_with_stoichiometry(structure, chain, stoichiometry):
    return True

def recursively_add_chains_to_structure(structure, chains, stoichiometry=None):
    if len(list(structure.get_chains())) >= MAX_CHAINS_IN_STRUCTURE:
        return structure, 0

    # Keep track of the structure that was determined as the best by lowest RMSD
    best_structure = None
    best_rmsd = float('inf')

    for chain in chains:
        if stoichiometry is not None and not is_compatible_with_stoichiometry(structure, chain, stoichiometry):
            continue

        # Copy the two objects to not modify the originals
        structure_atoms = list(structure.copy().get_atoms())
        chain_atoms = list(chain.copy().get_atoms())

        # Ensure that the lists of atoms are of the same length
        imposition_length = min(len(structure_atoms), len(chain_atoms))
        fixed_atoms = structure_atoms[:imposition_length]
        moving_atoms = chain_atoms[:imposition_length]

        # Superimpose and move the atoms in the moving_atoms list to the best rotation/translation
        superimposer = Superimposer()
        superimposer.set_atoms(fixed_atoms, moving_atoms)
        superimposer.apply(moving_atoms)

        if not has_clashes_with_structure(structure, moving_atoms):
            structure_copy = structure.copy()
            structure_copy.add(chain)
            resulting_structure, rmsd = recursively_add_chains_to_structure(structure_copy, chains)

            # Update the best structure found so far if the best substructure
            # and current superimposed structure have a lower RMSD

            # ??? Do we need to sum the total rmsd or would superimposing take that into account?
            total_rmsd = rmsd + superimposer.rms
            if total_rmsd < best_rmsd:
                best_structure = resulting_structure
                best_rmsd = total_rmsd

    return best_structure, best_rmsd




        # best_rmsd = float('inf')
        # best_superimposition = None
        # best_superimposed_moving_atoms = None

        # window_length = len(chain_atoms)
        # sliding_window_atoms = structure_atoms
        # static_atoms = chain_atoms
        # sliding_window_is_structure = True
        # if len(chain_atoms) > len(structure_atoms):
        #     window_length = len(structure_atoms)
        #     sliding_window_atoms = chain_atoms
        #     static_atoms = structure_atoms
        #     sliding_window_is_structure = False

        # number_of_windows = len(sliding_window_atoms) - len(static_atoms) + 1

        # for start_index in range(number_of_windows):
        #     atoms_in_window = sliding_window_atoms[start_index : start_index + window_length]

        #     superimposer = Superimposer()

        #     fixed_atoms, moving_atoms = (atoms_in_window, static_atoms) \
        #                                 if sliding_window_is_structure else (static_atoms, atoms_in_window)
        #     superimposer.set_atoms(fixed_atoms, moving_atoms)

        #     if start_index % 1000 == 0:
        #         print(start_index / len(sliding_window_atoms))

        #     if superimposer.rms < best_rmsd:
        #         best_rmsd = superimposer.rms
        #         best_superimposition = superimposer
        #         best_superimposed_moving_atoms = moving_atoms

        # print('applying')
        # best_superimposition.apply(best_superimposed_moving_atoms)