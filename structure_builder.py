from Bio.PDB import Structure, Model, NeighborSearch, PDBParser, Superimposer

from printer import prnt

MAX_CHAINS_IN_STRUCTURE = 20
MINIMUM_RMSD = 10

def testing():
    structure_id = 'complex0'
    model_id = 'model0'

    prnt('Starting testing', 'here', 'hello')

    structure = Structure.Structure(structure_id)
    model = Model.Model(model_id)
    structure.add(model)

    parser = PDBParser()
    pdb_str = parser.get_structure('6mgh', '6gmh.pdb')

    chains = list(pdb_str.get_chains())
    model.add(chains[0])

    return structure, chains


def has_clashes_with_structure(structure, atoms, clash_distance=2, minimum_atoms_for_clash=10):
    searcher = NeighborSearch(list(structure.get_atoms()))

    for atom in atoms:
        clashes = searcher.search(atom.get_coord(), clash_distance, 'A')
        if len(clashes) >= minimum_atoms_for_clash:
            return True

    return False


# The structure and chain are compatible if adding this chain id onto the structure
# is still within the boundaries set in the stoichiometry. If the chain id is not
# in the stoichiometry, it is considered as not valid to be in the structure.
def is_compatible_with_stoichiometry(structure, chain, stoichiometry):
    chain_id = chain.id
    chains = list(structure.get_chains())
    current_chain_count = len([chain for chain in chains if chain.id == chain_id])
    return chain_id in stoichiometry and \
           current_chain_count + 1 <= stoichiometry[chain_id]


# Recursively check all possible structures given the chains that are compatible with
# the stoichiometry and do not have clashes after superimposing.
def recursively_add_chains_to_structure(structure, chains, stoichiometry=None, rmsd=0):
    if stoichiometry is None and len(list(structure.get_chains())) >= MAX_CHAINS_IN_STRUCTURE:
        return structure, rmsd

    # Keep track of the structure that was determined as the best by lowest RMSD
    best_structure = None
    best_rmsd = float('inf')

    for chain in [chain.copy() for chain in chains]:
        # Do not attempt to use this chain if it is not compatible with given stoichiometry
        if stoichiometry and not is_compatible_with_stoichiometry(structure, chain, stoichiometry):
            continue

        structure_atoms = list(structure.get_atoms())
        chain_atoms = list(chain.get_atoms())

        # Ensure that the lists of atoms are of the same length
        #
        #       TODO - Figure out a better way to do this
        #
        imposition_length = min(len(structure_atoms), len(chain_atoms))
        fixed_atoms = structure_atoms[:imposition_length]
        moving_atoms = chain_atoms[:imposition_length]

        # Superimpose and transform the atoms in the moving_atoms list to the best rotation/translation
        superimposer = Superimposer()
        superimposer.set_atoms(fixed_atoms, moving_atoms)
        superimposer.apply(moving_atoms)

        # Continue building the structure with the added chain if it does not clash
        if superimposer.rms < MINIMUM_RMSD and not has_clashes_with_structure(structure, moving_atoms):
            structure_copy = structure.copy()

            # Add the chain to the model with a new chain id of all chain ids in the model
            model = list(structure_copy.get_models())[0]
            chain.id = ''.join(map(lambda chain: chain.id, list(model.get_chains()))) + chain.id
            model.add(chain)

            resulting_structure, resulting_rmsd = recursively_add_chains_to_structure(structure_copy, chains, stoichiometry, superimposer.rms)

            # Update the best structure found so far if the best substructure
            # and current superimposed structure have a lower RMSD
            if resulting_rmsd < best_rmsd:
                best_structure = resulting_structure
                best_rmsd = resulting_rmsd

    # Return the best structure and RMSD if one is found. If all chains were
    # incompatible or had clashes, return the last seen structure
    return (best_structure, best_rmsd) if best_structure else (structure, rmsd)