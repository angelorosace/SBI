from Bio.PDB import NeighborSearch, Superimposer
from copy import deepcopy

MAX_CHAINS_IN_STRUCTURE = 10

# Recursively check all possible structures given the interacting chains that are
# compatible with the stoichiometry and do not have clashes after superimposing.
def recursively_add_chains_to_structure(model, interactions, stoichiometry=None, rmsd=0):
    if stoichiometry is None and len(list(model.get_chains())) >= MAX_CHAINS_IN_STRUCTURE:
        return model, rmsd

    # Keep track of the model that was determined as the best by lowest RMSD
    best_model = None
    best_rmsd = float('inf')

    # Check each chain in the model for the best additonal chain, if there is one that doesn't clash
    for chain in [chain.copy() for chain in model.get_chains()]:

        # Check each interaction structure (from the input files) that the chain is found in
        for structure in [struct.copy() for struct in interactions[chain.id]]:
            matching_chain = list(filter(lambda c: c.id == chain.id, structure.get_chains()))[0]
            non_matching_chain = list(filter(lambda c: c.id != chain.id, structure.get_chains()))[0]

            # Do not attempt to add the chain if it is not compatible with given stoichiometry
            if stoichiometry and not is_compatible_with_stoichiometry(model, non_matching_chain, stoichiometry):
                continue

            rmsd = superimpose_chain(chain, matching_chain, non_matching_chain)

            # As long as the rotated/translated chain does not clash with anything else in the
            # model, add it and recursively check other chains
            if not has_clashes_with_structure(model, non_matching_chain.get_atoms()):
                resulting_model, resulting_rmsd = recursively_add_chains_to_structure(
                    add_chain_to_model(model, non_matching_chain),
                    interactions,
                    stoichiometry,
                    rmsd
                )

                # If this model is better than any seen so far, save it
                if resulting_rmsd < best_rmsd:
                    best_model = resulting_model
                    best_rmsd = resulting_rmsd

    # Return the best model and RMSD if one is found. If all chains were
    # incompatible or had clashes, return the last seen model
    return (best_model, best_rmsd) if best_model else (model, rmsd)


# The model and chain are compatible if adding this chain id onto the model
# is still within the boundaries set in the stoichiometry. If the chain id is not
# in the stoichiometry, it is considered as not valid to be in the model.
def is_compatible_with_stoichiometry(structure, chain, stoichiometry):
    chain_id = chain.id
    chains = list(structure.get_chains())
    current_chain_count = len([chain for chain in chains if chain.id == chain_id])
    return chain_id in stoichiometry and \
           current_chain_count + 1 <= stoichiometry[chain_id]


# Superimpose the interacting structure on the chain that matches what is in the model
def superimpose_chain(chain, matching_chain, non_matching_chain):
    moving_atoms = list(matching_chain.get_atoms())
    fixed_atoms = list(chain.get_atoms())

    # After superimposing the matching chains, apply the rotation and translation
    # on the part of the interacting structure that does not match what is on the model
    superimposer = Superimposer()
    superimposer.set_atoms(fixed_atoms, moving_atoms)
    superimposer.apply(list(non_matching_chain.get_atoms()))
    return superimposer.rms


# Checks distance between all alpha carbons in the structure and chain atoms for
# any clashes (where the two atoms are within a certain distance).
def has_clashes_with_structure(structure, atoms, clash_distance=2):
    is_alpha_carbon = lambda atom: atom.get_name() == 'CA'
    structure_alpha_carbon_atoms = list(filter(is_alpha_carbon, structure.get_atoms()))
    chain_alpha_carbon_atoms = list(filter(is_alpha_carbon, atoms))

    searcher = NeighborSearch(structure_alpha_carbon_atoms)
    for atom in chain_alpha_carbon_atoms:
        clashes = searcher.search(atom.get_coord(), clash_distance, 'A')
        if len(clashes) > 0:
            return True

    return False


# Returns a copy of the model with the chain added
def add_chain_to_model(model, chain):
    model_copy = deepcopy(model)
    model_object = list(model_copy.get_models())[0]
    model_object.add(chain)
    return model_copy