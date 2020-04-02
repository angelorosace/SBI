from Bio.PDB import NeighborSearch, Superimposer
from copy import deepcopy
from printer import prnt
from sim_check import chain_similarity
from datetime import datetime

seen_interactions = []

# Recursively check all possible structures given the interacting chains that are
# compatible with the stoichiometry and do not have clashes after superimposing.
def recursively_add_chains_to_structure(model, interactions, max_chains, stoichiometry=None, rmsd=0):
    prnt('Entering recursion with model of %d chains, %d interaction chains' % (len(list(model.get_chains())), len(interactions.keys())))

    if len(list(model.get_chains())) >= max_chains:
        prnt('Too many chains in the model (>=%d)' % max_chains)
        return model, rmsd

    # Keep track of the model that was determined as the best by lowest RMSD
    best_model = None
    best_rmsd = float('inf')

    # Check each chain in the model for the best additonal chain, if there is one that doesn't clash
    for chain in [chain.copy() for chain in model.get_chains()]:
        prnt('Checking chain ID %s' % chain.id)

        # Check each interaction structure (from the input files) that the chain is found in
        for structure in [struct.copy() for struct in interactions[chain.id]]:
            matching_chain = list(filter(lambda c: c.id == chain.id, structure.get_chains()))[0]
            non_matching_chain = list(filter(lambda c: c.id != chain.id, structure.get_chains()))[0]

            prnt('Matches chain ID %s with non-matching chain ID %s' % (matching_chain.id, non_matching_chain.id))

            # Do not attempt to add the chain if it is not compatible with given stoichiometry
            if stoichiometry and not is_compatible_with_stoichiometry(model, non_matching_chain, stoichiometry):
                prnt('Chain ID %s is not compatible with stoichiometry' % non_matching_chain.id)
                continue

            if not chain_similarity(matching_chain, non_matching_chain):
                interaction = ''.join(sorted([matching_chain.get_id(), non_matching_chain.get_id()]))
                if interaction not in seen_interactions:
                    seen_interactions.append(interaction)
                    rmsd = superimpose_chain(chain, matching_chain, non_matching_chain)

                    # As long as the rotated/translated chain does not clash with anything else in the
                    # model, add it and recursively check other chains
                    if not has_clashes_with_structure(model, non_matching_chain.get_atoms()):
                        resulting_model, resulting_rmsd = recursively_add_chains_to_structure(
                            add_chain_to_model(model, non_matching_chain),
                            interactions,
                            max_chains,
                            stoichiometry,
                            rmsd
                        )

                        prnt('Resulting RMSD of ', resulting_rmsd,'for model with chains',\
                             list(map(lambda c: c.id, resulting_model.get_chains())))

                        # If this model is better than any seen so far, save it
                        if resulting_rmsd < best_rmsd:
                            prnt('Resulting RMSD is lower than previous best RMSD of %f' % best_rmsd)
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

    prnt('Superimposing %d atoms' % len(moving_atoms))

    # After superimposing the matching chains, apply the rotation and translation
    # on the part of the interacting structure that does not match what is on the model
    superimposer = Superimposer()
    superimposer.set_atoms(fixed_atoms, moving_atoms)
    superimposer.apply(list(non_matching_chain.get_atoms()))
    return superimposer.rms


# Checks distance between all alpha carbons in the structure and chain atoms for
# any clashes (where the two atoms are within a certain distance).
def has_clashes_with_structure(structure, atoms, clash_distance=2):
    is_alpha_carbon_or_phosphate = lambda atom: atom.get_name() == 'CA' or atom.get_name() == "P"
    structure_alpha_carbon_atoms = list(filter(is_alpha_carbon_or_phosphate, structure.get_atoms()))
    chain_alpha_carbon_atoms = list(filter(is_alpha_carbon_or_phosphate, atoms))

    prnt('Checking alpha carbon clashes (<%f A) between %d atoms in the structure and %d atoms in the chain'\
        % (clash_distance, len(structure_alpha_carbon_atoms), len(chain_alpha_carbon_atoms)))

    searcher = NeighborSearch(structure_alpha_carbon_atoms)
    for atom in chain_alpha_carbon_atoms:
        clashes = searcher.search(atom.get_coord(), clash_distance, 'A')
        if len(clashes) > 0:
            prnt('Found %d clash(es)' % len(clashes))
            return True

    prnt('No clashes found')
    return False


# Returns a copy of the model with the chain added
def add_chain_to_model(model, chain):
    prnt('Adding chain ID %s to model' % chain.id)
    model_copy = deepcopy(model)
    model_object = list(model_copy.get_models())[0]
    model_object.add(chain)
    return model_copy

