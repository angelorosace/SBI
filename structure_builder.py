from Bio.PDB import Structure, Model, NeighborSearch, PDBParser, Superimposer
from copy import deepcopy

MAX_CHAINS_IN_STRUCTURE = 20
MAXIMUM_RMSD = 10

def testing():
    # structure_id = 'complex0'
    model_id = 'model0'

    # structure = Structure.Structure(structure_id)
    model = Model.Model(model_id)
    # structure.add(model)

    parser = PDBParser()

    ac = parser.get_structure('AC', 'AC.pdb')
    bc = parser.get_structure('BC', 'BC.pdb')
    bd = parser.get_structure('BD', 'BD.pdb')

    interactions = {
        'A': {ac},
        'B': {bc, bd},
        'C': {ac, bc},
        'D': {bd}
    }

    return ac, interactions

    # pdb_str = parser.get_structure('6gmh', '6gmh.pdb')

    # chains = list(pdb_str.get_chains())
    # model.add(chains[0])

    # d = {}

    # for c1 in chains:
    #     d[c1.id] = set()

    #     for c2 in chains:
    #         if c1 != c2:
    #             try:
    #                 if has_clashes_with_structure(c1, c2.get_atoms(), 8):
    #                     d[c1.id].add(c2.id)
    #             except:
    #                 pass

    # return model, chains, d


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


# The structure and chain are compatible if adding this chain id onto the structure
# is still within the boundaries set in the stoichiometry. If the chain id is not
# in the stoichiometry, it is considered as not valid to be in the structure.
def is_compatible_with_stoichiometry(structure, chain, stoichiometry):
    chain_id = chain.id
    chains = list(structure.get_chains())
    current_chain_count = len([chain for chain in chains if chain.id == chain_id])
    return chain_id in stoichiometry and \
           current_chain_count + 1 <= stoichiometry[chain_id]


def recurse(model, interactions, rmsd=0):
    print('\n')
    if len(list(model.get_chains())) >= MAX_CHAINS_IN_STRUCTURE:
        return model, rmsd

    best_model = None
    best_rmsd = float('inf')

    for chain in [chain.copy() for chain in model.get_chains()]:
        possible_structures = interactions[chain.id]

        print('Checking chain %s' % chain.id, possible_structures)

        for structure in [s.copy() for s in possible_structures]:
            matching_chain = list(filter(lambda c: c.id == chain.id, structure.get_chains()))[0]
            non_matching_chain = list(filter(lambda c: c.id != chain.id, structure.get_chains()))[0]
            moving_atoms = list(matching_chain.get_atoms())
            fixed_atoms = list(chain.get_atoms())

            superimposer = Superimposer()
            superimposer.set_atoms(fixed_atoms, moving_atoms)
            superimposer.apply(list(non_matching_chain.get_atoms()))

            if not has_clashes_with_structure(model, non_matching_chain.get_atoms()):
                m = deepcopy(model)
                mdl = list(m.get_models())[0]
                mdl.add(non_matching_chain)

                print('Adding %s to' % non_matching_chain.id, list(mdl.get_chains()))
                print('Current RMSD %f' % superimposer.rms)

                resulting_model, resulting_rmsd = recurse(m, interactions, superimposer.rms)

                print('Resulting RMSD %f' % resulting_rmsd)
                print('Resulting model', list(resulting_model.get_chains()))

                if resulting_rmsd < best_rmsd:
                    best_model = resulting_model
                    best_rmsd = resulting_rmsd

    return (best_model, best_rmsd) if best_model else (model, rmsd)

# Recursively check all possible structures given the chains that are compatible with
# the stoichiometry and do not have clashes after superimposing.
def recursively_add_chains_to_structure(structure, chains, interactions, stoichiometry=None, rmsd=0):
    if stoichiometry is None and len(list(structure.get_chains())) >= MAX_CHAINS_IN_STRUCTURE:
        return structure, rmsd

    # Keep track of the structure that was determined as the best by lowest RMSD
    best_structure = None
    best_rmsd = float('inf')

    for chain in [chain.copy() for chain in chains]:
        # Do not attempt to use this chain if it is not compatible with given stoichiometry
        if stoichiometry and not is_compatible_with_stoichiometry(structure, chain, stoichiometry):
            continue

        # structure_atoms = list(structure.get_atoms())
        # chain_atoms = list(chain.get_atoms())

        # Ensure that the lists of atoms are of the same length
        #
        #       TODO - Figure out a better way to do this
        #
        # imposition_length = min(len(structure_atoms), len(chain_atoms))
        # fixed_atoms = structure_atoms[:imposition_length]
        # moving_atoms = chain_atoms[:imposition_length]

        # Superimpose and transform the atoms in the moving_atoms list to the best rotation/translation
        superimposer = Superimposer()
        superimposer.set_atoms(fixed_atoms, moving_atoms)
        superimposer.apply(moving_atoms)

        # Continue building the structure with the added chain if it does not clash
        # if superimposer.rms < MAXIMUM_RMSD and not has_clashes_with_structure(structure, moving_atoms):
        if not has_clashes_with_structure(structure, moving_atoms):
            print('adding %s - %f' % (chain.id, superimposer.rms))

            structure_copy = structure.copy()

            # Add the chain to the model with a new chain id of all chain ids in the model
            # model = list(structure_copy.get_models())[0]
            model = structure_copy
            chain.id = ''.join(map(lambda chain: chain.id, list(model.get_chains()))) + chain.id
            model.add(chain)

            resulting_structure, resulting_rmsd = recursively_add_chains_to_structure(structure_copy, chains, interactions, stoichiometry, superimposer.rms)

            # Update the best structure found so far if the best substructure
            # and current superimposed structure have a lower RMSD
            if resulting_rmsd < best_rmsd:
                best_structure = resulting_structure
                best_rmsd = resulting_rmsd

    # Return the best structure and RMSD if one is found. If all chains were
    # incompatible or had clashes, return the last seen structure
    return (best_structure, best_rmsd) if best_structure else (structure, rmsd)