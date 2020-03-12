from Bio.PDB import Structure, Model, Chain, NeighborSearch, PDBParser, Superimposer

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

def has_clashes_with_structure(structure, chain, clash_distance=2, minimum_atoms_for_clash=10):
    searcher = NeighborSearch(list(structure.get_atoms()))

    for atom in chain.get_atoms():
        clashes = searcher.search(atom.get_coord(), clash_distance, 'A')
        if len(clashes) >= minimum_atoms_for_clash:
            return True

    return False

def recursively_add_chains_to_structure(structure, chains):
    if not len(chains):
        return structure

    superimposer = Superimposer()
    
    superimposer.set_atoms(list(structure.get_atoms()), list(chains[0].get_atoms())) 

    print(superimposer.rms)

    # chain = chains[0]
    # if not has_clashes_with_structure(structure, chain):







