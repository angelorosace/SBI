from argument_parser import parse_arguments

from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import PDBParser

from structure_builder import testing, recursively_add_chains_to_structure, recurse

if __name__ == '__main__':
    # inputs = parse_arguments()

    m, i = testing()

    model, rmsd = recurse(m, i)

    print(list(model.get_chains()))

    io=PDBIO()
    io.set_structure(model)
    io.save('out.pdb')