from Bio.PDB import *

parser = PDBParser(PERMISSIVE=1, QUIET=1)
#struct = parser.get_structure("A", "C:/Users/Angelo/IdeaProjects/SBI/6gmh_interactions/AB.pdb")
#print([model for model in struct])

files = ["C:/Users/Angelo/IdeaProjects/SBI/draft/3kuy_AB.pdb"]
allpdb = {parser.get_structure(filename, filename) for filename in files}
for pdb in allpdb:
    for chain in pdb.get_chains():
        print(chain.get_id())
        for atoms in chain.get_atoms():
            print(atoms.get_parent().get_resname())
            exit()
