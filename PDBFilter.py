from Bio.PDB import *
import os


class PDBFilter:

    """
    Class used to provide a filter for reading pdb files

    Attributes

    ----------
    parser : None
        A variable to store an instance of the PDBParser Object
    struct : List
        A list used to store the structure objects of the different proteins in the input.
    residue_map : dict
        A dictionary used for mapping residues names (keys) to their one-letter name (values).
    pdb_seq_dict : dict
        A dictionary used to store the chains (keys) of the input pdb files and their sequences (values).
    pdb_atom_dict: dict
        A dictionary used to store the chains (keys) of the input pdb files and their atoms (values).

    Methods

    -------
    build_structures(pdb_list=List)
        Extracts protein structures from the pdb files in the pbd_list.
    get_pdb_name_from_path(path=str)
        Extracts pbd name from file path
    get_seq_data(structures=List)
        Populates the pdb_seq_dict with data for the chains of the structures and their sequences.
    get_atom_data(structures=List)
        Populates the pdb_atom_dict with data for the chains of the structures and their atoms.
     write_pairings_files(pdb_list)
        Writes files containing the atoms for all the pairs of chains for each protein in the input list
    pair_chain(pdb=str)
        Pairs up all the different chains within the same protein and writes a file for each pair.
        The file contains the atoms for the two chains in the pair.
    write_interaction_file(pdb=str,chain1=List,chain2=List)
        Checks if the repository for storing the interaction files exists. If it exists it writes the files and stores
        them into the directory, if not it first creates the directory and then creates the files.
    write_interaction_file_aux(pdb=str,chain1=List,chain2=List)
        Formats the data and writes it in the output interaction files.
    """
    parser = None
    structs = []
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
        "TYR": "Y"
    }
    pdb_seq_dict = {}
    pdb_atom_dict = {}

    def __init__(self, pdb_list):
        self.parser = PDBParser(QUIET=1)
        structs = self.build_structures(pdb_list)
        self.get_seq_data(structs)
        self.get_atom_data(structs)

    def build_structures(self, pdb_list):
        structures = []
        for i, pdb in enumerate(pdb_list):
            structures.append(self.parser.get_structure(self.get_pdb_name_from_path(pdb), pdb))
        return structures

    @staticmethod
    def get_pdb_name_from_path(path):
        lst = path.split("/")
        tmp_name = lst[-1]
        name = tmp_name.split(".")
        return name[0]

    def get_seq_data(self, structures):
        for i, structure in enumerate(structures):
            tmp_dict = {}
            for model in structure:
                for chain in model:
                    seq = ""
                    for residue in chain:
                        if residue.get_resname() in self.residue_map.keys():
                            seq += self.residue_map[residue.get_resname()]
                        tmp_dict[chain.get_id()] = seq
            self.pdb_seq_dict[structure.get_id()] = tmp_dict

    def get_atom_data(self, structures):
        for i, structure in enumerate(structures):
            tmp_dict = {}
            for model in structure:
                j = 0
                for k, chain in enumerate(model):
                    atoms = []
                    for i, residue in enumerate(chain):
                        for atom in residue:
                            coords = atom.get_coord()
                            j += 1
                            atoms.append(["ATOM", j, atom.get_name(), atom.get_parent().get_resname(), chain.get_id(), i, coords[0], coords[1], coords[2], atom.get_occupancy(), atom.get_bfactor()])
                    tmp_dict[k] = atoms
            self.pdb_atom_dict[structure.get_id()] = tmp_dict

    def write_pairings_files(self, pdb_names):
        for pdb_name in pdb_names:
            self.pair_chain(pdb_name)

    def pair_chain(self, pdb):
        for i, chain in enumerate(self.pdb_atom_dict[pdb]):
            j = i+1
            while j < len(self.pdb_atom_dict[pdb]):
                atoms_first_chain = self.pdb_atom_dict[pdb][chain]
                atoms_second_chain = self.pdb_atom_dict[pdb][j]
                self.write_interaction_file(pdb, atoms_first_chain, atoms_second_chain)
                j += 1
        print("The pairing for all the chains in  the %s pdb file have been generated" % pdb)

    def write_interaction_file(self, pdb, chain1, chain2):
        if os.path.exists(pdb + '_interactions'):
            if not os.path.isfile(pdb + '_interactions/' + chain1[0][4] + chain2[0][4] + '.pdb'):
                self.write_interaction_file_aux(pdb, chain1, chain2)
        else:
            os.makedirs(pdb + '_interactions')
            self.write_interaction_file_aux(pdb, chain1, chain2)

    @staticmethod
    def write_interaction_file_aux(pdb, chain1, chain2):
        file = open(pdb + '_interactions/' + chain1[0][4] + chain2[0][4] + '.pdb', 'w')
        for atom in chain1:
            stri = ""
            for elem in atom:
                stri += str(elem) + "\t"
            file.write(stri)
            file.write("\n")
        for atom in chain2:
            stri = ""
            for elem in atom:
                stri += str(elem) + "\t"
            file.write(stri)
            file.write("\n")
        file.close()
