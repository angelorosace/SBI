import argparse
import os

def parse_arguments():
    parser = argparse.ArgumentParser(description='SBI Complex Builder')

    parser.add_argument('--input_directory', action='store', default=None, dest='input_directory', help='Directory location of PDB files')
    parser.add_argument('--output_directory', action='store', default=None, dest='output_directory', help='Directory location to save output files')
    args = parser.parse_args()

    pdb_input_files = get_pdb_files_in_directory(args.input_directory)
    if not len(pdb_input_files):
        raise Exception('%s does not contain any .pdb files' % args.input_directory)

    return pdb_input_files, args.output_directory

def get_pdb_files_in_directory(directory):
    if not os.path.isdir(directory):
        raise Exception('%s is not a directory' % directory)

    return [file for file in os.listdir(directory) if file.endswith('.pdb')]