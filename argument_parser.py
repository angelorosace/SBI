import argparse
import os
import csv

def parse_arguments():
    parser = argparse.ArgumentParser(description='SBI Complex Builder')
    parser.add_argument('--input_directory', action='store', default=None, dest='input_directory', help='Directory of chain interactions')
    parser.add_argument('--input_pdb', action='store', default=None, dest='input_pdb', help='Location of PDB to test')
    parser.add_argument('--output_directory', action='store', default=None, dest='output_directory', help='Directory location to save output files')
    parser.add_argument('--verbose', action='store_true', default=False, dest='verbose', help='The program should print messages during execution')
    parser.add_argument('--stoichiometry', action='store', default=None, dest='stoichiometry', help='Location of TSV file used as stoichiometry of the complex')
    parser.add_argument('--chain-limit', action='store', default=10, dest='chain_limit', help='Maximum number of chains in the final complex')
    args = parser.parse_args()

    if args.input_pdb is None and args.input_directory is None:
        raise Exception('One of --input_pdb and --input_directory must be given')

    if args.input_pdb is not None and args.input_directory is not None:
        raise Exception('Cannot run with both a PDB file and a directory of chain interactions')

    if args.output_directory is None:
        raise Exception('Output directory is required')

    inputs = {
        'output_directory': args.output_directory,
        'verbose': args.verbose,
        'stoichiometry': None,
        'chain_limit': int(args.chain_limit)
    }

    if args.stoichiometry is not None:
        if not os.path.isfile(args.stoichiometry):
            raise Exception('%s does not exist' % args.stoichiometry)

        inputs['stoichiometry'] = {}
        with open(args.stoichiometry, 'r') as tsv_file:
            reader = csv.reader(tsv_file)
            for line in reader:
                split_line = line[0].split('\t')
                if len(split_line) != 2:
                    raise Exception('Stoichiometry file is incorrectly formatted. Each line must be the chain id followed by the count delimited by a tab')

                chain_id = split_line[0]
                chain_count = int(split_line[1])
                inputs['stoichiometry'][chain_id] = chain_count

        inputs['chain_limit'] = sum(inputs['stoichiometry'].values())

    if args.input_pdb is not None:
        if not os.path.isfile(args.input_pdb):
            raise Exception('%s does not exist' % args.input_pdb)

        inputs['pdb'] = args.input_pdb

    if args.input_directory is not None:
        chain_interaction_files = get_pdb_files_in_directory(args.input_directory)
        if not len(chain_interaction_files):
            raise Exception('%s does not contain any .pdb files' % args.input_directory)

        inputs['chain_interactions'] = chain_interaction_files

    if not os.path.isdir(args.output_directory):
        os.mkdir(args.output_directory)

    return inputs

def get_pdb_files_in_directory(directory):
    if not os.path.isdir(directory):
        raise Exception('%s is not a directory' % directory)

    return [file for file in os.listdir(directory) if file.endswith('.pdb')]