from argument_parser import parse_arguments
from printer import set_is_verbose, prnt

from structure_builder import testing

if __name__ == '__main__':
    inputs = parse_arguments()
    if inputs['verbose']:
      set_is_verbose()

    testing()

    prnt('done for', 'now')