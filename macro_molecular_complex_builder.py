from argument_parser import parse_arguments
from printer import set_is_verbose

if __name__ == '__main__':
    inputs = parse_arguments()
    if inputs['verbose']:
      set_is_verbose()