# SBI

## Tutorial

### Software Requirements

The following are required to be installed before running the program:

```
Python    >= 3.0
Biopython == 1.76
argparse  == 1.4.0
```

### Installation

To use this program, you may either clone this Git repo and call `macro_molecular_complex_builder.py` directly or by installing the `dist/macro_molecular_complex_builder-0.0.1.tar.gz` using the following command:

```
pip install dist/macro_molecular_complex_builder-0.0.1.tar.gz
```

### Required Command Line Arguments

```
--input_directory
      Location of the directory where the list of interaction PDB files are.

--output_file
      Location of where the output PDB model should be written.
```

### Optional Command Line Arguments

```
--chain-limit
      Default = 10
      Maximum number of chains to use without stoichiometry. This setting is ignored when stoichiometry is given.

--stoichiometry
      Default = None
      File path to TSV indicating the total number of each chain in the final complex. See /.../ for an example.

--verbose
      Default = False
      Flag to determine if the program should output runtime information.
```

### Example

In the directory `/4g83` (among others), there is a list of interaction PDB files. These interaction files can be run through the program to generate a protein model by running the following command:

`python -m macro_molecular_complex_builder --input_directory ./4g83 --output_file out.pdb --verbose`

When running with stoichiometry, the following command can be run:

`python -m macro_molecular_complex_builder --input_directory ./4g83 --output_file out.pdb --stoichiometry ./example/stoichiometry.tsv`