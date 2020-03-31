# SBI

## Tutorial

### Software Requirements

The following are required to be installed before running the program:

```
Python    >= 3.0
Biopython == 1.76
argparse  == 1.4.0
```

### Required Command Line Arguments

```
--input_directory
      ...

--output_directory
      ...
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

In the directory `/example`, there is a list of interaction PDB files. These interaction files can be run through the program to generate a protein model by running the following command:

`python3 macro_molecular_complex_builder.py --input_directory ./example --output_directory ./output --verbose`

When running with stoichiometry, the following command can be run:

`python3 macro_molecular_complex_builder.py --input_directory ./example --output_directory ./output --stoichiometry ./example/stoichiometry.tsv`