# *![chem cat logo](demo_images/chem_cat_logo.svg)* ChemCat
Concatenate chemical structures easily

## Assembling chemical structures has never been so simple
Building input files is part of the theoretical chemist's daily work.

In numerous research projects, the most painful step is assembling chemical bricks together.
In practice, input files are carefully concatenated together using various GUI-based programs (free or not) in a time-consuming and not so controlled way.

This tool aims to provide with a simple yet powerful, command-line based solution for concatenating input files in a controlled way, without additional hassle. Indeed, the assemblage depends only on simple and chemically relevant parameters.

### Parameters

- 5 atomic landmarks: `base1`, `base2`, `new1`, `new2`, `new3`

![atomic landmarks](demo_images/landmarks.svg)

- Bond vector (base1 -> new1)

![bond vector](demo_images/bond_vector.svg)

- Bond angle (base1, new1, new2)

![bond angle](demo_images/bond_angle.svg)

- Bond dihedral angle (base2, base1, new1, new2)

![bond dihedral angle](demo_images/bond_dihedral.svg)

- New molecular dihedral angle (base1, new1, new2, new3)

![new molecular dihedral angle](demo_images/new_dihedral.svg)

And that's it! Only the bare minimum (but complete) set of convenient parameters for merging two structures together.

Definitely, assembling chemical structures has never been so simple! Finally a `cat` command for chemistry...

## A convenient builder for adsorption

**TODO** `Simple interface for adding adsorbates + add images`

## A powerful Python3 API

**TODO** `As a python3 module...`
```Python3
import chemcat
```

## Installation

Simply download this repository, and look for the `chemcat.py` script. It provides with a command-line interface and can be used as a Python3 module.

In the future, a `setup.py` script will be provided for an even easier installation:
Setuptool installation:
- Download this repository
- Execute `python3 ChemCat/setup.py`
- Voilà!

## Q&A
### Is this compatible with periodic boundary conditions?
Yes! You can combine input files with or without periodic boundary conditions. 

Note that the combined unit cell will be defined by ASE, and it will most likely use the unit cell from the first input file.

### Are VASP output files compatible with Molden?
Indeed, Molden does not properly handle VASP files written in Cartesian coordinates (with the `Cartesian` keyword). Therefore, by default, output files are written in direct coordinates if PBC are detected.

As a consequence, output files generated with this script should be compatible with Molden!

Note that you can force writing cartesian coordinates by specifying the `--force-cartesian` option for the command-line interface.

### I want periodicity in one direction only, can I use this script?
Probably not. ASE does not currently handle PBC in some directions only. If such systems are detected, a warning message will be displayed...

### How can I adsorb a molecule on a bridge or fcc/hcp site?
Simply provide with a bond vector not aligned with the Z-axis (just remember to define the bond angle accordindly).

Note: for the command line interface, simply use the `-bx` or `-by` options (see `python3 chemcat.py -h` for help)

### The final VASP file is no longer compatible with my POTCAR
By default, ASE is conserving the concatenation order. But you can change that behavior and re-sort all the atomic indexes to group them by atomic type, simply use the `--sort` option for the command-line interface.

### Can I have my final structure in the VASP5 file format?
Yes! ASE supports the VASP5 file format. Simply use the `--vasp5` option for the command-line interface.
