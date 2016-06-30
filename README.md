# ANISOFLOW

ANISOFLOW is still in development!!

## Documentation

To get an early documentation pdf, download it from [here](https://www.overleaf.com/read/trycqnfcynsp)

## Synopsis

The main philosophy of ANISOFLOW is to capture the anisotropies of the medium into the flow, for that purpose the program has three ways to carry the calculus. The first one is the most usual method used to solve the groundwa- ter equation, a scheme of seven blocks on finite difference, used by MODFLOW and the most popular programs to analyze aquifers. The second one is the finite difference equation proposed by Li, et al (2014) using a scheme of 19 blocks. And the last one is a finite difference scheme of 19 blocks we have proposed to attack the anisotropy of the mediums.

## Example

Some examples of anisotropic flow are provided in **ex/** folder.

## Installation

You have to have installed:
* **Fortran** and **C** compiler
* **BLAS** libraries
* **LACPACK** libraries
* **MPICH** libraries
* **PETSc 3.7** libraries
* **HDF5** parallel libraries connected with PETSc libraries (optional and recommend)

Then, just compile the source files in **src/** folder with **make ANISOFLOW**

## Tests

To test ANISOFLOW just run the bash files on the examples folders:
* **ex/Blesset/Run.sh**

## Contributors

Santiago Ospina De Los Ríos, Kevin Alverto Pérez, Oscar David Álvarez-Villa y Andrea Gómez

## License

ANISOFLOW is Copyright, 2015-2016 **National University of Colombia**, **Gotta Ingeniería S.A.S.**, and the **ANISOFLOW Development Team**, licensed under terms of the **GNU General Public License (GPL)**. This includes all software, documentation, and associated materials.

ANISOFLOW is free software, you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation.
