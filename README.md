# Action of Riemann surface automorphisms on homology.
This repository contains Sage code developed for the paper [Explicit Homology Representation for Finite Groups Acting on Riemann Surfaces](), written by myself and S. Allen Broughton. 

Given a generating vector for a group $G$ acting on a Riemann surface $S$ with quotient genus 0, the function `homology_representation` computes the representation of $G$ acting on $H_1(S; R)$ for appropriate rings $R$. To use the function, download `hom_reps.py` and in Sage run `load("/.../hom_reps.py")`.

`Tutorial.ipynb` demonstrates how the code may be used. The folder `Companion_Files` contains files used in the tutorial, including the method of BRR, precomputed files listing all topologically inequivalent generating vectors acting in genus $2 \leq \sigma \leq 12$, and the comparison of runtimes referred to in the paper. 

## Requirements
I shall not give a complete list of system requirements, but as a rough guideline:
* the bulk of the code was written in [Sage](https://www.sagemath.org/) version 9.7,
* notebooks were written in Jupyter,
* to run the comparison against the BRR code the additional GAP package `kbmag` will need to be installed in the version of GAP Sage uses.

### Remarks
This code has been tested succesfully on Sage version 9.7 and 9.8, but is known to fail in version 10.7 due to changes in the GAP interface. If there if sufficient demand then this code should be versioned to allow consistency.

## Acknowledgements
Not all the code contained in this repository is written by me. I had no part in writing `polyB.sage`; this was programmed by [Antonino Behn and Anita Rojas](https://github.com/rojas-ani/sage-routines) for the paper [*Adapted hyperbolic polygons and symplectic representations for group actions on Riemann surfaces*](https://doi.org/10.1016/j.jpaa.2012.06.030) by Antonio Behn, Rubí E. Rodríguez and Anita M. Rojas. I have included it here for ease of access, as I use it for comparison.

