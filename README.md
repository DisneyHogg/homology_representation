# Action of Riemann surface automorphisms on homology.
This repository contains Sage code developed for the paper [Explicit Homology Representation for Finite Groups Acting on Riemann Surfaces](), written by myself and S. Allen Broughton. 

Given a generating vector for a group $G$ acting on a Riemann surface $S$ with quotient genus 0, the function `homology_representation` computes the representation of $G$ acting on $H_1(S; R)$ for appropriate rings $R$. To use the function, download `hom_reps.py` and in Sage run `load("/.../hom_reps.py")`.

## Requirements
I shall not give a complete list of system requirements, but as a rough guideline:
* the bulk of the code was written in [Sage](https://www.sagemath.org/) version 9.4,
* notebooks were written in Jupyter. 

