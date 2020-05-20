# Switching_Cov_Coupled
Switching with covariate code. This has a seed set such that we can prove it produces identical results to the original code in the Switching_coupled repository when no covariate data is used. 



### Genstates_coupled:
random number generation is aligned at the very start, and
at the start of each iteration, but get out of sync after the generation
of a set of switching times: although the pregeneration of exponentials
and uniforms is now being done in the same order, one version does a
fixed set of them and the other only does as many as it needs.

The solution, in the _coupled versions of the code (attached), is to use
the state of the RNG before each "group" to decide its state after, by
generating a random seed. Before that is used, we generate the necessary
exponentials and uniforms, in a kind of branching off of the RNG.
