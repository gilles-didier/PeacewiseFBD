# PiecewiseFBD
Probability densities of fossil ages and tree topologies under the skyline FBD model

Software includes
 - 'getPost'
	computes the divergence time distibution associated to a given clade from a set of possible trees, the fossil stratigraphic intervals and provides coda files for MCMC convergence diagnostics
 - 'getAIC'
	computes the divergence time distribution associated to a given set of nodes of a single tree and the fossil stratigraphic intervals
 - 'testAIC'
	computes the extinction time distibutions associated to a given set of clades from a set of possible trees, the fossil stratigraphic intervals and provides coda files for MCMC convergence diagnostics

type

	> make all
	in a console opened on the src directory for compiling standalone software.


Directory "src" contains the C sources of the standalone software

Directory "data" contains the dataset studied in "Testing extinction events and temporal shifts in diversification and fossilization rates":
    - 'Empirical_dataset_tree.newick'
      contains 100 equiparsimonious trees of Cotylosauria
	- 'Empirical_dataset_fossils.csv'
	   contains the corresponding fossil ages
	- 'Simulated_dataset_tree.newick'
      contains the simulated tree
   - 'Simulated_dataset_fossils.csv'
       contains the  corresponding fossil ages


A complete description of the options of the programs is given below.

------------
 getPost 
------------

--------------------------
REQUIREMENT

	'getPost' requires the gsl libraries.

--------------------------
COMPILING

	Just type
	> make getPost
	in a console opened on the directory containing the source files to build the binary.

--------------------------
DESCRIPTION

	'getPost' computes the divergence time distibution associated to a given clade from a set of possible trees, the fossil ages, and the diversification rates


--------------------------
--------------------------
MANUAL
--------------------------
--------------------------


--------------------------

NAME
	getPost - provides coda files for MCMC convergence diagnostics and posterior distribution of the model specification parameters
	
SYNOPSIS
	getPost [OPTIONS] <tree(s)> <fossil ages> <model specification> [output File]

DESCRIPTION
	Provides coda files for MCMC convergence diagnostics and posterior distribution of the model specification parameters of  <model specification> from the tree(s) in file <tree(s)> (in Newick format) with the fossil ranges provided in <fossil ages> (in csv format). The two files <output File>.out and <output File>.ind are to be read by the R package coda.

	Options are
	-w <speciation rate width> <extinction rate width> <fossilization rate width> <sampling probability width>
		set the widths of the sliding windows used for sampling the speciation, extinction and fossilization rates and of the sampling probability during the MCMC
	-i <speciation rate upper bound> <extinction rate  upper bound> <fossilization rate upper bound> <sampling probability upper bound>
		set the upper bounds of the interval in which the speciation, extinction and fossilization rates and the sampling probability are uniformly drawn at the beginning of the MCMC (lower bounds are all 0)
	-b <number>
		set the number of iterations for the burning stage
	-g <number>
		set the thinning parameter (only one iteration in <number> is considered)  
	-f <proportion>
		set the proportion of moves in the parameters in the MCMC proposal (other moves are drawn uniformly amont the fossil ages)
	-a <speciation proportion> <extinction proportion> <fossilization proportion>
		set the relation proportion of moves in the speciation, the extinction and the fossilzation rates (thus in the sampling probability)
	-s <number>
		set the number of samples required
	-r <number>
		set the random seed
	-h
		display help

EXAMPLE
