This package implements a scalable, parallel implementation of the fast algorithm for fitting 
a model of biological computation on large data set. 

Note: All files have name like "....VI.h/cpp" which we have updated new optimization algorithm for computating the parameters of model

Submitted publication
------------------------------------------------------------------------------------------------------------------------------------------
Fast Computational Software for Biological Database via Stochastic Variational Algorithm in Parallel Environment 
Dang Thanh Tung and Hirohisa Kishino

Submitted journal: https://academic.oup.com/mbe 

Abstract 
-----------------------------------------------------------------------------------------------------------------------------------------
Bayesian mixture models for modeling across site variation of the substitution process are now used in a wide variety of applications in phylogenetic reconstruction. Although Monte-Carlo Markov chain (MCMC) sampling techniques make approximate inference possible for both finite and infinite mixture models, the computational burden is prohibitive on the large modern data sets. To overcome this problem, we developed new algorithms for fast and accurate inference of the model underlying the PhyloBayes MPI program approaching variational Bayesian procedures. Variational frameworks convert the problem of approximating posterior distributions into solving a sequence of unconstrained optimization problems. We analyzed empirical large-scale datasets to compare time estimates produced by variational algorithm with those reported by using MCMC approaches in PhyloBayes MPI. We demonstrated that variational methods achieve accuracy competitive with Markov chain Monte Carlo approaches while requiring orders of magnitude less computational time.
