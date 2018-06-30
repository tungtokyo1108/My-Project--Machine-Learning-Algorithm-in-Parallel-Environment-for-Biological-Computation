This package implements a scalable, parallel implementation of the fast algorithm for fitting 
a model of biological computation on large data set. 

Note: All files have name like "....VI.h/cpp" which we have updated new optimization algorithm for computating the parameters of model

Preparing publication
------------------------------------------------------------------------------------------------------------------------------------------
Title: Fast Computational Software for Biological Database via Stochastic Variational Algorithm in Parallel Environment 

Authors: Tung Dang and Hirohisa Kishino

Preferred journal: https://academic.oup.com/mbe 

Our presentation (included guides of software): https://drive.google.com/file/d/16bKeL7sQQtv82Tw17euQ-y6lZP9MlHzx/view 

Abstract 
-----------------------------------------------------------------------------------------------------------------------------------------
The pattern of molecular evolution varies among gene sites and genes in a genome. By taking into account the complex heterogeneity of evolutionary processes among sites in a genome, Bayesian infinite mixturemodels of genomic evolution enable robust phylogenetic inference. With large modern data sets, however,the computational burden of Markov chain Monte Carlo sampling techniques becomes prohibitive. Here, we have developed a variational Bayesian procedure to speed up the widely used PhyloBayes MPI program, which deals with the heterogeneity of amino acid propensity. Rather than sampling from the posterior distribution, the procedure approximates the (unknown) posterior distribution using a manageable distribution called the variational distribution. The parameters in the variational distribution are estimated by minimizing Kullback-Leibler divergence. To examine performance, we analyzed three large data sets consisting of mitochondrial, plastid-encoded, and nuclear proteins. Our variational method accurately approximated the Bayesian phylogenetic tree, mixture proportions, and the amino acid propensity of each component of the mixture while using orders of magnitude less computational time.
