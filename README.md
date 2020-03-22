This GitHub repo gathers all the input files and scripts related to our study entitled "Phylogeographic and phylodynamic approaches to epidemiological hypothesis testing" (now on [biorXiv](https://www.biorxiv.org/content/10.1101/788059v1)): BEAST XML files of the continuous phylogeographic and skygrid-GLM analyses, as well as R scripts and related files needed to run all the landscape phylogeographic testing analyses. Continuous phylogeographic and phylodynamic (skygrid-GLM) inferences were performed with the Bayesian methods implemented in the open-source program [BEAST](http://github.com/beast-dev/beast-mcmc). Subsequent dispersal statistics estimation and landscape phylogeographic analyses were implemented and performed with R functions available in the open-source package "[seraphim](http://evolve.zoo.ox.ac.uk/Evolve/Seraphim.html)".

## System requirements

- BEAST: the program runs on any operating system and requires the installation of the BEAGLE library for fast computations (https://github.com/beagle-dev/beagle-lib). Version of BEAST used in the present study: 1.10.4
- "seraphim": the R package was mainly tested on Unix operating systems and requires the preliminary installation of several other R packages listed in the manual provided with the package

## Installation guide

- instructions to install BEAST: http://beast.community/installing
- instructions to install BEAGLE: https://github.com/beagle-dev/beagle-lib
- instructions to install "seraphim": https://github.com/sdellicour/seraphim

(Typical install time on a "normal" desktop computer: very variable)

## Demo and instructions for use

The file `R_script_analyses.r` gathers all the R scripts used to performed post hoc landscape phylogeographic analyses. In addition, there are also several existing tutorial allowing to perform the same analyses:
- instructions to run a continuous phylogeographic analysis in BEAST: tutorial available [here](https://beast.community/workshop_continuous_diffusion_yfv). Output: time-scaled annotated phylogenetic trees. Run time: several weeks, depending on the computational resources.
- instructions to estimate dispersal statistics with "seraphim": tutorial available [here](https://github.com/sdellicour/seraphim/blob/master/tutorials/Estimating_dispersal_statistics.pdf) (related example files are available in the directory with the same name). Output: estimated dispersal statistics (graphics, mean/median values, HPD intervals). Run time: < 1 hour.
- instructions to use "seraphim" to perform relaxed random walk (RRW) simulations along posterior trees: tutorial available [here](https://github.com/sdellicour/seraphim/blob/master/tutorials/RRW_simulations_along_trees.pdf) (related example files are available in the directory with the same name). Output: spatially-annotated posterior trees along which a new RRW diffusion process has been re-simulated. These simulations are used to generate a null dispersal model, which is further exploited to test the impact of environmental factors on the dispersal location and velocity of viral lineages, as well as to test the impact of migratory bird flyways on the dispersal frequency of viral lineages. Run time: several hours when ran in parallel on several CPUs.
- instructions to test the impact of environmental factors on lineage dispersal locations with "seraphim": tutorial available [here](https://github.com/sdellicour/seraphim/blob/master/tutorials/Impact_on_dispersal_direction.pdf) (related example files are available in the directory with the same name). Output: distribution of *E* statistics (see the text for further detail) and statistical support (Bayes factors) for each environmental factor. Run time: < 1 day.
- instructions to test the impact of environmental factors on lineage dispersal velocity with "seraphim": tutorial available [here](https://github.com/sdellicour/seraphim/blob/master/tutorials/Impact_on_dispersal_velocity.pdf) (related example files are available in the directory with the same name). Output: distribution of *Q* statistics (see the text and Appendix for further detail) and statistical support (Bayes factors) for each environmental factor. Run time: several days when ran in parallel on several CPUs.

