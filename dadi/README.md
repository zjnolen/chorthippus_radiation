∂a∂i
================

-   [Introduction](#introduction)
-   [SFS Production](#sfs-production)
-   [Model Fitting](#model-fitting)
-   [Model Selection](#model-selection)
-   [Parameter Uncertainties](#parameter-uncertainties)
-   [References](#references)

Introduction
============

δaδi is a demographic inference tool developed by Gutenkunst et al. (2009) that uses diffusion approximation to fit an estimated SFS to your data SFS. Using an optimization algorithm it explores various parameters in your demographic model, allowing it to find the maximum likelihood parameter values. This maximum likelihood can be used to select the best fitting demographic model from a set of candidate models, and these parameter estimates can be used to infer demographic history.

Fitting demographic models allows multiple hypotheses to be tested. For instance, we could compare a model with and without migration between two populations to see which is more likely. We can also compare migration models for different population pairs to see which have higher rates of migration. The important thing is to first choose models that are appropriate to answer your question. For this manuscript, we use the four simple models: (A) Divergence without gene flow (no\_mig), (B) Divergence with Continous Gene Flow (sym\_mig), (C) Divergence with Ancestral Gene Flow (anc\_sym\_mig), and (D) Divergence with Gene Flow in Secondary Contact, and two more complex ones: (E) Isolation with Migration and Growth (IMG) and (F) Secondary Contact with Growth (SCG).

Each of these models form a hierarchy, with the simpler models being special cases of the more complex ones. For instance, in the case of no\_mig, we can see this as a special case of sym\_mig where the migration rate is zero.

Once a set of models have all been fit, we can test our hypotheses using [model selection](#model-selection) and [parameter estimation](#parameter-uncertainties).

For additional resources, there are also example scripts available in the [repository for δaδi](https://bitbucket.org/gutenkunstlab/dadi/src/master/examples). The analyses in these manuscript were run under version 2.0.5.

SFS Production
--------------

We built 2D SFS for all population level comparisons in the Biguttulus group, along with one to compare the subspecies of the *Pseudochorthippus* genus. The protocol for this can be found in the [SFS section](../sfs) of this repository.

Model Fitting
-------------

To fit our models we used the [dadi_pipeline](https://github.com/dportik/dadi_pipeline) developed by Portik et al. (2017), running three independent optimizations for each model and population comparison to ensure convergence. This pipeline performs several rounds of replicate optimizations, carrying forward the best fit parameters from each round to the next to approach a global optimum.

We implemented this pipeline in the script [`demo_model_run.py`](demo_model_run.py), which runs one dadi pipeline optimization for a given taxa pair and model. The settings for each model are defined in this script and presented in the Molecular Ecology paper this repository is built for.

Our script uses a modified version of the Models_2D.py and Optimize_Functions.py scripts in the dadi pipeline, both of which are in this repository. Modifications are noted in comments containing 'ZN' in the Optimize_Functions.py script. As we used different optimizers for each model (optimize_log_fmin for models without growth and optimize_log for models with growth), we manually changed this in Optimize_Functions.py depending on the models being run. Changing this easily is now built into dadi pipeline.

This script will also create a 'results_summary' file, appending each dadi pipeline run upon completion to easily assess model convergence from one file.

Model Selection
---------------

We performed model selection using the Godambe Information Matrix (GIM) adjusted likelihood ratio test described by Coffman et al. (2016). This method performs a likelihood ratio test and adjusts to account for linkage present in the SFS using bootstraps. Bootstrap production is described in [the sfs section of this repository](../sfs#bootstrapping).

We used the [`GIM_LRT.py`](GIM_LRT.py) script to perform this test, which reads in parameter estimates from our results_summary files.

Parameter Uncertainties
-----------------------

To determine parameter uncertainties, we used the GIM_uncert function in dadi to calculate standard deviation around our parameter estimates. We then converted the parameters into real values, scaling them to mutation rate. We then converted the uncertainties around each parameter into real values using principles of the [propogation of uncertainty](https://en.wikipedia.org/wiki/Propagation_of_uncertainty). The calculations for this are done in the [`calc_params_uncert.py`](calc_params_uncert.py) in this folder, which reads in parameters from our results_summary files.

References
==========

Coffman, Alec J., Ping Hsun Hsieh, Simon Gravel, and Ryan N. Gutenkunst. 2016. “Computationally Efficient Composite Likelihood Statistics for Demographic Inference.” *Molecular Biology and Evolution* 33 (2): 591–93. doi:[10.1093/molbev/msv255](https://doi.org/10.1093/molbev/msv255).

Gutenkunst, Ryan N., Ryan D. Hernandez, Scott H. Williamson, and Carlos D. Bustamante. 2009. “Inferring the Joint Demographic History of Multiple Populations from Multidimensional SNP Frequency Data.” *PLOS Genetics* 5 (10): e1000695. doi:[10.1371/journal.pgen.1000695](https://doi.org/10.1371/journal.pgen.1000695).

Portik, Daniel M., Adam D. Leaché, Danielle Rivera, Michael F. Barej, Marius Burger, Mareike Hirschfeld, Mark-Oliver Rödel, David C. Blackburn, and Matthew K. Fujita. 2017. “Evaluating Mechanisms of Diversification in a Guineo-Congolian Tropical Forest Frog Using Demographic Model Selection.” *Molecular Ecology* 26 (19): 5245–63. doi:[10.1111/mec.14266](https://doi.org/10.1111/mec.14266).
