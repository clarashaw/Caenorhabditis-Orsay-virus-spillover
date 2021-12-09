# Caenorhabditis-Orsay-virus-spillover

This repository contains the code and data used in the manuscript "Developing an empirical model for spillover and emergence: Orsay virus host range in Caenorhabditis" submitted to Proceedings B.

The R code "TitrationResults.R" should be used with data file "2021_11_29_MJOrsayTM_FCK_FCH_GAA.manip.csv" to determine the TCID50 of our viral filtrate.
The R code "PhyloMixedModel.R" should be used with data file "spillover.blocks.csv" for the analysis of our spillover data.
The R code "PhyloMixedModelTransmission.R" should be used with the data file "TransmissionData.csv" for the analysis of transmission following spillover.

Both "PhyloMixedModel.R" and "PhyloMixedModelTransmission.R" use a tree file ("caenorhabditis_35species_1167orthos.treefile") from Stevens et al. 2020 Current Biology obtained from https://zenodo.org/record/3571457#.YJvu76hKg2x
