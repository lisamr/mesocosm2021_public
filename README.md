# mesocosm2021_public
Code for models and figures analyzing diversity-disease patterns in greenhouse mesocosms

## Explanation of scripts in 'Analysis_scripts' folder

- betabinoom_custom_brms.R: code necessary to run models in `brms` with a beta-binomial likelihood
- four_trts_com_disdensity.R: seperate models for each of the four community assembly treatments; response variable = community-level density of diseased inviduals 
- four_trts_com_disprevalence.R: seperate models for each of the four community assembly treatments; response variable = community-level prevalence of diseased inviduals 
- four_trts_spp_disprevalence.R: seperate models for each of the four community assembly treatments; response variable = species-level prevalence of diseased inviduals 
- host_competence_calculation.R: simple calculations of mean and SD of host competence, estimated from monospecific trays
- splevel_models.R: Model assessing species-level disease prevalence from all trays togehter. See results under section "Drivers of species-level disease prevalence".

Contact lrosenthal@ucdavis.edu for questions.
