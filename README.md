# Societal risk-benefit trade-offs of potential adverse events following immunisation with COVID-19 vaccines

## Authors

Frank Sandmann,1,2 Nicholas G. Davies,1 Edwin van Leeuwen,1,2 William Waites,1,3 Centre for the Mathematical Modelling of Infectious Diseases COVID-19 working group, Peter White,2,4 Mary Ramsay,5 Anna Vassall,1 W John Edmunds,1 Mark Jit1,6


1 Centre for Mathematical Modelling of Infectious Diseases, London School of Hygiene and Tropical Medicine, London, UK

2 Statistics, Modelling and Economics Department, National Infection Service, Public Health England, London, UK

3 School of Informatics, University of Edinburgh, Edinburgh, UK

4 MRC Centre for Global Infectious Disease Analysis and NIHR Health Protection Research Unit in Modelling and Health Economics, Imperial College London, London, UK

5 Immunisation and Countermeasures Department, National Infection Service, Public Health England, London, UK

6 School of Public Health, University of Hong Kong, Hong Kong SAR, China


## Summary of project

This repository contains the code required to explore whether the benefits of the COVID-19 vaccine may outweigh the harm from a wider societal perspective, taking into account mild disease, transmission and the economic impact of both adverse events and COVID-19. 

We extended a previously developed age-structured epidemiological-economic model of COVID-19 mass vaccination in the UK (https://doi.org/10.1016/S1473-3099(21)00079-7) to explore the impact of the adverse events following immunisation (AEFIs) observed with the AZD1222 (Oxford/AstraZeneca) vaccine. 

In a subsequent threshold analysis, we used a wide range of hypothetical combinations of rates and severities of AEFIs. 

We traded-off the risk of AEFIs in vaccinated individuals with the community-wide benefit of vaccination over ten years. 

All outcomes were expressed in terms of quality-adjusted life years, QALYs, the net health value (defined as the difference in QALYs and costs, with costs converted into QALYs at £20,000/QALY), and the benefit-risk ratio (of the vaccine over the harm from AEFIs).


## Overview of files

The main script of the analysis runs 3 vaccination scenarios, 3 physical distancing scenarios, and the threshold analysis of 7x14 scenarios of severity and risk of adverse events. On a current laptop (Core i7 2.8GHz CPU, 32GB RAM) this takes around 2-3 hours to run. The analysis requires covidM to be installed (originally developed by Nick Davies et al.). The version used here can be found in the folder "covidm_for_fitting". The file path in the main script will need to be updated to the directory of this project where it was saved.

The scripts contain the following:

- "0_AEFI-COVID19-vacc_main.R" : The main file of the analysis to run the simulations, plot visualisations, save results.

- "1_AEFI-COVID19-vacc_funs.R" : Custom-made functions for the analysis, including wrapper functions using covidM and the cost-effectiveness analysis.

- "2_AEFI-COVID19-vacc_dataUK.R" : Includes the required input data used for the UK; all open access and in the public domain.

- The file "nationallifetables3yearuk.xlsx" contains publicly available data of life tables for the UK from ONS (available from https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/lifeexpectancies/bulletins/nationallifetablesunitedkingdom/2017to2019 ; updated in September every year)


## Linked publications

The paper can be found at:
[url to be added]


## Funding statement

...

The authors had sole responsibility for the study design, data collection, data analysis, data interpretation, and writing. 

The views expressed in this publication are those of the authors and not necessarily those of EC, MRC, NHS, NIHR, PHE or the UK Department of Health and Social Care.
