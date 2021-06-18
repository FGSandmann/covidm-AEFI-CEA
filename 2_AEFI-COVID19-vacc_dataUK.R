# Value of COVID-19 vaccine - AEFI
# UK data-file

################################
## set-up
################################

# data wrangling for the econ analysis
library(dplyr)


################################
## data for epi model (covidm)
################################

# include health burden process as in covidm for consistency
# Health burden processes
probs = fread(
  "Age,Prop_symptomatic,IFR,Prop_inf_hosp,Prop_inf_critical,Prop_critical_fatal,Prop_noncritical_fatal,Prop_symp_hospitalised,Prop_hospitalised_critical
10,0.66,8.59E-05,0.002361009,6.44E-05,0.5,0,0,0.3
20,0.66,0.000122561,0.003370421,9.19E-05,0.5,9.47E-04,0.007615301,0.3
30,0.66,0.000382331,0.010514103,0.000286748,0.5,0.001005803,0.008086654,0.3
40,0.66,0.000851765,0.023423527,0.000638823,0.5,0.001231579,0.009901895,0.3
50,0.66,0.001489873,0.0394717,0.001117404,0.5,0.002305449,0.018535807,0.3
60,0.66,0.006933589,0.098113786,0.005200192,0.5,0.006754596,0.054306954,0.3
70,0.66,0.022120421,0.224965092,0.016590316,0.5,0.018720727,0.150514645,0.3
80,0.66,0.059223786,0.362002579,0.04441784,0.5,0.041408882,0.332927412,0.3
100,0.66,0.087585558,0.437927788,0.065689168,0.5,0.076818182,0.617618182,0.3")


reformat = function(P)
{
  # 70-74,3388.488  75-79,2442.147  80-84,1736.567  85-89,1077.555  90-94,490.577  95-99,130.083  100+,15.834
  x = c(P[1:7], weighted.mean(c(P[8], P[9]), c(3388.488 + 2442.147, 1736.567 + 1077.555 + 490.577 + 130.083 + 15.834)));
  return (rep(x, each = 2))
}

## alternative data for UK:
##Prop_critical_fatal assumed 32% (based on 20k pts https://www.bmj.com/content/369/bmj.m1985/) 
##Prop_hospitalised_critical assumed 17% (based on 20k pts https://www.bmj.com/content/369/bmj.m1985/) 
probs[, "Prop_critical_fatal"]        <- 0.32
probs[, "Prop_hospitalised_critical"] <- 0.17

# Ensemble estimates, Table S3, p. 22/23 https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2918-0/MediaObjects/41586_2020_2918_MOESM1_ESM.pdf
P.death_nonicu <- c(0.00003,
                    0.00001, 0.00001, 
                    0.00003, 0.00006, 
                    0.00013, 0.00024, 0.00040, 0.00075, 
                    0.00121, 0.00207, 0.00323, 0.00456, 
                    0.01075, 0.01674, 
                    0.03203, 
                    0.08292)

## ONS 75-79 vs 75+: https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/bulletins/annualmidyearpopulationestimates/mid2019estimates
P.death_nonicu[16] <- ((2325296/5687895) * P.death_nonicu[16]   + (1-2325296/5687895) * P.death_nonicu[17])
P.death_nonicu <- P.death_nonicu[1:16]

# instead of this assumption, use slightly higher IFRs informed by REACT3, Table 3 
#P.death_nonicu[14:15] <- 0.0313 # if using this, overshooting observed deaths; rather only change 75+ to be conservative
P.death_nonicu[16] <- 0.116

# update prop symptomatic with ratio admisssion/cases over deaths/cases in UK on 15.07.2020
#probs[, "Prop_symp_hospitalised"] <- probs[, Prop_symp_hospitalised]/(18.5/2.3)*(44.7/15.4) 
#probs[, "Prop_symp_hospitalised"] <- P.death_nonicu * (44.7/15.4) 

P.symp_hospitalised <- c(0.010, 0.002, 0.002, 0.003, 
                         0.005, 0.011, 0.015, 0.020, 
                         0.027, 0.044, 0.062, 0.073, 
                         0.080, 0.084, 0.109, 0.454)

P.icu_symp     = P.symp_hospitalised * reformat(probs[, Prop_hospitalised_critical]);
P.nonicu_symp  = P.symp_hospitalised * reformat(probs[, (1 - Prop_hospitalised_critical)]);

#Ip = pre-symptomatic (100% infectious)
#Is = symptomatic (100% infectious)
#Ia = symptomatic/subclinical (50% infectious)
# report = "ipo" means report the incidence, prevalence, and "outcidence" of the subprocesses
burden_processes_UK = list(
  
  list(source = "Ip", type = "multinomial", names = c("to_icu", "to_nonicu", "to_nonhosp"), report = c("", "", ""),
       prob   = matrix(c(P.icu_symp, P.nonicu_symp, 1 - P.icu_symp - P.nonicu_symp), nrow = 3, ncol = 16, byrow = T),
       delays = matrix(c(cm_delay_gamma(7, 7, 60, 0.25)$p, cm_delay_gamma(7, 7, 60, 0.25)$p, cm_delay_skip(60, 0.25)$p),
                       nrow = 3, byrow = T)),
  
  list(source = "to_icu", type = "multinomial", names = "icu", report = "pi",
       prob   = matrix(1, nrow = 1, ncol = 16, byrow = T),
       delays = matrix(cm_delay_gamma(10, 10, 60, 0.25)$p, nrow = 1, byrow = T)),
  
  list(source = "to_nonicu", type = "multinomial", names = "nonicu", report = "pi",
       prob   = matrix(1, nrow = 1, ncol = 16, byrow = T),
       delays = matrix(cm_delay_gamma(8, 8, 60, 0.25)$p, nrow = 1, byrow = T)),
  
  list(source = "new_EEa", type = "multinomial", names = c("new_infections"), report = c("i"),
       prob   = matrix(1, nrow = 1, ncol = 16, byrow = T),
       delays = matrix(cm_delay_skip(60, 0.25)$p, nrow = 1)),
  
  list(source = "E", type = "multinomial", names = c("death", "null"), report = c("o", ""),
       prob   = matrix(c(P.death_nonicu, 1 - P.death_nonicu), nrow = 2, ncol = 16, byrow = T),
       delays = matrix(c(cm_delay_gamma(22, 22, 60, 0.25)$p, cm_delay_skip(60, 0.25)$p), nrow = 2, byrow = T))
)


##############################
## B) Data for econ model
##############################

# COSTS:
# HOSPITALISATION
#- non-ICU hospitalisations: NHS reference costs 2018/19 (latest data), using ICD-10 codes for other viral pneumonia: J12.8 (Other viral pneumonia), J12.9 (Viral pneumonia, unspecified). The base HRG code per hospitalisation for other viral pneumonia is DZ11 (PD14 for age <=18 years).
#- ICU hospitalisations: NHS reference costs 2018/19 (latest data), Adult Critical Care (activity-weighted HRG codes XC01Z- XB07Z) and Paediatric Critical Care (activity-weighted HRG codes XB01Z- XB07Z). Using prevalence to account for bed-day values, and assuming 10 days in ICU in line with previous paper (Davies et al).
#- alternative costs: hospitalisation from 2009 flu pandemic

# PRIMARY CARE
#- GP visits (unit costs plus proxy of cases based on flusurvey)
#- helpline calls (unit costs plus proxy of cases based on flusurvey)

# VACCINATION
#- administration (assumed service payment for flu)
#- AEs (assumed costs equivalent to 1 GP visit)
#- vaccine based on cost price (cf news reports)
#- government funding for R&D


##############################
# costs NHS perspective (hospitalisation, vaccination, primary care)
##############################

##############################
# costs hospitalisation
##############################
# NHS reference costs 2018/19 using ICD-10/HRG codes for other viral pneumonia
# ICD	J128	Other viral pneumonia	        DZ11	a	PD14	p	PM45	PM_Neut_Cancer
# ICD	J129	Viral pneumonia, unspecified	DZ11	a	PD14	p	PM45	PM_Neut_Cancer
# a Base HRG; p Age 18 years and under

# act-weighted hosp episode adults: 1770.38; <=18y: 1482.60
# act-weighted ICU/day adults: 1504.47; <=18y: 2246.43
# assumed 10 days ICU above in subprocess (cf previous paper Nick et al)

cost_hosp_m <- c(rep(1482.60,3), sum(1482.60/5*4, 1770.38/5*1), rep(1770.38,12)) %>%
  setNames(c(params$pop[[1]]$group_names))

cost_icu_m <- c(rep(2246.43,3), sum(2246.43/5*4, 1504.47/5*1), rep(1504.47,12)) %>%
  setNames(c(params$pop[[1]]$group_names))

##############################
# alternative costs hospitalisation (2009 costs pandemic flu)
##############################

# Lau, Hauck et al 2019: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6491983/
#data.frame("group"  = c("0-4", "5-14", "15-24", "25-44", "45-64", "65+"), 
#           "mean" = c(2123.30, 2293.38, 2237.61, 2191.37, 2785.67, 3299.10))
cost_hosp2_m <- c(2123.30, rep(2293.38,2), rep(2237.61,2), rep(2191.37,4), rep(2785.67,4), rep(3299.10,3)) %>%
  setNames(c(params$pop[[1]]$group_names))

# HS index for inflating costs from 2011 to 2018/19 (to align with ref costs)
infl_HSidx <- c("2010"=3.0, "2011"=2.1,  "2012"=1.7,  "2013"=1.1,  "2014"=0.9, 
                "2015"=1.3, "2016"=0.35, "2017"=2.12, "2018"=1.16, "2019"=2.31)

cost_hosp2_m <- cost_hosp2_m * (1+sum(infl_HSidx[which(names(infl_HSidx) %in% c(2011:2019))])/100)



##############################
# costs PPE in hospital
##############################

#previous study on MERS-CoV in 2015 estimated the additional costs on 
#enhanced PPE equipment (mask, gown, gloves, goggles) at £2.50 per patient visit, 
#with six patient visits per day. Accounting for the additional time of an estimated 
#15 minutes to put on and take off the PPE as well as disposal plus documentation 
#per patient visit, at 6 visits per patient per day came at additional £29.50 for nurses 
#and £45 for physicians, and the total costs per patient at £119 (uprated to 2019 value)
#[https://www.sciencedirect.com/science/article/pii/S0195670116305813#tbl2fna]. 

# assume eg 1:8 nursing ratio
#https://www.nice.org.uk/guidance/sg1/resources/safe-staffing-for-nursing-in-adult-inpatient-wards-in-acute-hospitals-61918998469

cost_PPE_m <- 1/8 * 119 * (1+sum(infl_HSidx[which(names(infl_HSidx) %in% c(2015:2019))])/100)
#cost_PPE_m <- rep(119* (1+sum(infl_HSidx[which(names(infl_HSidx) %in% c(2015:2019))])/100), 
#                  length(params$pop[[1]]$group_names)) %>%
#  setNames(c(params$pop[[1]]$group_names))



##############################
## costs calls NHS111 and GP visits
##############################

# PHE synromic surveillance data concentrates on trends; 
#gives "total nr of calls", 
#percentage of "potential COVID19 calls" in total (10-40% due to covid-19-like symptoms)
#pot covid19 calls by age as % of daily calls per age (which are not presented....)

#However, while "potential covid19 calls" in total at 5% for weeks now,
#flusurvey calls in total (for n=~4k) at about 100 per 1,000 participants with fever/cough (so at least 2x as high)

# instead of 5% and 10% for GP vists and NHS111 calls use 0.5% and 1% to be conservative..
# may be worth trying to get better data from syndromic surveillance? 
# similar to Kawasaki disease, which is not rated much differently by kids but possibly by parents https://link.springer.com/article/10.1007/s00431-017-2937-5
# may not needed to be included..?

#library(digitize)
#df_flusurvey <- digitize(paste0(cm_path, "Flusurvey_June2020.png"))
# flusurvey data for 15 weeks from week 9-23, England, rate of HC contact per 1,000 participants
#"A total of 4,103 participants completed the weekly COVID-19 surveillance survey in week 23, of which 134 (3.3%) reported fever or cough"
df_flusurvey <- data.frame("GPvisit" = c(107.90, 93.47, 99.65, 47.42, 39.17,
                                         38.48, 43.29, 49.48, 48.79, 40.54, 
                                         47.42, 51.54, 46.73, 35.05, 44.67),
                           "GPcalls" = c(43.67, 48.85, 62.06, 85.63, 87.93,
                                         99.42, 104.02, 117.81, 144.25, 178.73,
                                         181.03, 174.13, 170.68, 137.93, 170.68),
                           "NHS111" = c(34.95, 37.24, 116.90, 174.21, 171.34, 
                                        172.49, 135.24, 140.40, 140.40, 149.57, 
                                        118.62, 154.15, 117.47, 101.43, 96.27)
)

colMeans((df_flusurvey/10))
# use only GP visits and NHS111 to avoid overlap; also, values seem high (eg twice as high as compared to NHS111 official trends),
# so could go with 0.5% and 1% (probably also age-dependent)?


##############################
# GP visits
##############################

prob_GPvisits_per_case = 0.05

cost_GP_m <- rep(39 * prob_GPvisits_per_case, length(params$pop[[1]]$group_names)) %>%
  setNames(c(params$pop[[1]]$group_names))


##############################
# NHS111 calls
##############################
# costs for helpline calls (Turner et al., 2012 on NHS111 costs from 4 pilot sites of mean GBP8 per call; the estimate of 12.26 is also reported by them and reflecting the actual costs); cost values at 2011 levels
# inflate to 2019

prob_calls_per_case = 0.1

cost_calls_m <- rep(12.26 * 
                      (1+sum(infl_HSidx[which(names(infl_HSidx) %in% c(2011:2019))])/100) * 
                      prob_calls_per_case, 
                    length(params$pop[[1]]$group_names)) %>%
  setNames(c(params$pop[[1]]$group_names))


##############################
## costs adverse events following immunisation (AEFI)
##############################

# frequent assumption to use costs of 1 GP visit for AE; para 2.3.6.3 
# https://www.nice.org.uk/guidance/ng103/evidence/economic-modelling-report-pdf-6532084909
# unit costs GP visit in 2019: 39£

pr_vaccAEFI <- 2*0.1 # assumption of 10%; somewhat similar to seasonal flu

cost_vaccAEFI_m <- rep(39 * pr_vaccAEFI, length(params$pop[[1]]$group_names)) %>%
  setNames(c(params$pop[[1]]$group_names))


##############################
## costs of vaccine administration
##############################

# Item of Service fee of £12.58 per COVID vaccine dose in 2020; assuming administered via GP, and not accounting for costs of GP visit (eg as done as part of other visit)
# https://www.england.nhs.uk/coronavirus/wp-content/uploads/sites/52/2020/03/C0856_COVID-19-vaccineletter_9-Novrevb.pdf#page=5
# could be half if administered through community systems

cost_vaccAdmin_m <- rep(2*12.58, length(params$pop[[1]]$group_names)) %>%
  setNames(c(params$pop[[1]]$group_names))


##############################
## costs of vaccination R&D (public gov funds)
##############################

#-government funding for R&D
# "In total, Sharma reported in a daily Downing Street press conference that the UK government has committed £250 million pounds toward the development of a COVID-19 vaccine."
#https://www.gov.uk/government/news/uk-pledges-544-million-to-find-coronavirus-vaccine
#(47000000+84000000)
cost_vaccRD <- 250000000


##############################
## costs of vaccination
##############################

# explored as part of sensitivity analysis (threshold pricing)


##############################
## costs of lockdown (GDP loss)
##############################
#https://www.ons.gov.uk/economy/grossdomesticproductgdp/timeseries/abmi/qna

GDP_UK_2019Q4 <- 523917000000
GDP_UK_2019daily <- GDP_UK_2019Q4/as.numeric(as.Date("2019-12-31") - as.Date("2019-10-01"))

c_PD_min_trigger <- 1000
c_PD_min         <- GDP_UK_2019daily*0.02
c_PD_lockdown    <- GDP_UK_2019daily*0.05
c_lockdown_high  <- GDP_UK_2019daily*0.1#0.2


##############################
# QALYs lost (NHS perspective)
##############################

#-QALY loss per symptomatic case: 0.008 (based on ILI for 2009 H1N1 pandemic influenza; Van Hoek et al. 2011) 
#-QALY loss per non-fatal hospitalisation: 0.0201 (based on COVID-19; Halpin et al. 2020) #0.018 (based on seasonal flu; Baguelin et al. 2015)
#-QALY loss per non-fatal ICU: 0.15 (based on ICU value)
#-QALY loss per long COVID: 0.034 (based on the ratio of DALY weights applied to the value of an ILI case)
#-QALY loss per fatality: based on most recent life expectancy for UK from ONS as 3-year average over 2016-2018, and age-/sex-specific QALY norms based on EQ-5D-3L for the UK (Ara and Brazier, 2010)
#-QALY loss per vaccination (adverse events): assumed 1 QALD at a chance of 10%; roughly following another study on flu vaccination but could be updated later when more info available https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5302117/
#Discounting beyond year 1 at 3.5%.


#-------------------------------
# QALYs lost per symptomatic episode (assume ILI-values; alternatively use utilities and days..?)
#-------------------------------
# https://www.bmj.com/content/345/bmj.e4445

QALY_ILI_m <- rep(0.008, length(params$pop[[1]]$group_names)) %>%
  setNames(c(params$pop[[1]]$group_names))

#-------------------------------
# QALYs lost per symptomatic episode (assume ILI-values; alternatively use utilities and days..?)
#-------------------------------
#disability weights for post-acute consequences vs moderate community cases (0.219/0.051=4.29),{Wyper, 2020 #129} 
#multiplied by assumed QALYs lost per symptomatic case (4.29*0.008=0.034 QALYs), and assuming 
#a proportion of 10% of cases experiencing long COVID symptoms.{Greenhalgh, 2020 #128}

QALY_longCOVID_m <- rep(0.219/0.051*0.008*0.1, length(params$pop[[1]]$group_names)) %>%
  setNames(c(params$pop[[1]]$group_names))

QALY_ILI_m = QALY_ILI_m+QALY_longCOVID_m

#-------------------------------
# QALYs lost per hospitalisation
#-------------------------------

# Halpin et al., (2020); COVID-19  hospitalisations https://doi.org/10.1002/jmv.26368
#(68*(0.061*(6.5+71)/365.25)+32*(0.155*(12+71)/365.25))/(68+32)

QALY_hosp_m <- rep(0.0201, length(params$pop[[1]]$group_names)) %>%
  setNames(c(params$pop[[1]]$group_names))


# explore in sensitivity analysis QoL decrement for ICU treatment
# in UK at least 0.10, but if extending beyond 1 year could be higher (older study difference of 0.15)
# (LE may also be lower than for others who are healty, so lifetables would skew estimates..)
# https://ccforum.biomedcentral.com/articles/10.1186/cc12745 0.10 (12 mo.)
# https://ccforum.biomedcentral.com/articles/10.1186/cc8848#Tab1 0.15 Y1, 0.11 Y2.5, 0.14 Y5 assume 0.1*5=0.5 (undisc)

QALY_ICU_m <- rep(0.15, length(params$pop[[1]]$group_names)) %>%
  setNames(c(params$pop[[1]]$group_names))

#-------------------------------
# QALYs lost per vaccination (adverse events)
# assumed 1 QALD lost per vaccinee at a risk of 10%; cf. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5302117/
#-------------------------------

QALY_vaccAEFI_m <- rep(1/365.25*2*0.1, length(params$pop[[1]]$group_names)) %>%
  setNames(c(params$pop[[1]]$group_names))

#-------------------------------
# QALYs lost per premature death
#-------------------------------

##### optional adjustment for changes in mortality and comorbidity (Briggs et al)
# adjust standardised mortality rate in life-expectancy
SMR <- 1.25 # taking reports of 20-25% excess deaths

# get most recent life expectancy for UK from ONS as 3-year average over 2017-2019
## https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/lifeexpectancies/datasets/nationallifetablesunitedkingdomreferencetables
LE_UK <- readxl::read_excel(paste0(econ_path, "nationallifetables3yearuk.xlsx"), 
                          sheet="2017-2019", range = "A7:F108") %>%
  setNames(c("age", "mx_M", "qx_M", "lx_M", "dx_M", "ex_M")) %>%
  bind_cols(readxl::read_excel(paste0(econ_path, "nationallifetables3yearuk.xlsx"), 
                              sheet="2017-2019", range = "H7:L108") %>%
              setNames(c("mx_F", "qx_F", "lx_F", "dx_F", "ex_F")) )

# adjust by standardised mortality rate (SMR)
LE_UK$adjSurvM <- c(100000, c(LE_UK$lx_M * exp(-LE_UK$qx_M * SMR))[-101])
LE_UK$adjSurvF <- c(100000, c(LE_UK$lx_F * exp(-LE_UK$qx_F * SMR))[-101])

suppressWarnings({
  LE_UK  <- LE_UK %>% 
    mutate(YrsAliveM = c(c((adjSurvM+lag(adjSurvM) )/2)[-1],
                         last(c((adjSurvM+lag(adjSurvM) )/2)[-1]) ),
           YrsAliveF = c(c((adjSurvF+lag(adjSurvF) )/2)[-1],
                         last(c((adjSurvF+lag(adjSurvF) )/2)[-1]) ) )
  
  LE_UK$adjLEmale   <- sapply(seq_len(nrow(LE_UK)), 
                              function(x){ sum(LE_UK$YrsAliveM[x:nrow(LE_UK)])/LE_UK$YrsAliveM[x] })
  LE_UK$adjLEfemale <- sapply(seq_len(nrow(LE_UK)), 
                              function(x){ sum(LE_UK$YrsAliveF[x:nrow(LE_UK)])/LE_UK$YrsAliveF[x] })
})


# use age groups
LE_UK <- LE_UK %>%
  dplyr::select(age, male=adjLEmale, female=adjLEfemale) %>%
  dplyr::mutate(age = ifelse(age %in% 0:4, "0-4",
                      ifelse(age %in% 5:9, "5-9",
                      ifelse(age %in% 10:15, "10-14",
                      ifelse(age %in% 15:19, "15-19",
                      ifelse(age %in% 20:24, "20-24",
                      ifelse(age %in% 25:29, "25-29",
                      ifelse(age %in% 30:34, "30-34",
                      ifelse(age %in% 35:39, "35-39",
                      ifelse(age %in% 40:44, "40-44",
                      ifelse(age %in% 45:49, "45-49",
                      ifelse(age %in% 50:54, "50-54",
                      ifelse(age %in% 55:59, "55-59",
                      ifelse(age %in% 60:64, "60-64",
                      ifelse(age %in% 65:69, "65-69",
                      ifelse(age %in% 70:74, "70-74",
                             "75+"))))))))))))))) ) %>%
  dplyr::group_by(age) %>%
  dplyr::summarise(male   = mean(male),
                   female = mean(female)) %>%
  dplyr::mutate(age = factor(age, levels= params$pop[[1]]$group_names)) %>%
  dplyr::arrange(age)

#LE_UK <- as.data.table(LE_UK)

## function for discounting life years (quality-adjusted or not)
discount_LY <- function(target_value, timeframe, comorb_reduc_val, disc_rate, female=TRUE) {
  
  # adjust for round numbers; "ceiling" for decimals
  x <- ifelse(timeframe-floor(timeframe)==0, 
              ceiling(timeframe)+1, 
              ceiling(timeframe))
  
  # create table
  y <- data.frame(Year = seq_len(x), Value = target_value )
  
  # add age-/sex-specific QALY value depending on UK norms of EQ-5D (Ara and Brazier, 2010)
  if(female==TRUE) { sex = 0} else { sex = 1}
  
  y <- y %>% dplyr::mutate(Value = 0.9508566 + 0.0212126*sex - 0.0002587*Year - 0.0000332*Year^2)
  
  # adjust for comorbidities, if wanted to
  if(!is.null(comorb_reduc_val)){ 
    y <- y %>% dplyr::mutate(Value = Value * comorb_reduc_val)
  }
  
  # discount after 1st year
  y <- y %>% dplyr::mutate(Disc_value = Value * (1+disc_rate)^-(Year-1) )
  
  # correct for fraction in last year
  y[x,"Disc_value"] <- (timeframe-floor(timeframe))*y[x,"Value"] * (1+disc_rate)^-y[x-1,"Year"]
  
  #print result
  return( round(sum(y$Disc_value), 2) )
}


# print QALYs
QALY_UK = LE_UK

# discount quality-adjusted life-expectancy (QALE); different utilities by sex
for(j in c("female","male")){
  if(j == "female") {fem=TRUE} else {fem=FALSE}
  
  # discounting
  QALY_UK[[paste0(j, "_QALE_disc", sep="")]] <- 
    QALY_UK[, j] %>% unlist() %>% 
    purrr::map_dbl(function(x) discount_LY(x, x, 0.90, 0.035, female=fem) )
}

# mean QALYs total by age
QALY_UK <- QALY_UK %>%
  transmute(group  = age,
            undisc = (female           +male)    /2, 
            disc   = (female_QALE_disc +male_QALE_disc)/2 )


