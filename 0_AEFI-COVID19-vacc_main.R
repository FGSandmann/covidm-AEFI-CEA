# Value of COVID-19 vaccine - AEFI analysis
# original code of "covidm" by ND, CEA code adapted by FGS
# version 2 based on revision for TLID, 21.01.2020

# adaptation for AEFI analysis, 17.03.2021

####################
## 0) set-up
####################

# UPDATE THIS TO FILEPATH OF PROJECT/OTHER SCRIPT FILES
econ_path <- "~/42_Adverse-events/02_Code/" 

# load epi model (covidm)
cm_path = paste0(econ_path, "covidm_for_fitting/")
cm_force_rebuild = FALSE;
cm_build_verbose = TRUE;
cm_version = 2; 
source(paste0(cm_path, "/R/covidm.R"))

# data wrangling for the econ analysis
library(dplyr)

# Build simpler lockdown function
library(stringr)

# source functions for lookdown, epi, econ models
source(paste0(econ_path, "1_AEFI-COVID19-vacc_funs.R"))


##############################
# 1) Data for UK (set-up for national analysis)
##############################

# build parameters for UK over 10 years
params = cm_parameters_SEI3R(cm_uk_locations("UK", 0), 
                             date_start = "2020-01-01", date_end = "2029-12-31",
                             deterministic = TRUE )

# customise
for(i in 1:length(params$pop)){
  
  # basic demography of births, deaths, aging (assuming deaths=births; disabled aging; most recent data from 2019)
  #https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/vitalstatisticspopulationandhealthreferencetables
  b_rate = 713000 / sum( sapply(1:length(params$pop), function(x) sum(params$pop[[x]]$size)))
  params$pop[[i]]$B <- c(b_rate / 365.25, rep(0, 15))
  params$pop[[i]]$D <- rep(b_rate / 365.25, 16)
  #params$pop[[i]]$A <- rep(0 / (5 * 365.25), 16)
  
  # age-dependent ratio of infection:cases (based on Davies et al, Nature paper)
  params$pop[[i]]$y <- c(0.2904047, 0.2904047, 0.2070468, 0.2070468, 0.2676134, 
                         0.2676134, 0.3284704, 0.3284704, 0.3979398, 0.3979398,
                         0.4863355, 0.4863355, 0.6306967, 0.6306967, 0.6906705, 0.6906705)
  
  # susecptibility (based on Davies et al, Nature paper)
  params$pop[[i]]$u <- c(0.3956736, 0.3956736, 0.3815349, 0.3815349, 0.7859512,
                         0.7859512, 0.8585759, 0.8585759, 0.7981468, 0.7981468,
                         0.8166960, 0.8166960, 0.8784811, 0.8784811, 0.7383189, 0.7383189)
  
  # scale u (susceptibility) to achieve desired R0
  current_R0 = cm_calc_R0(params, i); # calculate R0 in population i of params
  target_R0  = 2.7
  params$pop[[i]]$u = params$pop[[i]]$u * target_R0 / current_R0
  
  # natural waning
  params$pop[[i]]$wn <- rep(1/(365.25/52*45), length(params$pop[[i]]$group_names) )
  
  ## Set seeds to control start of outbreak
  params$pop[[i]]$dist_seed_ages = cm_age_coefficients(20, 50, 5 * (0:length(params$pop[[i]]$size))) # infections start in individuals aged 20-50
  
  # 5 new infections each day for 7 days every month
  params$pop[[i]]$seed_times = sapply(30*0:(12* (as.numeric(gsub("-.*$", "", params$time1))-as.numeric(gsub("-.*$", "", params$date0)) )),
                                      function(x) x+rep(0:6, each = 5)) %>% as.vector #rep(0:6, each = 5)) %>% as.vector
  
  # schedule not currently working with observer; moved summer/winter holidays into observer below
}


##############################
# 2) implement lockdown function
##############################

# school holidays to be used with lockdown function (assuming future winter/summer holidays at similar days as in 2020)
school_holidays <- c(sapply(c(seq(as.numeric(gsub("-.*$", "", params$date0)), 
                                  as.numeric(gsub("-.*$", "", params$time1)))), function(x) 
                                    c(as.numeric(as.Date(paste0(x,   c("-07-22"))) - as.Date(params$date0)), 
                                      as.numeric(as.Date(paste0(x,   c("-09-02"))) - as.Date(params$date0)), 
                                      as.numeric(as.Date(paste0(x,   c("-12-20"))) - as.Date(params$date0)), 
                                      as.numeric(as.Date(paste0(x+1, c("-01-03"))) - as.Date(params$date0)) ))) %>%
  setNames(paste0("t_", 1:length(.)))

# To use the observer need to resource the backend with the observer compiled in.
fun_lockdown <- function(lockdown_ini     = 30, 
                         contact_ini      = c(0.9, 0.1,   0.1, 0.1), #home, work, school, and other
                         lockdown_sub     = 20,
                         contact_sub      = c(1,   0.1,   0.1, 0.1),
                         lockdown_sub2    = 50,
                         contact_sub2     = c(1,   0.1,   0.1, 0.1),
                         lockdown_end     = 500/55977178*100000,
                         contact_end      = c(1,   0.67, 0.67, 0.67),
                         contact_holidays = c(1,   0.67,    0, 0.67),
                         t_noFurtherLockdowns = 1e+06,
                         t_hols = school_holidays) { # school holidays (summer/winter for 10 years)
  
  cm_source_backend(
    user_defined = 
      list(model_v2 = list(cpp_observer = 
                             control_lockdown(lockdown_ini, contact_ini,
                                              lockdown_sub, contact_sub,
                                              lockdown_sub2, contact_sub2,
                                              lockdown_end, contact_end,
                                              contact_holidays,
                                              t_noFurtherLockdowns,
                                              t_hols["t_1"],  t_hols["t_2"],  t_hols["t_3"],  t_hols["t_4"],  t_hols["t_5"], 
                                              t_hols["t_6"],  t_hols["t_7"],  t_hols["t_8"],  t_hols["t_9"],  t_hols["t_10"], 
                                              t_hols["t_11"], t_hols["t_12"], t_hols["t_13"], t_hols["t_14"], t_hols["t_15"], 
                                              t_hols["t_16"], t_hols["t_17"], t_hols["t_18"], t_hols["t_19"], t_hols["t_20"], 
                                              t_hols["t_21"], t_hols["t_22"], t_hols["t_23"], t_hols["t_24"], t_hols["t_25"],  
                                              t_hols["t_26"], t_hols["t_27"], t_hols["t_28"], t_hols["t_29"], t_hols["t_30"],  
                                              t_hols["t_31"], t_hols["t_32"], t_hols["t_33"], t_hols["t_34"], t_hols["t_35"],  
                                              t_hols["t_36"], t_hols["t_37"], t_hols["t_38"], t_hols["t_39"], t_hols["t_40"]
                             ))))
}


# source UK input data on costs and QALYs (needs params as input)
source(paste0(econ_path, "2_AEFI-COVID19-vacc_dataUK.R"))


# use updated data for UK
params$processes = burden_processes_UK

# changed to match 130k admissions and ~95% reduction in Rt during initial lockdown
lockdown_trigger <- 30 #30/100000*sum(params$pop[[1]]$size) 


##############################
# change assumptions on costs and QALYs of AEFIs for this analysis
##############################

# function for barplots on aggregated net benefits (absolute and incremental)
econ_fun_barplot <- function(data_set, scenario_name = "A", #placeholder
                             determinist = TRUE, 
                             incl_GDPloss = FALSE, 
                             saveOutcomes = FALSE, nr_iter = 1, 
                             GDP_low = c_PD_min, GDP_medium = c_PD_lockdown, 
                             PD_low_trigger = 1000, 
                             cost_vaccDose_B  = 10, # similar to 2xAZD vaccine
                             cost_vaccDose_C  = 10, # similar to 2xAZD vaccine 
                             cost_Freezers_include = FALSE, #not for AZ vaccine..
                             discB = 0.035, discC = 0.035, lambda_val=20000){
  
  res_df <- econ_mod(data_set, n_agegrs = 16, n_years = 10, n_scens = 3, #faster if specified
                     discount_rate_b  = discB, discount_rate_c = discC,
                     c_PD_min_trigger = PD_low_trigger,
                     c_PD_min         = GDP_low/16,     # daily losses per age group (assumed equal)
                     c_PD_lockdown    = GDP_medium/16,  # daily losses per age group (assumed equal)
                     run_determ       = determinist,
                     cost_vacc_B      = cost_vaccDose_B,
                     cost_vacc_C      = cost_vaccDose_C,
                     costs_lockdown   = incl_GDPloss,
                     add_Freezercost  = cost_Freezers_include,
                     save_by_agegr    = TRUE,         # results for all ages and years
                     n_iter           = nr_iter,      # number of iterations run the model
                     LEs_mort=FALSE)
  
  # sum costs and QALYs for plotting
  res_age <- res_df %>% 
    group_by(group, year, run, scenario, compartment) %>% 
    summarise_if(is.numeric, sum, na.rm=TRUE) %>%
    dplyr::select(-undisc, -t, -value, -undisc2) %>%
    mutate(ScenA = scenario_name) %>%
    dplyr::group_by(group, year, ScenA, scenario, compartment) %>% 
    mutate(disc_MB = ifelse(grepl("QALY", compartment), disc*lambda_val, disc),
           disc_HB = ifelse(grepl("costs", compartment), disc/lambda_val, disc),
           parm    = ifelse(grepl("QALYs", compartment), "health_loss", 
                     ifelse(grepl("costs_lockdown", compartment), "expenses_lockdown", 
                     ifelse(grepl("costs", compartment), "expenses_HC",
                            NA))),
           disc_MB = ifelse(disc_MB==0, NA, disc_MB),
           disc_HB = ifelse(disc_HB==0, NA, disc_HB) ) %>%
    dplyr::select(-disc)
  
  # get NMB and NHB
  res_age2 <- res_age %>%
    dplyr::select(-disc_HB) %>%
    dplyr::group_by(group, year, 
                    ScenA, scenario, run, parm) %>% 
    summarise_if(is.numeric, sum, na.rm=TRUE) %>% 
    tidyr::spread(parm, disc_MB) %>%
    mutate(indicator = "monetary benefit") %>%
    bind_rows(res_age %>%
                dplyr::select(-disc_MB) %>%
                dplyr::group_by(group, year, 
                                ScenA, scenario, run, parm) %>% 
                summarise_if(is.numeric, sum, na.rm=TRUE) %>% 
                tidyr::spread(parm, disc_HB) %>%
                mutate(indicator = "health benefit"))
  
  # include GDP loss proxy?
  if(incl_GDPloss) {
    
    res_age2 <- res_age2 %>%
      mutate(NB_HC       = -health_loss - expenses_HC,     # negative health_loss due to QALY losses!
             NB_lockdown = -health_loss - expenses_lockdown) %>%
      tidyr::gather(parm, disc, NB_HC:NB_lockdown) %>%
      mutate(disc = ifelse(disc==0, NA, disc) )
  } else {
    res_age2 <- res_age2 %>%
      mutate(NB_HC       = -health_loss - expenses_HC) %>% # negative health_loss due to QALY losses!
      tidyr::gather(parm, disc, NB_HC) %>%
      mutate(disc = ifelse(disc==0, NA, disc) )
  }
  
  # sum costs and QALYs for plotting
  res_df <- res_df[, lapply(.SD, sum, na.rm=TRUE), by = .(scenario, run, compartment), .SDcols=c("undisc", "disc") ]
  
  # barplot of monetary/health benefit of parameters
  # all factors into costs/QALYs
  res_df2 <- res_df %>% ungroup() %>%
    mutate(ScenA = scenario_name) %>%
    dplyr::group_by(ScenA, scenario, run, compartment) %>% 
    dplyr::select(-undisc) %>%
    summarise_if(is.numeric, sum, na.rm=TRUE) %>% 
    dplyr::group_by(ScenA, scenario, compartment) %>% 
    mutate(disc_MB = ifelse(grepl("QALY", compartment), disc*lambda_val, disc),
           disc_HB = ifelse(grepl("costs", compartment), disc/lambda_val, disc),
           parm = ifelse(grepl("QALYs", compartment), "health_loss", 
                  ifelse(grepl("costs_lockdown", compartment), "expenses_lockdown", 
                  ifelse(grepl("costs", compartment), "expenses_HC",
                         NA))),
           disc_MB = ifelse(disc_MB==0, NA, disc_MB),
           disc_HB = ifelse(disc_HB==0, NA, disc_HB) ) %>%
    dplyr::select(-disc)
  
  #aggregated NB
  res_df3 <- res_df2 %>%
    dplyr::select(-disc_HB) %>%
    dplyr::group_by(ScenA, scenario, run, parm) %>% 
    summarise_if(is.numeric, sum, na.rm=TRUE) %>% 
    dplyr::group_by(ScenA, scenario, run, parm) %>% 
    tidyr::spread(parm, disc_MB) %>%
    mutate(indicator = "monetary benefit") %>%
    bind_rows(res_df2 %>%
                dplyr::select(-disc_MB) %>%
                dplyr::group_by(ScenA, scenario, run, parm) %>% 
                summarise_if(is.numeric, sum, na.rm=TRUE) %>% 
                dplyr::group_by(ScenA, scenario, run, parm) %>% 
                tidyr::spread(parm, disc_HB) %>%
                mutate(indicator = "health benefit"))
  
  # include GDP loss proxy?
  if(incl_GDPloss) {
    
    res_df3 <- res_df3 %>%
      mutate(NB_HC       = -health_loss - expenses_HC,     # negative health_loss due to QALY losses!
             NB_lockdown = -health_loss - expenses_lockdown) %>%
      tidyr::gather(parm, disc, NB_HC:NB_lockdown) %>%
      mutate(disc = ifelse(disc==0, NA, disc) )
  } else {
    res_df3 <- res_df3 %>%
      mutate(NB_HC       = -health_loss - expenses_HC) %>% # negative health_loss due to QALY losses!
      tidyr::gather(parm, disc, NB_HC) %>%
      mutate(disc = ifelse(disc==0, NA, disc) )
  }
  
  # natural outcomes (eg for CE-plane)
  if(saveOutcomes==TRUE) {
    res_df4 <- res_df %>% ungroup() %>%
      dplyr::filter(scenario == "runA") %>%
      dplyr::rename(disc_base = disc, undisc_base = undisc) %>%
      dplyr::select(-scenario) %>%
      full_join(res_df %>% ungroup,
                by=c("run", "compartment")) %>%
      dplyr::mutate(incr_undisc = undisc - undisc_base,
                    incr_disc   = disc   - disc_base) %>%
      dplyr::select(-disc_base, -undisc_base)
    
    return(list(res_df3, res_df4, res_age2, res_age))
  } else {
    return(list(res_df3))
  }
}



# function for vaccination scenarios and epi curves
epi_fun_vaccplot <- function(data_set, scen_name = "A", cap_text = NULL,
                             label_facets = c(runA = "V0: no vaccination", 
                                              runB = "V1: lower VE estimate,\n45-week protection", 
                                              runC = "V2: higher VE estimate,\n3-year protection"),
                             vacc_start_date=as.Date("2020-12-08"),#+28, # start on 08-12-2020, but protection after 28 days
                             define_yaxis_limit=NA){
  
  res_df <- data_set[data_set$compartment == "cases", .(value = sum(value)), by = .(scenario, run, population, t)]
  df_p2  <- data_set[data_set$compartment == "obs0" & data_set$group == "0-4" & data_set$value == 5]
  
  # check whether there was any lockdown; ignore warning message of empty rows
  if(nrow(data_set[data_set$compartment == "obs0" & data_set$group == "0-4" & data_set$value %in% 1:3]) == 0) {
    df_p3 <- data_set[data_set$compartment == "obs0" & data_set$group == "0-4" & data_set$value == 5] %>% mutate(t=-1)
  } else {
    df_p3 <- data_set[data_set$compartment == "obs0" & data_set$group == "0-4" & data_set$value %in% 1:3]
  }
  
  ggplot(res_df %>% mutate(facet_top = "A")) + 
    geom_rect(data = 
                cbind.data.frame(df_p2 %>% 
                                   group_by(scenario) %>%
                                   dplyr::filter(c(TRUE,diff(t)!=1)) %>%
                                   dplyr::rename(xMin = t),
                                 df_p2 %>% 
                                   group_by(scenario) %>%
                                   dplyr::filter(c(diff(t)!=1, TRUE)) %>% ungroup() %>%
                                   dplyr::select(xMax = t)) %>%
                mutate(facet_top = "A"),
              aes(xmin = xMin, xmax = xMax,
                  ymin = -Inf, ymax = Inf), colour=NA, alpha=0.5, fill = 'lightblue') +
    geom_rect(data = 
                cbind.data.frame(df_p3 %>% 
                                   group_by(scenario) %>%
                                   dplyr::filter(c(TRUE,diff(t)!=1)) %>%
                                   dplyr::rename(xMin = t),
                                 df_p3 %>% 
                                   group_by(scenario) %>%
                                   dplyr::filter(c(diff(t)!=1, TRUE)) %>% ungroup() %>%
                                   dplyr::select(xMax = t)) %>%
                mutate(facet_top = "A"),
              aes(xmin = xMin, xmax = xMax,
                  ymin = -Inf, ymax = Inf), colour=NA, alpha=0.2, fill = 'pink') +
    
    geom_vline(data=res_df[!grepl("A|no vaccination", res_df$scenario), ] %>%
                 mutate(facet_top = "A"), 
               mapping = aes(xintercept =  as.numeric(as.Date(vacc_start_date)-as.Date(params$date0))),
               linetype="longdash", colour = 'grey60') +
    geom_line(aes(t, value, group = scenario)) +
    ggh4x::facet_nested(~ facet_top + scenario, 
                        labeller = as_labeller(c(A    = scen_name, label_facets)) ) +
    scale_x_continuous(breaks = seq(0,
                                    as.numeric(as.Date(params$time1)-as.Date(params$date0)), 
                                    by = 365),
                       labels = c(seq(0,
                                      as.numeric(as.Date(params$time1)-as.Date(params$date0)), 
                                      by = 365)/365),
                       expand = c(0, 0.2), limits = c(0, NA)) +
    scale_y_continuous(label = scales::comma,
                       expand = c(0, 0.2), n.breaks = 6,
                       limits = c(0, define_yaxis_limit)) +
    labs(x="time (year)", y="incidence of new symptomatic cases", caption = cap_text) +
    theme_bw() +
    theme(#axis.title  = element_text(size=14), 
      strip.text  = element_text(size=14),
      panel.spacing.x=unit(0.75, "lines"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      #plot.margin = unit(c(t=0.5, r=0.5, b=0.5, l=2), "lines"),
      axis.text.x = element_text(#angle = 45, hjust = 1,
                                 size=12) )
}


# function for estmating costs and QALYs of different AEFI
scenario_AEFI <- function(epi_data,
                          scen_name       = 0, 
                          QALYs_AEFI      = NULL,
                          costs_AEFI      = NULL,
                          doses_nr        = 1,     # risk of AEFI only for 1 dose; otherwise assuming twice the QALYs lost from death etc...
                          pr_AEFI         = 12345, 
                          base_case       = FALSE,
                          threshold_value = 20000,
                          saveNatOutcomes = TRUE, 
                          includeEconomy  = TRUE){ 
  
  # AEFI base case for AZ
  if(base_case == TRUE){
    
    # QALYs lost on AEFI per vaccinated individual by age
    # mix of local and systemic AEFI from AZ trials for the prime dose (18-55 years: 87%, 56-69: 75%, 70+: 63%)
    # and rare but fatal AEFIs reported recently (1 in 50k in <40, and 1 in 100k in >=40; PHE leaflet based on MHRA data)
    QALY_vaccAEFI_m <<- 
      c(sev_AEFI$QALYs[[1]] * c(rep(0, 3), rep(0.87, 8), rep(0.75, 3), rep(0.63, 2)) +
          sev_AEFI$QALYs[[7]] * c(rep(0, 3), rep(1/50000, 5), rep(1/100000, 8))) %>%
      setNames(c(params$pop[[1]]$group_names))
    
    
    # costs on AEFI per vaccinated individual by age; assuming low and high estimates
    cost_vaccAEFI_m <<- 
      c(sev_AEFI$costs[[1]] * c(rep(0, 3), rep(0.87, 8), rep(0.75, 3), rep(0.63, 2)) +
          sev_AEFI$costs[[7]] * c(rep(0, 3), rep(1/50000, 5), rep(1/100000, 8))) %>%
      setNames(c(params$pop[[1]]$group_names))
  
  } else{
    
    # change assumptions on costs and QALYs of AEFIs for this analysis
    cost_vaccAEFI_m <<- 
      (rep(costs_AEFI*doses_nr, length(params$pop[[1]]$group_names)) * pr_AEFI) %>%
      setNames(c(params$pop[[1]]$group_names))
    
    if(length(QALYs_AEFI)==1) {
      QALY_vaccAEFI_m <<- 
        (rep(QALYs_AEFI*doses_nr, length(params$pop[[1]]$group_names)) * pr_AEFI) %>%
        setNames(c(params$pop[[1]]$group_names))
    } else(
      QALY_vaccAEFI_m <<- 
        (QALYs_AEFI*doses_nr * pr_AEFI) %>%
        setNames(c(params$pop[[1]]$group_names))
    )

    # add scenario name
    if(length(QALYs_AEFI)>1) {QALYs_AEFI = 20}
  }
  
  
  # run econ model
  df <- econ_fun_barplot(data_set      = epi_data, 
                         scenario_name = pr_AEFI,
                         lambda_val    = threshold_value,
                         incl_GDPloss  = includeEconomy,
                         saveOutcomes  = saveNatOutcomes)
  
  # add scenario name and save results
  if(length(QALYs_AEFI)>1) {QALYs_AEFI = 20}
  
  df[[1]] <- df[[1]] %>% mutate(epi_rate = scen_name, 
                                costs_AEFI   = costs_AEFI, QALYs_AEFI = QALYs_AEFI, 
                                lambda_value = threshold_value,
                                costs_AEFI_m = mean(cost_vaccAEFI_m), QALYs_AEFI_m = mean(QALY_vaccAEFI_m))
  df[[2]] <- df[[2]] %>% mutate(epi_rate     = scen_name,  ScenA = pr_AEFI, 
                                costs_AEFI   = costs_AEFI, QALYs_AEFI = QALYs_AEFI,
                                lambda_value = threshold_value, 
                                costs_AEFI_m = mean(cost_vaccAEFI_m), QALYs_AEFI_m = mean(QALY_vaccAEFI_m))
  df[[3]] <- df[[3]] %>% mutate(epi_rate     = scen_name,  ScenA = pr_AEFI, 
                                costs_AEFI   = costs_AEFI, QALYs_AEFI = QALYs_AEFI,
                                lambda_value = threshold_value, 
                                costs_AEFI_m = mean(cost_vaccAEFI_m), QALYs_AEFI_m = mean(QALY_vaccAEFI_m))
  df[[4]] <- df[[4]] %>% mutate(epi_rate     = scen_name,  ScenA = pr_AEFI, 
                                costs_AEFI   = costs_AEFI, QALYs_AEFI = QALYs_AEFI,
                                lambda_value = threshold_value, 
                                costs_AEFI_m = mean(cost_vaccAEFI_m), QALYs_AEFI_m = mean(QALY_vaccAEFI_m))
  
  # return all lists
  return(df)
}

# probability/risk of AEFI
prob_AEFI  = c(0,
               1/1000000,
               1/100000,
               1/10000,
               0.001, 0.01, seq(0.1,0.9, 0.1))

# severity of AEFI
sev_AEFI = list("QALYs" = list(0.0015,                   # ~1 day of flu
                               round(1/365.25,4),        # ~2 days of flu (1 QALD lost) 
                               round(3.5/365.25,4),      # ~1 week of flu
                               #similar to 0.018 QALYs lost from hospital admission due to seasonal influenza or COVID-19
                               round(7/365.25,4),        # ~2 weeks of flu
                               0.05,                     # ~4 weeks of flu
                               1,                        # 1 QALY lost pp
                               unlist(QALY_UK[, "disc"]) ),                      # 20 QALYs lost by premature death at age ~30-39
                "costs" = list(3,                        # (39*0.05+1) ~1 day of flu; 39 GP visit + paracetamol etc(39*0.05+1)
                               3,                        # (39*0.05+1) ~2 days of flu; 39 GP visit + paracetamol etc
                               40,                       # (39+1)  ~1 week of flu (possibly 2 GP visits)
                               (2*40 + 2618),            # weekly GP visits, non-elective long stay bronchopneumonia 
                               (4*40 + 2618),            # weekly GP visits, non-elective long stay bronchopneumonia 
                               (12*40 + 2618 + 10*1467), # 1 QALY lost pp; assume 10 days of ICU for adults
                               (40    + 2618 + 10*1467)) # 20 QALYs lost by premature death at age ~30-39
)

# proxy of 20 QALY loss for those aged 30-40; may account for impact in others too in younger/older ages (eg relatives?)
#mean(unlist(QALY_UK[QALY_UK$group %in% c("30-34","35-39"), "disc"]))
#range(unlist(QALY_UK[, "disc"]))

# explore lower and higher costs per QALY, too (Claxton et al. estimates vs VSL from UK Treasury)
threshold_vals_list = c(15000, 20000, 60000)


# delay of vaccination roll-out in ages <50 year (assumed 2 months if AZ wasn't used anymore)
# general issue of short immunity assumptions (for pessimistic scenario) and the model not having been fitted; 
# thus can only be exploratory and provide some indications (= delay is bad, but by how much)
# also, trying to not link it too closely to UK situation to keep the study relevant..
delay_vals_list = c(0, 56)

########################################
########################################
# deterministic analysis
########################################
########################################


#####
## 1) no-lockdown scenario
#####

fun_lockdown(lockdown_ini  = 6*1e7/55977178*100000,
             lockdown_sub  = 6*1e7/55977178*100000,
             lockdown_sub2 = 6*1e7/55977178*100000)


# run epi model
res_epi_raw999 = epi_mod(VE_B_disease   = rep(0.56, 16), # run vaccination scenarios  # ~10 sec
                         VE_B_infection = rep(0.13, 16), 
                         VE_C_disease   = rep(0.7,  16), 
                         VE_C_infection = rep(0.7,  16), 
                         vacc_delay_time = 0, 
                         run_determ     = TRUE)

res_epi_raw999$epi_rate = 999


# plot epi curves
yaxis_limit_epicurves = 100000

scen_epi_pvacc999 <- epi_fun_vaccplot(data_set           = res_epi_raw999, 
                                      vacc_start_date    = as.Date("2020-12-08"),#+28 for protection after 28 days; not here for plotting though
                                      scen_name          = "no lockdown",
                                      define_yaxis_limit = yaxis_limit_epicurves)


# base case scenario (only run for "no delay" and 20k per QALY)
df_AEFI <- scenario_AEFI(epi_data  = res_epi_raw999, 
                         threshold_value = 20000, 
                         doses_nr  = 1, # risk of AEFI only for 1 dose; otherwise assuming twice the QALYs lost from death etc
                         pr_AEFI   = 54321, 
                         base_case = TRUE, # use mix of observed AEFIs for AZ
                         scen_name = paste0(999, "_", 54321))

# save
scen_econ_AEFI999     <- data.table::rbindlist(df_AEFI[seq(1,length(df_AEFI),4)])
scen_econRes_AEFI999  <- data.table::rbindlist(df_AEFI[seq(2,length(df_AEFI),4)])
scen_econAge_AEFI999  <- data.table::rbindlist(df_AEFI[seq(3,length(df_AEFI),4)])
scen_econAge2_AEFI999 <- data.table::rbindlist(df_AEFI[seq(4,length(df_AEFI),4)])


# threshold analysis on rates and severity of AEFI; cycle through all pre-specified values
system.time({
  df_AEFI <- sapply(delay_vals_list, function(zz){ #1501/60= ~25 min
    
    # run vaccination scenarios  # ~10 sec
    df_base <- epi_mod(VE_B_disease   = rep(0.56, 16), 
                       VE_B_infection = rep(0.13, 16), 
                       VE_C_disease   = rep(0.7,  16), 
                       VE_C_infection = rep(0.7,  16), 
                       vacc_delay_time = zz, 
                       run_determ     = TRUE)
    
    sapply(threshold_vals_list, function(z){ #317/60= ~9.2 min
      sapply(prob_AEFI, function(x){ 
        sapply(1:length(sev_AEFI$QALYs), function(y){ 
          
          scenario_AEFI(epi_data = df_base, 
                        threshold_value = z, 
                        doses_nr = 1, # risk of AEFI only for 1 dose; otherwise assuming twice the QALYs lost from death etc
                        pr_AEFI  = x, 
                        QALYs_AEFI = sev_AEFI$QALYs[[y]], 
                        costs_AEFI = sev_AEFI$costs[[y]], 
                        scen_name  = paste0(999, "_", zz))
        })
      })
    })
  })
})


# save
scen_econ_AEFI999     <- rbind(scen_econ_AEFI999, data.table::rbindlist(df_AEFI[seq(1,length(df_AEFI),4)]),
                               fill=TRUE )
scen_econRes_AEFI999  <- rbind(scen_econRes_AEFI999, data.table::rbindlist(df_AEFI[seq(2,length(df_AEFI),4)]),
                               fill=TRUE )
scen_econAge_AEFI999  <- rbind(scen_econAge_AEFI999, data.table::rbindlist(df_AEFI[seq(3,length(df_AEFI),4)]),
                               fill=TRUE )
scen_econAge2_AEFI999 <- rbind(scen_econAge2_AEFI999, data.table::rbindlist(df_AEFI[seq(4,length(df_AEFI),4)]),
                               fill=TRUE )




#####
## 2) initial lockdown only scenario until 21 June
#####

fun_lockdown(lockdown_ini  = lockdown_trigger,
             lockdown_sub  = 100,
             lockdown_sub2 = 100,
             t_noFurtherLockdowns = as.numeric(as.Date("2021-07-19")-as.Date(params$date0))
             )


# run epi model
res_epi_raw998 = epi_mod(VE_B_disease   = rep(0.56, 16), # run vaccination scenarios  # ~10 sec
                         VE_B_infection = rep(0.13, 16), 
                         VE_C_disease   = rep(0.7,  16), 
                         VE_C_infection = rep(0.7,  16), 
                         vacc_delay_time = 0, 
                         run_determ     = TRUE)

res_epi_raw998$epi_rate = 998


# plot epi curves
scen_epi_pvacc998 <- epi_fun_vaccplot(data_set           = res_epi_raw998, 
                                      vacc_start_date    = as.Date("2020-12-08"),#+28 for protection after 28 days; not here for plotting though
                                      scen_name          = "initial lockdowns only",
                                      define_yaxis_limit = yaxis_limit_epicurves)


# base case scenario (only run for "no delay" and 20k per QALY)
df_AEFI <- scenario_AEFI(epi_data  = res_epi_raw998, 
                         threshold_value = 20000, 
                         doses_nr  = 1, # risk of AEFI only for 1 dose; otherwise assuming twice the QALYs lost from death etc
                         pr_AEFI   = 54321, 
                         base_case = TRUE, # use mix of observed AEFIs for AZ
                         scen_name = paste0(998, "_", 54321))

# save
scen_econ_AEFI998     <- data.table::rbindlist(df_AEFI[seq(1,length(df_AEFI),4)])
scen_econRes_AEFI998  <- data.table::rbindlist(df_AEFI[seq(2,length(df_AEFI),4)])
scen_econAge_AEFI998  <- data.table::rbindlist(df_AEFI[seq(3,length(df_AEFI),4)])
scen_econAge2_AEFI998 <- data.table::rbindlist(df_AEFI[seq(4,length(df_AEFI),4)])



# threshold analysis on rates and severity of AEFI; cycle through all pre-specified values
system.time({
  df_AEFI <- sapply(delay_vals_list, function(zz){ #1501/60= ~25 min
    
    # run vaccination scenarios  # ~10 sec
    df_base <- epi_mod(VE_B_disease   = rep(0.56, 16), 
                       VE_B_infection = rep(0.13, 16), 
                       VE_C_disease   = rep(0.7,  16), 
                       VE_C_infection = rep(0.7,  16), 
                       vacc_delay_time = zz, 
                       run_determ     = TRUE)
    
    sapply(threshold_vals_list, function(z){ #317/60= ~9.2 min
      sapply(prob_AEFI, function(x){ 
        sapply(1:length(sev_AEFI$QALYs), function(y){ 
          
          scenario_AEFI(epi_data = df_base, 
                        threshold_value = z, 
                        doses_nr = 1, # risk of AEFI only for 1 dose; otherwise assuming twice the QALYs lost from death etc
                        pr_AEFI  = x, 
                        QALYs_AEFI = sev_AEFI$QALYs[[y]], 
                        costs_AEFI = sev_AEFI$costs[[y]], 
                        scen_name  = paste0(998, "_", zz))
        })
      })
    })
  })
})


# save
scen_econ_AEFI998     <- rbind(scen_econ_AEFI998, data.table::rbindlist(df_AEFI[seq(1,length(df_AEFI),4)]),
                               fill=TRUE )
scen_econRes_AEFI998  <- rbind(scen_econRes_AEFI998, data.table::rbindlist(df_AEFI[seq(2,length(df_AEFI),4)]),
                               fill=TRUE )
scen_econAge_AEFI998  <- rbind(scen_econAge_AEFI998, data.table::rbindlist(df_AEFI[seq(3,length(df_AEFI),4)]),
                               fill=TRUE )
scen_econAge2_AEFI998 <- rbind(scen_econAge2_AEFI998, data.table::rbindlist(df_AEFI[seq(4,length(df_AEFI),4)]),
                               fill=TRUE )






#####
## 3) intermittent lockdown scenarios
#####

# run for 100 per 100,000; # stricter rates not realistic as no political will for another lockdown
scen_incidence_rate <- 100

# run for all rates
system.time({for(i in 1:length(scen_incidence_rate)){ # 1920/60 = ~32 min
  
  fun_lockdown(lockdown_ini  = lockdown_trigger,
               lockdown_sub  = scen_incidence_rate[i],
               lockdown_sub2 = scen_incidence_rate[i])
  
  # run epi model
  res_epi_raw = epi_mod(VE_B_disease   = rep(0.56, 16), # run vaccination scenarios  # ~10 sec
                       VE_B_infection  = rep(0.13, 16), 
                       VE_C_disease    = rep(0.7,  16), 
                       VE_C_infection  = rep(0.7,  16), 
                       vacc_delay_time = 0, 
                       run_determ      = TRUE)
  
  res_epi_raw$epi_rate = scen_incidence_rate[i]
  
  assign(paste0("res_epi_raw", scen_incidence_rate[i]), res_epi_raw, envir=.GlobalEnv)
  
  
  # plot epi curves
  scen_epi_pvacc <- epi_fun_vaccplot(
    data_set           = res_epi_raw, 
    vacc_start_date    = as.Date("2020-12-08"),#+28 for protection after 28 days; not here for plotting though
    scen_name          = paste0("incidence trigger of physical distancing: ", scen_incidence_rate[i], " per 100,000"),
    define_yaxis_limit = yaxis_limit_epicurves)
  
  assign(paste0("scen_epi_pvacc", scen_incidence_rate[i]), scen_epi_pvacc, envir=.GlobalEnv)
  
  
  # base case scenario (only run for "no delay" and 20k per QALY)
  df_AEFI <- scenario_AEFI(epi_data  = res_epi_raw, 
                           threshold_value = 20000, 
                           doses_nr  = 1, # risk of AEFI only for 1 dose; otherwise assuming twice the QALYs lost from death etc
                           pr_AEFI   = 54321, 
                           base_case = TRUE, # use mix of observed AEFIs for AZ
                           scen_name = paste0(scen_incidence_rate[i], "_", 54321))
  
  # save 
  assign(paste0("scen_econ_AEFI", scen_incidence_rate[i]), 
         data.table::rbindlist(df_AEFI[seq(1,length(df_AEFI),4)]), envir=.GlobalEnv)
  assign(paste0("scen_econRes_AEFI", scen_incidence_rate[i]), 
         data.table::rbindlist(df_AEFI[seq(2,length(df_AEFI),4)]), envir=.GlobalEnv)
  assign(paste0("scen_econAge_AEFI", scen_incidence_rate[i]), 
         data.table::rbindlist(df_AEFI[seq(3,length(df_AEFI),4)]), envir=.GlobalEnv)
  assign(paste0("scen_econAge2_AEFI", scen_incidence_rate[i]), 
         data.table::rbindlist(df_AEFI[seq(4,length(df_AEFI),4)]), envir=.GlobalEnv)
  
  
  # threshold analysis on rates and severity of AEFI; cycle through all pre-specified values
  df_AEFI <- sapply(delay_vals_list, function(zz){ #317/60= ~9.2 min
    
    # run vaccination scenarios  # ~10 sec
    df_base <- epi_mod(VE_B_disease   = rep(0.56, 16), 
                       VE_B_infection = rep(0.13, 16), 
                       VE_C_disease   = rep(0.7,  16), 
                       VE_C_infection = rep(0.7,  16), 
                       vacc_delay_time = zz,
                       run_determ     = TRUE)
    
    sapply(threshold_vals_list, function(z){
    sapply(prob_AEFI, function(x){ 
      sapply(1:length(sev_AEFI$QALYs), function(y){ 
        
        scenario_AEFI(epi_data = df_base, 
                      threshold_value = z,
                      doses_nr = 1,
                      pr_AEFI = x, 
                      QALYs_AEFI = sev_AEFI$QALYs[[y]], 
                      costs_AEFI = sev_AEFI$costs[[y]], 
                      scen_name  = paste0(scen_incidence_rate[i], "_", zz)) })
    })
    })
  })
  
  # save 
  assign(paste0("scen_econ_AEFI", scen_incidence_rate[i]), 
         rbind(get(paste0("scen_econ_AEFI", scen_incidence_rate[i])), 
               data.table::rbindlist(df_AEFI[seq(1,length(df_AEFI),4)]),
               fill=TRUE), envir=.GlobalEnv)
  assign(paste0("scen_econRes_AEFI", scen_incidence_rate[i]), 
         rbind(get(paste0("scen_econRes_AEFI", scen_incidence_rate[i])), 
               data.table::rbindlist(df_AEFI[seq(2,length(df_AEFI),4)]),
               fill=TRUE), envir=.GlobalEnv)
  assign(paste0("scen_econAge_AEFI", scen_incidence_rate[i]), 
         rbind(get(paste0("scen_econAge_AEFI", scen_incidence_rate[i])), 
               data.table::rbindlist(df_AEFI[seq(3,length(df_AEFI),4)]),
               fill=TRUE), envir=.GlobalEnv)
  assign(paste0("scen_econAge2_AEFI", scen_incidence_rate[i]), 
         rbind(get(paste0("scen_econAge2_AEFI", scen_incidence_rate[i])), 
               data.table::rbindlist(df_AEFI[seq(4,length(df_AEFI),4)]),
               fill=TRUE), envir=.GlobalEnv)
}
})


##### SAVE DATA - RECOMMENDED
save.image(file = paste0(econ_path, "AEFI_res_", format(Sys.time(), "%d%b%Y"), ".RData"))



#########################
## plot results
#########################

## combine results
scen_econ_abs_ttl  <- rbindlist(mget(c(ls(pattern = "scen_econ_AEFI\\d+")[grepl("999", ls(pattern = "scen_econ_AEFI\\d+"))], 
                                       ls(pattern = "scen_econ_AEFI\\d+")[grepl("998", ls(pattern = "scen_econ_AEFI\\d+"))], 
                                       ls(pattern = "scen_econ_AEFI\\d+")[!grepl("999|998", ls(pattern = "scen_econ_AEFI\\d+"))])))
scen_econ_res_ttl  <- rbindlist(mget(c(ls(pattern = "scen_econRes_AEFI\\d+")[grepl("999", ls(pattern = "scen_econRes_AEFI\\d+"))], 
                                       ls(pattern = "scen_econRes_AEFI\\d+")[grepl("998", ls(pattern = "scen_econRes_AEFI\\d+"))], 
                                       ls(pattern = "scen_econRes_AEFI\\d+")[!grepl("999|998", ls(pattern = "scen_econRes_AEFI\\d+"))])))
scen_econ_absAge_ttl  <- rbindlist(mget(c(ls(pattern = "scen_econAge_AEFI\\d+")[grepl("999", ls(pattern = "scen_econAge_AEFI\\d+"))], 
                                          ls(pattern = "scen_econAge_AEFI\\d+")[grepl("998", ls(pattern = "scen_econAge_AEFI\\d+"))], 
                                          ls(pattern = "scen_econAge_AEFI\\d+")[!grepl("999|998", ls(pattern = "scen_econAge_AEFI\\d+"))])))
scen_econ_resAge_ttl  <- rbindlist(mget(c(ls(pattern = "scen_econAge2_AEFI\\d+")[grepl("999", ls(pattern = "scen_econAge2_AEFI\\d+"))], 
                                          ls(pattern = "scen_econAge2_AEFI\\d+")[grepl("998", ls(pattern = "scen_econAge2_AEFI\\d+"))], 
                                          ls(pattern = "scen_econAge2_AEFI\\d+")[!grepl("999|998", ls(pattern = "scen_econAge2_AEFI\\d+"))])))


# separate epi scenario from delay time
scen_econ_abs_ttl = scen_econ_abs_ttl %>% tidyr::separate(epi_rate, c("epi_rate", "delay_time")) %>%
  mutate(epi_rate   = as.numeric(epi_rate), 
         delay_time = as.numeric(delay_time))
scen_econ_res_ttl = scen_econ_res_ttl %>% tidyr::separate(epi_rate, c("epi_rate", "delay_time")) %>%
  mutate(epi_rate   = as.numeric(epi_rate), 
         delay_time = as.numeric(delay_time))
scen_econ_absAge_ttl = scen_econ_absAge_ttl %>% tidyr::separate(epi_rate, c("epi_rate", "delay_time")) %>%
  mutate(epi_rate   = as.numeric(epi_rate), 
         delay_time = as.numeric(delay_time))
scen_econ_resAge_ttl = scen_econ_resAge_ttl %>% tidyr::separate(epi_rate, c("epi_rate", "delay_time")) %>%
  mutate(epi_rate   = as.numeric(epi_rate), 
         delay_time = as.numeric(delay_time))

# reorder scenarios to start with no-lockdown
order_scen <- c(unique(scen_econ_abs_ttl$epi_rate)[grepl("999", unique(scen_econ_abs_ttl$epi_rate))],
                unique(scen_econ_abs_ttl$epi_rate)[grepl("998", unique(scen_econ_abs_ttl$epi_rate))],
                unique(scen_econ_abs_ttl$epi_rate)[!grepl("999|998", unique(scen_econ_abs_ttl$epi_rate))] %>% sort())

# labels for figures
label_scen <-  c(Benefit     = "(A): health value (excl. costs)",
                 NB_HC       = "(B): net health value (incl. healthcare costs)",
                 NB_lockdown = "(C): net health value (incl. macro-economics costs)",
                 runA    = "V0: no vaccination", 
                 runB    = "V1: lower VE estimate,\n45-week protection", 
                 runC    = "V2: higher VE estimate,\n3-year protection",
                 `15000` = "£15,000 per QALY", 
                 `20000` = "£20,000 per QALY", 
                 `30000` = "£30,000 per QALY", 
                 `60000` = "£60,000 per QALY",
                 `1e+05` = "£100,000 per QALY", 
                 `1e-06` = "rare but fatal events\n(rate: 1 in 1 million, QALY loss: 4.9-23.0, costs: £17,328)", 
                 `0.001` = "intermediate events\n(rate: 1 in 1000, QALY loss: 1.0, costs: £17,768)", 
                 `0.9`   = "frequent but mild events\n(rate: 9 in 10, QALY loss: 0.0015, costs: £3)",
                 `999`   = "no\nlockdown",
                 `998`   = "initial lockdowns\nonly",
                 `0-19`  = "0-19 y.",
                 `0-14`  = "0-14 y.",
                 `15-29` = "15-29 y.",
                 `20-29` = "20-29 y.",
                 `30-39` = "30-39 y.",
                 `40-59` = "40-59 y.",
                 `40-49` = "40-49 y.",
                 `50-59` = "50-59 y.",
                 `60+`   = "60+ y.",
                 `12345` = "stop",
                 `54321` = "base case",
                 `0`     = "no delay",
                 `7`     = "1 week delay",
                 `14`    = "2 weeks delay",
                 `21`    = "3 weeks delay",
                 `28`    = "4 weeks delay",
                 `56`    = "2 months delay",
                 `84`    = "3 months delay",
                 `10`    = "PD trigger:\n10 per 100,000",
                 `20`    = "PD trigger:\n20 per 100,000",
                 `30`    = "PD trigger:\n30 per 100,000",
                 `40`    = "PD trigger:\n40 per 100,000",
                 `50`    = "PD trigger:\n50 per 100,000",
                 `60`    = "PD trigger:\n60 per 100,000",
                 `100`   = "PD trigger:\n100 per 100,000")

# labels for AEFI scenarios
label_AEFI <-  c(`0.0015` = "~1 day of symptoms (0.0015 QALYs lost)",
                 `0.0027` = "~2 days of symptoms (0.0027 QALYs lost)",
                 `0.0096` = "~1 week of symptoms (0.0096 QALYs lost)",
                 `0.0192` = "~2 weeks of symptoms (0.0192 QALYs lost)",
                 `0.0383` = "~4 weeks of symptoms (0.0383 QALYs lost)",
                 `0.05`   = "~5 weeks of symptoms (0.05 QALYs lost)",
                 `1`      = "serious, non-life-threatening (1 QALY lost)",
                 `20`     = "serious, life-threatening (5-23 QALYs lost)")

# custom theme for plotting
fun_theme = function(){
  
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18, vjust=3.5),
        strip.text.x = element_text(size=20),
        strip.text.y = element_text(size=16),
        legend.title = element_text(size=16),
        legend.text  = element_text(size=16),
        plot.caption = element_text(size=16),
        axis.line    = element_line(colour="grey80"),
        panel.border = element_blank(),
        axis.text.x  = element_text(size=16, angle = 45, vjust = 1, hjust = 1),
        axis.text.y  = element_text(size=18),
        plot.margin  = unit(c(t=0.5, r=0.5, b=1.5, l=2), "lines"))
}


# add health value (gross, ie without costs)
df_AEFI2 <- scen_econ_res_ttl %>%
  mutate(
    parm = ifelse(grepl("QALYs", compartment), "health_loss", 
           ifelse(grepl("costs_lockdown", compartment), "expenses_lockdown", 
           ifelse(grepl("costs", compartment), "expenses_HC",
                  NA))) ) %>%
  dplyr::filter(!grepl("costs", compartment)) %>%
  # get sum of QALYs
  dplyr::group_by(ScenA, scenario, run, parm, epi_rate, lambda_value, delay_time,
                  costs_AEFI, QALYs_AEFI, costs_AEFI_m, QALYs_AEFI_m) %>% 
  summarise_if(is.numeric, sum, na.rm=TRUE) %>% 
  dplyr::group_by(ScenA, scenario, run, parm, epi_rate, delay_time, lambda_value) %>% 
  tidyr::spread(parm, disc) %>%
  mutate(indicator = "gross health benefit") %>%
  mutate(expenses_HC       = 0, # set to 0 for net value estimates to ignore costs here
         expenses_lockdown = 0, # set to 0 for net value estimates to ignore costs here
         Benefit           = -health_loss) %>%  # negative due to estimating QALY losses!
  tidyr::gather(parm, disc, Benefit) %>%
  mutate(disc = ifelse(disc==0, NA, disc) ) %>%
  dplyr::select(-incr_undisc, -incr_disc) %>%
  full_join(scen_econ_abs_ttl)


# add health value by age (gross, ie without costs)
df_AEFI_age <- scen_econ_resAge_ttl %>%
  dplyr::select(-disc_MB) %>%
  dplyr::rename(disc = disc_HB) %>%
  mutate(
    parm = ifelse(grepl("QALYs", compartment), "health_loss", 
           ifelse(grepl("costs_lockdown", compartment), "expenses_lockdown", 
           ifelse(grepl("costs", compartment), "expenses_HC",
                  NA))) ) %>%
  dplyr::filter(!grepl("costs", compartment)) %>%
  # get sum of QALYs
  dplyr::group_by(ScenA, scenario, run, year, group, parm, epi_rate, lambda_value, delay_time,
                  costs_AEFI, QALYs_AEFI, costs_AEFI_m, QALYs_AEFI_m) %>% 
  summarise_if(is.numeric, sum, na.rm=TRUE) %>% 
  dplyr::group_by(ScenA, scenario, run, year, group, parm, epi_rate, delay_time, lambda_value) %>% 
  tidyr::spread(parm, disc) %>%
  mutate(indicator = "gross health benefit") %>%
  mutate(expenses_HC       = 0, # set to 0 for net value estimates to ignore costs here
         expenses_lockdown = 0, # set to 0 for net value estimates to ignore costs here
         Benefit           = -health_loss) %>%     # negative due to QALY losses!
  tidyr::gather(parm, disc, Benefit) %>%
  mutate(disc = ifelse(disc==0, NA, disc) ) %>%
  full_join(scen_econ_absAge_ttl) %>%
  mutate(expenses_lockdown = ifelse(indicator == "health benefit", 0, expenses_lockdown))

# incremental net value of vaccination versus "no-vaccination" (==runA)
# similar to how the CEA would have presented itself at the start
df_AEFI_incr = df_AEFI2 %>% ungroup() %>%
  dplyr::filter(scenario == "runA", !ScenA == 0) %>%
  dplyr::rename(disc_base = disc) %>%
  dplyr::select(ScenA, run, parm, indicator, lambda_value, epi_rate, delay_time, costs_AEFI, QALYs_AEFI, disc_base) %>%
  full_join(df_AEFI2 %>% dplyr::filter(!ScenA == 0) %>% ungroup, 
            by=c("ScenA", "run", "parm", "indicator", "delay_time", "lambda_value", "epi_rate", "costs_AEFI", "QALYs_AEFI")) %>%
  # benefit vs no-vaccine
  dplyr::mutate(disc = ifelse(disc<0, -(disc_base - disc), disc_base - disc)) %>%
  dplyr::select(-disc_base) %>%
  dplyr::filter(scenario != "runA") %>%
  mutate(epi_rate    = factor(epi_rate, levels = order_scen),
         QALYs_AEFI2 = factor(label_AEFI[paste0(QALYs_AEFI)], levels = label_AEFI)) %>%
  # indicator of likely AEFIs for COVID vaccines (as observed in phase 3 trials and reported in March 2021)
  mutate(disc2 = cut(disc, breaks = c(-Inf, -1e+9, -1e+8, -1e+7, -1e+6, -1e+5,
                                      0, 
                                      1e+5, 1e+6, 1e+7, 1e+8, 1e+9, 1e+12, Inf), 
                     labels = as.factor(c("<-100 mln.", "<-100 mln.", "<-10 mln.", "<-1 mln.", "<-100,000", 
                                          "<0 & >-100,000", ">0 & <100,000", 
                                          ">100,000", ">1 mln.", ">10 mln.", ">100 mln.", ">1 bln.", ">1 trln.")), 
                     include.lowest=TRUE, right=FALSE)
  ) %>%
  mutate(observed = ifelse((QALYs_AEFI %in% c(0.0015, 0.0027) & ScenA %in% paste0(seq(0.5, 0.9, 0.1))) |
                             (QALYs_AEFI == 20 & ScenA %in% c(1e-06, 1e-05)), TRUE, FALSE))


# prepare benefit-risk ratios
# add indicator of costs and QALYs, AEFI harm and vaccine benefit
scen_econ_res_ttl = scen_econ_res_ttl %>% 
  mutate( parm = ifelse(grepl("QALYs", compartment), "health_loss", 
                 ifelse(grepl("costs_lockdown", compartment), "expenses_lockdown", 
                 ifelse(grepl("costs", compartment), "expenses_HC",
                        NA))),
          parm3 = ifelse(grepl("costs_vaccAEFI", compartment), "AEFI costs", 
                  ifelse(grepl("QALYs_vaccAEFI", compartment), "AEFI QALYs", 
                  ifelse(!grepl("AEFI", compartment) & grepl("QALYs", compartment), "health_loss", 
                  ifelse(!grepl("AEFI", compartment) & grepl("costs_lockdown", compartment), "expenses_lockdown",
                  ifelse(!grepl("AEFI", compartment) & grepl("costs", compartment), "expenses_HC",
                         NA)))))
  ) 


# get health value, net value (healthcare costs), and net value (macro-economic costs) by category for AEFI
# health value (no costs)
df_AEFI2b <- scen_econ_res_ttl %>%
  dplyr::select(-incr_undisc, -incr_disc, -undisc)  %>%
  dplyr::group_by(ScenA, scenario, run, parm, parm3, epi_rate, delay_time, lambda_value) %>% 
  tidyr::spread(parm, disc) %>%
  mutate(indicator = "gross health benefit") %>%
  mutate(expenses_HC = 0,
         expenses_lockdown = 0,
         # keep values positive to trade-off harm against benefit (against no vaccine; next step below)
         Benefit       = health_loss) %>% 
  tidyr::gather(parm, disc, Benefit) %>%
  mutate(disc = ifelse(disc==0, NA, disc) ) %>%
  mutate( parmAEFI = ifelse(grepl("AEFI", compartment), "AEFI", "non_AEFI") )

# sum by binary AEFI indicator (leave separate to visualise costs/QALYs separtely too)
df_AEFI2bb <- df_AEFI2b %>%
  # get sum of QALYs
  dplyr::group_by(ScenA, scenario, run, indicator, parm, parmAEFI, #parm3, 
                  epi_rate, delay_time, lambda_value, costs_AEFI, QALYs_AEFI, costs_AEFI_m, QALYs_AEFI_m) %>% 
  summarise_if(is.numeric, sum, na.rm=TRUE)

# get net monetary/health value by AEFI
df_AEFI2c <- scen_econ_res_ttl %>% ungroup() %>%
  dplyr::select(-incr_undisc, -incr_disc, -undisc)  %>%
  mutate(disc_MB = ifelse(grepl("QALY", compartment),  disc*lambda_value, disc),
         disc_HB = ifelse(grepl("costs", compartment), disc/lambda_value, disc)
  ) %>%
  dplyr::select(-disc, -QALYs_AEFI_m, -costs_AEFI_m)


#aggregated NB
df_AEFI2d <- df_AEFI2c %>%
  # net monetary value
  dplyr::select(-disc_HB) %>%
  dplyr::rename(disc = disc_MB) %>%
  mutate(indicator = "monetary benefit") %>%
  # net health value
  bind_rows(
    df_AEFI2c %>%
      dplyr::select(-disc_MB) %>%
      dplyr::rename(disc = disc_HB) %>%
      mutate(indicator = "health benefit")
  ) %>%
  mutate( parmAEFI = ifelse(grepl("AEFI", compartment), "AEFI", "non_AEFI") ) %>%
  # get absolute net values (losses and benefits), similar to health value, for trading off harms vs benefits
  dplyr::mutate(disc = ifelse(is.na(disc), 0, abs(disc)) )


# results by binary AEFI indicator
df_AEFI2dd <- df_AEFI2d %>%
  dplyr::group_by(ScenA, scenario, epi_rate, delay_time, run, indicator, parm, parmAEFI, #parm3, 
                  costs_AEFI, QALYs_AEFI, lambda_value) %>% 
  summarise_if(is.numeric, sum, na.rm=TRUE) %>% 
  tidyr::spread(parm, disc) %>%
  mutate(expenses_lockdown = ifelse(is.na(expenses_lockdown), 0, expenses_lockdown) ) %>%
  mutate(NB_HC       = -health_loss - expenses_HC,     # negative health_loss due to QALY losses!
         NB_lockdown = -health_loss - expenses_lockdown) %>%
  tidyr::gather(parm, disc, NB_HC:NB_lockdown) %>%
  # get absolute net values (losses and benefits), similar to health value, for trading off harms vs benefits
  dplyr::mutate( disc = abs(disc))


# combine health value with net values
df_AEFI2e = df_AEFI2bb %>%
  full_join(df_AEFI2dd)


# get added value against no vaccination
df_AEFI_incr4 = df_AEFI2e %>% ungroup() %>%
  # value of vaccine vs no-vaccine
  dplyr::filter(scenario == "runA") %>% #, ScenA == 0) %>%
  dplyr::rename(disc_base = disc) %>%
  dplyr::select(run, parm, ScenA, parmAEFI, indicator, lambda_value, delay_time, epi_rate, 
                costs_AEFI, QALYs_AEFI, disc_base) %>%
  full_join(df_AEFI2e %>% 
              dplyr::filter(!scenario == "runA") %>% 
              ungroup, 
            by=c("run", "parm", "ScenA", "parmAEFI","indicator", "lambda_value", "delay_time", "epi_rate", 
                 "costs_AEFI", "QALYs_AEFI")) %>%
  # get positive values of harm from AEFI (vs no AEFI), vs positive values of vaccine (vs no vaccine)
  dplyr::mutate( disc = ifelse(parmAEFI=="AEFI", disc - disc_base, disc_base - disc)) %>%
  dplyr::select(-disc_base)


# get ratios
df_AEFI_incr5 = df_AEFI_incr4 %>% ungroup() %>% 
  dplyr::filter(!scenario == "runA", !ScenA == 0) %>% 
  dplyr::select(-health_loss, -expenses_HC, -expenses_lockdown) %>%
  tidyr::spread(parmAEFI, disc) %>% 
  dplyr::mutate(disc3 = ifelse(non_AEFI/AEFI=="Inf", non_AEFI, # take non_AEFI vaccine benefit as only AEFI can be 0
                               ifelse(non_AEFI/AEFI<1, -AEFI/non_AEFI,
                                      non_AEFI/AEFI)) ) %>%
  mutate(epi_rate    = factor(epi_rate, levels = order_scen),
         QALYs_AEFI2 = factor(label_AEFI[paste0(QALYs_AEFI)], levels = label_AEFI)) %>%
  mutate(rel = ifelse(disc3>6, 6, round(disc3,0)),
         rel2 = cut(disc3, breaks = c(-Inf, seq(-6,-1,1), seq(0,7,1), Inf), 
                    labels = as.factor(c(paste0("<", seq(-6,0,1), "x"),
                                         "<1x",
                                         paste0(">", seq(1,6,1), "x"), ">6x")), 
                    include.lowest=TRUE, right=FALSE)
  ) %>%
  mutate(observed = ifelse((QALYs_AEFI %in% c(0.0015, 0.0027) & ScenA %in% paste0(seq(0.5, 0.9, 0.1))) |
                             (QALYs_AEFI == 20 & ScenA %in% c(1e-06, 1e-05)), TRUE, FALSE))



########
# Fig. 1): bar plot of AEFI base case scenario (benefit versus harm)
########

p_AEFI = df_AEFI_incr5 %>%
  dplyr::filter(!ScenA == 0, 
                delay_time==54321,
                epi_rate %in% c(999, 998, 100),
                indicator %in% c("health benefit", "gross health benefit"), 
                lambda_value == 20000) %>%
  mutate(epi_rate    = factor(epi_rate, levels = order_scen)) %>%
  mutate(ScenA2 = ifelse(ScenA == 1e-06, "rare but fatal\n(1e-06, 20.0)",
                  ifelse(ScenA == 0.001, "intermediate\n(0.001, 1.0)",
                  ifelse(ScenA == 0.9,   "frequent but mild\n(0.9, 0.0015)",
                         ScenA))),
         ScenA2 = factor(ScenA2, levels = c("rare but fatal\n(1 in 1,000,000;\n20.0 QALYs lost)", 
                                            "intermediate\n(0.001, 1.0)",
                                            "frequent but mild\n(9 in 10;\n0.0015 QALYs lost)"))) %>%
  tidyr::gather(parm3, disc, AEFI:non_AEFI) %>% 
  mutate( parm2 = ifelse(grepl("^AEFI", parm3), "harm from\nAEFI", "benefit from\nvaccine") )

# reorder
p_AEFI = p_AEFI %>%
  mutate(parm2 = factor(parm2, levels = c(unique(p_AEFI$parm2)[2], 
                                          unique(p_AEFI$parm2)[1])),
         rel3 = ifelse(disc3>=100, paste0(">", signif(floor(disc3/10)*10,2), "x"), 
                       paste0(">", signif(floor(disc3),2), "x"))
  )


# plot
p1 = p_AEFI %>% 
  mutate(disc = disc/1000000) %>%
  ggplot(aes(fill=parm3, y=(disc), x=factor(parm2))) + #
  geom_bar(position="stack", stat="identity", show.legend = T) +
  scale_fill_manual(name = "", 
                    breaks = c("non_AEFI", "AEFI"),
                    values = c("AEFI" = "darkred",  "AEFI" = "red",      
                               "non_AEFI" = "#0080FF", "non_AEFI" = "lightblue"), 
                    labels = c("AEFI" = "harm", "non_AEFI" = "benefit") ) +  
  ggh4x::facet_nested(epi_rate ~ parm+scenario, #scales = "free_y",
                      labeller = as_labeller(c("rare but fatal\n(1e-06, 20.0)"    = "rare but fatal AEFI\n(1 in 1,000,000; 20.0 QALYs lost)",
                                               "intermediate\n(0.001, 1.0)"       = "intermediate AEFI\n(1 in 1,000; 1.0 QALYs lost)",
                                               "frequent but mild\n(0.9, 0.0015)" = "frequent but mild AEFI\n(9 in 10; 0.0015 QALYs lost)",
                                               label_scen)) ) +
  theme_bw() +
  scale_y_continuous(n.breaks = 10, expand = c(0, 0)) + #, labels = scales::comma
  expand_limits(y = c(0, NA)) +
  geom_text(aes(x=x_pos, y=disc, label = rel3, vjust = v_just), 
            data = p_AEFI %>% 
              group_by(run, parm, indicator, lambda_value, epi_rate, delay_time, scenario) %>% 
              mutate(disc   = disc/1000000,
                     x_pos  = ifelse(parm=="NB_lockdown" & epi_rate==100 & scenario=="runC", 1.75, 1.5),
                     v_just = ifelse(parm=="NB_lockdown" & epi_rate==100 & scenario=="runC", 1.25, -0.5)) %>%
              dplyr::filter(!parm3=="AEFI"), #vjust = -0.5, 
            angle  = 0,  size   = 8, lineheight = 0.8, colour = "black") + 
  geom_hline(yintercept = 0) +
  labs(x = "benefit of vaccine vs. harm from adverse events following immunisation (AEFI)", 
       y = "value in terms of health\n(expressed in million QALYs)") +
  fun_theme() + 
  theme(panel.spacing.y = unit(1.5, "lines"))

# add footnote
title <- ggdraw() + 
  draw_label(paste0("\nBenefit-risk ratios (`>?x`) indicate the benefit (and net benefit; dQALYs - dCosts/£20,000) of the vaccine over the harm (and net harm) of the AEFIs, discounted at 3.5% over 10 years.
                     V1: vaccination with 62% VE against symptomatic disease, 4% against asymptomatic infection (AstraZeneca-Oxford Uni. vaccine, low estimate); vaccine-induced protection of 45-weeks
                     V2: vaccination with 91% VE against symptomatic disease, 56% against asymptomatic infection (AstraZeneca-Oxford Uni. vaccine, high estimate); vaccine-induced protection of 3-years"),
             size=15, x = 1, y=0.78, hjust = 1) 


cowplot::plot_grid(p1, title, 
                   ncol = 1,
                   rel_heights = c(1, 0.1) )

# save
ggsave(paste0("Fig1-CoVAEFI_determ_", format(Sys.time(), "%d%b%Y"), "_v4.png"), path = econ_path, width = 57, height = 28, units = "cm")



########
# Fig. 3): heat map of AEFI threshold scenario analysis
########

# plot
p1 = df_AEFI_incr5 %>%
  dplyr::filter(!ScenA == 0, 
                delay_time == 0,
                indicator %in% c("health benefit" , 
                                 "gross health benefit"), 
                lambda_value == 20000) %>%
  ggplot(aes(x = factor(ScenA), y= QALYs_AEFI2, fill = rel2)) + 
  geom_tile(aes(color = disc3>0), size=1.5) + #rel>=6), size=1.5) +
  scale_fill_manual(name = "benefit of\nvaccine\nvs harm\nfrom AEFI",
                    values=c("#DC143C",  colorRampPalette(c("red","orange"))(5),
                             colorRampPalette(c("lightblue","blue"))(6) ),
                    na.value="#EEEEEE") +
  scale_colour_manual(name = "", 
                      values = c("FALSE" = "darkred",  "TRUE" = "#0080FF"),
                      labels = c("FALSE" = "negative", "TRUE" = "positive")) +
  ggh4x::facet_nested(epi_rate~parm + scenario, scales = "free",
                      labeller = as_labeller(c(label_scen)) ) +
  guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2, override.aes = list(fill = NA))) +
  theme_bw() +
  geom_point(aes(size=ifelse(observed, "dot", "no_dot")), shape = 4, colour= "white", show.legend = FALSE) +
  scale_size_manual(values=c(dot=4, no_dot=NA), guide="none") +
  labs(x = "risk of adverse event", y = "severity of adverse event") +
  fun_theme() +
  theme( axis.text.x  = element_text(size=14, angle = 45, vjust = 1, hjust = 1),
         legend.title = element_text(size = 18),
         legend.text  = element_text(size = 18) )



# add footnote
title <- ggdraw() + 
  draw_label(paste0("Benefit-risk ratios are positive when the benefit (and net benefit; dQALY_loss - dCosts/£20,000) of the vaccine exceeds the (net) harm of AEFIs, and negative when the (net) harm of AEFIs exceeds the (net) benefit of the vaccine.
                     V1: vaccination with 62% VE against symptomatic disease, 4% against asymptomatic infection (AstraZeneca-Oxford Uni. vaccine, low estimate); assumed vaccine-induced protection of 45-weeks
                     V2: vaccination with 91% VE against symptomatic disease, 56% against asymptomatic infection (AstraZeneca-Oxford Uni. vaccine, high estimate); assumed vaccine-induced protection of 3-years
                     X indicate the observed rates and severities of AEFIs of COVID-19 vaccines."),
             size=15, x = 1, y=0.78, hjust = 1) 

cowplot::plot_grid(p1, title, 
                   ncol = 1,
                   rel_heights = c(1, 0.05) )


# save
ggsave(paste0("Fig3-CoVAEFI_determ_", format(Sys.time(), "%d%b%Y"), "_v2.png"), path = econ_path, width = 70/3*3.0, height = 32.5, units = "cm")


########
# SFig. 1): heat map of AEFI threshold scenario analysis (different costs per QALY)
########

# plot
p1 = df_AEFI_incr5 %>%
  dplyr::filter(!ScenA == 0,
                delay_time == 0,
                indicator %in% c("health benefit" , 
                                 "gross health benefit"), 
                epi_rate == 100) %>%
  ggplot(aes(x = factor(ScenA), y= QALYs_AEFI2, fill = rel2)) + 
  #geom_tile(aes(color = disc>=6), size=1.5) + #rel>=6), size=1.5) +
  geom_tile(aes(color = disc3>0), size=1.5) + #rel>=6), size=1.5) +
  scale_fill_manual(#name = "benefit/risk-ratio:\nbenefit of vaccine\n(without AEFI)\nvs harm from AEFI",
    name = "benefit of\nvaccine\nvs harm\nfrom AEFI",
    values=c("#DC143C",  colorRampPalette(c("red","orange"))(5),
             colorRampPalette(c("lightblue","blue"))(6) ),
    na.value="#EEEEEE") +
  scale_colour_manual(name = "", 
                      values = c("FALSE" = "darkred",  "TRUE" = "#0080FF"),
                      labels = c("FALSE" = "negative", "TRUE" = "positive")) +
  ggh4x::facet_nested(epi_rate+lambda_value~parm + scenario, scales = "free",
                      labeller = as_labeller(c(label_scen)) ) +
  guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2, override.aes = list(fill = NA))) +
  theme_bw() +
  geom_point(aes(size=ifelse(observed, "dot", "no_dot")), shape = 4, colour= "white", show.legend = FALSE) +
  scale_size_manual(values=c(dot=4, no_dot=NA), guide="none") +
  labs(x = "risk of adverse event", y = "severity of adverse event") +
  fun_theme()


# add footnote
title <- ggdraw() + 
  draw_label(paste0("Benefit-risk ratios are positive when the benefit (and net benefit; dQALY_loss - dCosts/£20,000) of the vaccine exceeds the (net) harm of AEFIs, and negative when the (net) harm of AEFIs exceeds the (net) benefit of the vaccine.
                     V1: vaccination with 62% VE against symptomatic disease, 4% against asymptomatic infection (AstraZeneca-Oxford Uni. vaccine, low estimate); assumed vaccine-induced protection of 45-weeks
                     V2: vaccination with 91% VE against symptomatic disease, 56% against asymptomatic infection (AstraZeneca-Oxford Uni. vaccine, high estimate); assumed vaccine-induced protection of 3-years
                     X indicate the observed rates and severities of AEFIs of COVID-19 vaccines."),
             size=15, x = 1, y=0.78, hjust = 1) 

cowplot::plot_grid(p1, title, 
                   ncol = 1,
                   rel_heights = c(1, 0.05) )


# save
ggsave(paste0("FigS1-CoVAEFI_determ_", format(Sys.time(), "%d%b%Y"), "_v2.png"), path = econ_path, width = 70/3*3.0, height = 30, units = "cm")



########
# Fig. 2): bar plot of AEFI base case scenario by age
########

# incremental
p_incr_age = df_AEFI_age %>% ungroup() %>%
  dplyr::filter(scenario == "runA", 
                lambda_value == 20000) %>%
  dplyr::rename(disc_base = disc) %>%
  dplyr::select(group, year, ScenA, run, parm, indicator, lambda_value, delay_time, epi_rate, costs_AEFI, QALYs_AEFI, disc_base) %>%
  full_join(df_AEFI_age %>% ungroup() %>%
              dplyr::filter(!scenario == "runA", 
                            lambda_value == 20000) ) %>%
  dplyr::mutate( disc = disc - disc_base) %>%
  dplyr::select(-disc_base) %>% 
  dplyr::filter(lambda_value == 20000) %>% 
  mutate(epi_rate = factor(epi_rate, levels = order_scen),
         year     = 2019 + year) %>%
  mutate(group2 = ifelse(group %in% c("0-4", "5-9", "10-14"),"0-14",
                         ifelse(group %in% c("15-19", "20-24", "25-29"),"15-29",
                                ifelse(group %in% c("30-34", "35-39"),"30-39",
                                       ifelse(group %in% c("40-44", "45-49"),"40-49",
                                              ifelse(group %in% c("50-54", "55-59"),"50-59",
                                                     ifelse(group %in% c("60-64", "65-69", "70-74", "75+"),"60+",
                                                            group))))))
  ) %>%
  group_by(group2, year, ScenA, scenario, run, indicator, parm, epi_rate, delay_time, QALYs_AEFI, lambda_value) %>%
  summarise_if(is.numeric, sum, na.rm=TRUE)




# plot AZ risk
p_incr_age2 = p_incr_age %>%
  dplyr::filter(delay_time == 54321) %>% 
  dplyr::filter(!year == 2020) %>%
  # get background colour when positive or negative
  group_by(epi_rate, delay_time, indicator, parm, ScenA, group2) %>%
  mutate(pos = ifelse(disc<0, 1, 0)) %>%
  group_by(epi_rate, delay_time, indicator, parm, ScenA, group2) %>%
  mutate(pos2 = ifelse(mean(pos)==1, "negative", 
                ifelse(mean(pos)==0, "positive", 
                       "mixed")),
         pos2 = factor(pos2, levels = c("positive", "negative", "mixed")))



p1a = p_incr_age2 %>%
  dplyr::filter(!group2 == "60+",
                indicator %in% c("gross health benefit")) %>% #, "health benefit"
  ggplot(aes(x = factor(year), y=disc, fill = scenario)) +
  geom_rect(aes(fill = pos2), xmin = -Inf, xmax = Inf, ymin = -Inf,ymax = Inf, alpha = 0.01) +
  ggh4x::facet_nested(epi_rate~ #parm+#ScenA+
                        group2, scales = "fixed", labeller = as_labeller(c(label_scen)) ) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  theme_bw() +
  scale_x_discrete(labels = c(rbind(unique(p_incr_age2$year)[(2*1:ceiling(length(unique(p_incr_age2$year))/2))-1], "")) ) +
  scale_y_continuous(n.breaks = 3, labels = function(x) format(x, big.mark = ",", scientific = F)) +
  geom_hline(yintercept = 0) +
  labs(x = "time (calendar years)", 
       y = "health value in terms of QALYs\n(excl. costs)") +
      #y = "value (in terms of QALYs)") +
  fun_theme() +
  theme(legend.position = "bottom",
        strip.text.y = element_blank(),
        strip.background.y = element_blank()) +
  scale_fill_manual(name = "", 
                    breaks = c("positive",  "negative", "mixed", "runB",  "runC"),
                    values = c("positive" = "lightblue",  "negative" = "#DC143C", "mixed" = "orange",
                               "runB" = "#DC143C",  "runC" = "steelblue"), 
                    labels = c("positive" = "positive value (throughout)",  
                               "negative" = "negative value (throughout)", 
                               "mixed"    = "mixed value",
                               "runB"     = "V1: lower VE estimate,\n45-week protection", 
                               "runC"     = "V2: higher VE estimate,\n3-year protection")) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE#, 
                           #override.aes = list(alpha = c(0.2,0.2,0.2,1,1))
                           ))

p1b = p_incr_age2 %>%
  dplyr::filter(group2 == "60+",
                indicator %in% c("gross health benefit")) %>% #, "health benefit"
  ggplot(aes(x = factor(year), y=disc, fill = scenario)) +
  geom_rect(aes(fill = pos2), xmin = -Inf, xmax = Inf, ymin = -Inf,ymax = Inf, alpha = 0.01) +
  ggh4x::facet_nested(epi_rate~ #parm+#ScenA+
                        group2, scales = "fixed", labeller = as_labeller(c(label_scen)) ) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  theme_bw() +
  scale_x_discrete(labels = c(rbind(unique(p_incr_age2$year)[(2*1:ceiling(length(unique(p_incr_age2$year))/2))-1], "")) ) +
  scale_y_continuous(n.breaks = 3, labels = function(x) format(x, big.mark = ",", scientific = F)) +
  geom_hline(yintercept = 0) +
  labs(x = "", y = "") +
  fun_theme() +
  theme(legend.position = "bottom",
        plot.margin  = unit(c(t=0.5, r=0.5, b=1.5, l=0.5), "lines")) +
  scale_fill_manual(name = "", 
                    breaks = c("positive",  "negative", "mixed", "runB",  "runC"),
                    values = c("positive" = "lightblue",  "negative" = "#DC143C", "mixed" = "orange",
                               "runB" = "#DC143C",  "runC" = "steelblue"), 
                    labels = c("positive" = "positive value (throughout)",  
                               "negative" = "negative value (throughout)", 
                               "mixed"    = "mixed value",
                               "runB"     = "V1: lower VE estimate,\n45-week protection", 
                               "runC"     = "V2: higher VE estimate,\n3-year protection")) +
  guides(fill=FALSE) #guide_legend(nrow=2, byrow=TRUE, override.aes = list(alpha = c(0.2,0.2,0.2,1,1))))


p2a = p_incr_age2 %>%
  dplyr::filter(!group2 == "60+",
                indicator %in% c("health benefit"),
                parm == "NB_HC") %>%
  ggplot(aes(x = factor(year), y=disc, fill = scenario)) +
  geom_rect(aes(fill = pos2), xmin = -Inf, xmax = Inf, ymin = -Inf,ymax = Inf, alpha = 0.01) +
  ggh4x::facet_nested(epi_rate~ group2, scales = "fixed", labeller = as_labeller(c(label_scen)) ) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  theme_bw() +
  scale_x_discrete(labels = c(rbind(unique(p_incr_age2$year)[(2*1:ceiling(length(unique(p_incr_age2$year))/2))-1], "")) ) +
  scale_y_continuous(n.breaks = 3, labels = function(x) format(x, big.mark = ",", scientific = F)) +
  geom_hline(yintercept = 0) +
  labs(x = "time (calendar years)", 
       y = "net health value in terms of QALYs\n(incl. healthcare costs)") +
  fun_theme() +
  theme(legend.position = "bottom",
        strip.text.y = element_blank(),
        strip.background.y = element_blank()) +
  scale_fill_manual(name = "", 
                    breaks = c("positive",  "negative", "mixed", "runB",  "runC"),
                    values = c("positive" = "lightblue",  "negative" = "#DC143C", "mixed" = "orange",
                               "runB" = "#DC143C",  "runC" = "steelblue"), 
                    labels = c("positive" = "positive value (throughout)",  
                               "negative" = "negative value (throughout)", 
                               "mixed"    = "mixed value",
                               "runB"     = "V1: lower VE estimate,\n45-week protection", 
                               "runC"     = "V2: higher VE estimate,\n3-year protection")) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE, override.aes = list(alpha = c(0.2,0.2,0.2,1,1))))

p2b = p_incr_age2 %>%
  dplyr::filter(group2 == "60+",
                indicator %in% c("health benefit"),
                parm == "NB_HC") %>%
  ggplot(aes(x = factor(year), y=disc, fill = scenario)) +
  geom_rect(aes(fill = pos2), xmin = -Inf, xmax = Inf, ymin = -Inf,ymax = Inf, alpha = 0.01) +
  ggh4x::facet_nested(epi_rate~ group2, scales = "fixed", labeller = as_labeller(c(label_scen)) ) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  theme_bw() +
  scale_x_discrete(labels = c(rbind(unique(p_incr_age2$year)[(2*1:ceiling(length(unique(p_incr_age2$year))/2))-1], "")) ) +
  scale_y_continuous(n.breaks = 3, labels = function(x) format(x, big.mark = ",", scientific = F)) +
  geom_hline(yintercept = 0) +
  labs(x = "", y = "") +
  fun_theme() +
  theme(legend.position = "bottom",
        plot.margin  = unit(c(t=0.5, r=0.5, b=1.5, l=0.5), "lines")) +
  scale_fill_manual(name = "", 
                    breaks = c("positive",  "negative", "mixed", "runB",  "runC"),
                    values = c("positive" = "lightblue",  "negative" = "#DC143C", "mixed" = "orange",
                               "runB" = "#DC143C",  "runC" = "steelblue"), 
                    labels = c("positive" = "positive value (throughout)",  
                               "negative" = "negative value (throughout)", 
                               "mixed"    = "mixed value",
                               "runB"     = "V1: lower VE estimate,\n45-week protection", 
                               "runC"     = "V2: higher VE estimate,\n3-year protection")) +
  guides(fill=FALSE) #guide_legend(nrow=2, byrow=TRUE, override.aes = list(alpha = c(0.2,0.2,0.2,1,1))))



p3a = p_incr_age2 %>%
  dplyr::filter(!group2 == "60+",
                indicator %in% c("health benefit"),
                parm == "NB_lockdown") %>%
  ggplot(aes(x = factor(year), y=disc, fill = scenario)) +
  geom_rect(aes(fill = pos2), xmin = -Inf, xmax = Inf, ymin = -Inf,ymax = Inf, alpha = 0.01) +
  ggh4x::facet_nested(epi_rate~ group2, scales = "fixed", labeller = as_labeller(c(label_scen)) ) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  theme_bw() +
  scale_x_discrete(labels = c(rbind(unique(p_incr_age2$year)[(2*1:ceiling(length(unique(p_incr_age2$year))/2))-1], "")) ) +
  scale_y_continuous(n.breaks = 3, labels = function(x) format(x, big.mark = ",", scientific = F)) +
  geom_hline(yintercept = 0) +
  labs(x = "time (calendar years)", y = "net health value in terms of QALYs\n(incl. macro-economic costs)") +
  fun_theme() +
  theme(legend.position = "bottom",
        strip.text.y = element_blank(),
        strip.background.y = element_blank()) +
  scale_fill_manual(drop = FALSE, name = "", 
                    breaks = c("positive",  "negative", "mixed", "runB",  "runC"),
                    values = c("positive" = "lightblue",  "negative" = "#DC143C", "mixed" = "orange",
                               "runB" = "#DC143C",  "runC" = "steelblue"), 
                    labels = c("positive" = "positive value (throughout)",  
                               "negative" = "negative value (throughout)", 
                               "mixed"    = "mixed value",
                               "runB"     = "V1: lower VE estimate,\n45-week protection", 
                               "runC"     = "V2: higher VE estimate,\n3-year protection")) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE, override.aes = list(alpha = c(0.2,0.2,0.2,1,1))))

p3b = p_incr_age2 %>%
  dplyr::filter(group2 == "60+",
                indicator %in% c("health benefit"), #c("monetary benefit"),
                parm == "NB_lockdown") %>%
  ggplot(aes(x = factor(year), y=disc, fill = scenario)) +
  geom_rect(aes(fill = pos2), xmin = -Inf, xmax = Inf, ymin = -Inf,ymax = Inf, alpha = 0.01) +
  ggh4x::facet_nested(epi_rate~ group2, scales = "fixed", labeller = as_labeller(c(label_scen)) ) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  theme_bw() +
  scale_x_discrete(labels = c(rbind(unique(p_incr_age2$year)[(2*1:ceiling(length(unique(p_incr_age2$year))/2))-1], "")) ) +
  scale_y_continuous(n.breaks = 3, labels = function(x) format(x, big.mark = ",", scientific = F)) +
  geom_hline(yintercept = 0) +
  labs(x = "", y = "") +
  fun_theme() +
  theme(legend.position = "bottom",
        plot.margin  = unit(c(t=0.5, r=0.5, b=1.5, l=0.5), "lines")) +
  scale_fill_manual(drop = FALSE, name = "", 
                    breaks = c("positive",  "negative", "mixed", "runB",  "runC"),
                    values = c("positive" = "lightblue",  "negative" = "#DC143C", "mixed" = "orange",
                               "runB" = "#DC143C",  "runC" = "steelblue"), 
                    labels = c("positive" = "positive value (throughout)",  
                               "negative" = "negative value (throughout)", 
                               "mixed"    = "mixed value",
                               "runB"     = "V1: lower VE estimate,\n45-week protection", 
                               "runC"     = "V2: higher VE estimate,\n3-year protection")) +
  guides(fill=FALSE) #guide_legend(nrow=2, byrow=TRUE, override.aes = list(alpha = c(0.2,0.2,0.2,1,1))))




# plot vertically all 3
p1 = cowplot::plot_grid(p1a + 
                          guides(fill=FALSE) +
                          theme(axis.title.x = element_blank(),
                                axis.text.x  = element_blank(),
                                plot.caption = element_blank()), 
                        p1b + 
                          guides(fill=FALSE) +
                          theme(axis.title.x = element_blank(),
                                axis.text.x  = element_blank(),
                                plot.caption = element_blank()), 
                        nrow = 1, 
                        rel_widths = c(1, 0.3) )

p2 = cowplot::plot_grid(p2a + 
                          guides(fill=FALSE) +
                          theme(axis.title.x = element_blank(),
                                axis.text.x  = element_blank(),
                                plot.caption = element_blank()), 
                        p2b + 
                          guides(fill=FALSE) +
                          theme(axis.title.x = element_blank(),
                                axis.text.x  = element_blank(),
                                plot.caption = element_blank()), 
                        nrow = 1, 
                        rel_widths = c(1, 0.3) )

p3 = cowplot::plot_grid(p3a, p3b, 
                        nrow = 1, 
                        align = 'h', axis = 'tb', 
                        rel_widths = c(1, 0.3) )


# add footnote
title <- ggdraw() + 
  draw_label(paste0("Health values (and net health values; dQALYs - dCosts/£20,000), discounted at 3.5% over 10 years, are negative in case of health losses (and excess costs) compared to no-vaccination (A)
                     V1: vaccination with 62% VE against symptomatic disease, 4% against asymptomatic infection (AstraZeneca-Oxford Uni. vaccine, low estimate); vaccine-induced protection of 45-weeks
                     V2: vaccination with 91% VE against symptomatic disease, 56% against asymptomatic infection (AstraZeneca-Oxford Uni. vaccine, high estimate); vaccine-induced protection of 3-years"),
             size=15, x = 1, y=0.78, hjust = 1) 

p3 = cowplot::plot_grid(p3, title, 
                        ncol = 1,
                        rel_heights = c(1, 0.05) )

# plot together
cowplot::plot_grid(p1, p2, p3, 
                   ncol = 1, align = 'h', axis = 'l', 
                   labels = c("a)", "b)", "c)"),  label_size = 20,
                   rel_heights = c(0.75, 0.75, 1) )

# save
ggsave(paste0("Fig2-AZrisk-CoVAEFI_determ_", format(Sys.time(), "%d%b%Y"), "_v6.png"), path = econ_path, width = 50, height = 70/4*3, units = "cm")



########
# SFig. 2): heat map of AEFI threshold scenario analysis (absolute values)
########

# incremental net value of vaccination versus "no-vaccination"
# similar to how the CEA would have presented itself at the start
p1 = df_AEFI_incr %>%
  dplyr::filter(indicator %in% c("health benefit", "gross health benefit"), 
                #!parm == "NB_lockdown",
                !ScenA == 54321,
                delay_time == 0,
                lambda_value == 20000) %>%
  ggplot(aes(x = factor(ScenA), y= QALYs_AEFI2, fill = disc2)) + 
  geom_tile(aes(color = disc>=1), size=1.5) + #rel>=6), size=1.5) +
  scale_fill_manual(name = "health value\nof vaccine\n(compared to\nno vaccine)\nin terms of\nQALYs",
                    values=c("#DC143C", colorRampPalette(c("red","orange"))(4),
                             colorRampPalette(c("lightblue","blue"))(9) ),
                    na.value="#EEEEEE") +
  scale_colour_manual(name = "", 
                      values = c("FALSE" = "darkred",  "TRUE" = "#0080FF"),
                      labels = c("FALSE" = "negative", "TRUE" = "positive")) +
  ggh4x::facet_nested(epi_rate~parm + scenario, scales = "free",
                      labeller = as_labeller(c(label_scen)) ) +
  guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2, override.aes = list(fill = NA))) +
  theme_bw() +
  geom_point(aes(size=ifelse(observed, "dot", "no_dot")), shape = 4, show.legend = FALSE) +
  scale_size_manual(values=c(dot=4, no_dot=NA), guide="none") +
  labs(x = "rate of adverse event", y = "severity of adverse event") +
  fun_theme()


# add footnote
title <- ggdraw() + 
  draw_label(paste0("Net health values (dQALYs - dCosts/£20,000; discounted at 3.5% over 10 years) are negative in case of health losses and excess costs of vaccination compared to no-vaccination (A).
                     V1: vaccination with 62% VE against symptomatic disease, 4% against asymptomatic infection (AstraZeneca-Oxford Uni. vaccine, low estimate); vaccine-induced protection of 45-weeks
                     V2: vaccination with 91% VE against symptomatic disease, 56% against asymptomatic infection (AstraZeneca-Oxford Uni. vaccine, high estimate); vaccine-induced protection of 3-years
                     X indicate the known rates and severities of AEFIs of COVID-19 vaccines."),
             size=15, x = 1, y=0.78, hjust = 1) 

cowplot::plot_grid(p1, title, 
                   ncol = 1,
                   rel_heights = c(1, 0.05) )

# save
ggsave(paste0("FigS2-CoVAEFI_determ_", format(Sys.time(), "%d%b%Y"), "_v3.png"), path = econ_path, width = 70, height = 30, units = "cm")



########
# SFig. 3): heat map of AEFI threshold scenario analysis (incremental against no delay)
########

# incremental net value of vaccination versus "no-vaccination"
# similar to how the CEA would have presented itself at the start
df_AEFI_incr_delay = df_AEFI2 %>% ungroup() %>%
  dplyr::filter(scenario == "runA", ScenA == 0, delay_time==0) %>%
  dplyr::rename(disc_base = disc) %>%
  dplyr::select(run, parm, indicator, lambda_value, epi_rate, 
                costs_AEFI, QALYs_AEFI, disc_base) %>%
  full_join(df_AEFI2 %>% dplyr::filter(!ScenA == 0) %>% ungroup, 
            by=c("run", "parm", "indicator", 
                 "lambda_value", "epi_rate", "costs_AEFI", "QALYs_AEFI")) %>%
  dplyr::mutate(disc = ifelse(disc<0, -(disc_base - disc), disc_base - disc)) %>%
  dplyr::select(-disc_base) %>%
  dplyr::filter(scenario != "runA") %>%
  mutate(epi_rate    = factor(epi_rate, levels = order_scen),
         QALYs_AEFI2 = factor(label_AEFI[paste0(QALYs_AEFI)], levels = label_AEFI)) %>%
  # indicator of likely AEFIs for COVID vaccines (as observed in phase 3 trials and reported in March 2021)
  mutate(disc2 = cut(disc, breaks = c(-Inf, -1e+9, -1e+8, -1e+7, -1e+6, -1e+5,
                                      0, 
                                      1e+5, 1e+6, 1e+7, 1e+8, 1e+9, 1e+12, Inf), 
                     labels = as.factor(c("<-100 mln.", 
                                          "<-100 mln.", "<-10 mln.", "<-1 mln.", "<-100,000", 
                                          "<0 & >-100,000", ">0 & <100,000", 
                                          ">100,000", ">1 mln.", ">10 mln.", ">100 mln.", ">1 bln.", ">1 trln.")), 
                     include.lowest=TRUE, right=FALSE)
  ) %>%
  mutate(observed = ifelse((QALYs_AEFI %in% c(0.0015, 0.0027) & ScenA %in% paste0(seq(0.5, 0.9, 0.1))) |
                             (QALYs_AEFI == 20 & ScenA %in% c(1e-06, 1e-05)), TRUE, FALSE))


df_AEFI_incr_delay = df_AEFI2 %>% ungroup() %>%
  dplyr::filter(!scenario == "runA", 
                delay_time==0) %>%
  dplyr::rename(disc_base = disc) %>%
  dplyr::select(scenario, ScenA, 
                run, parm, indicator, lambda_value, epi_rate, 
                costs_AEFI, QALYs_AEFI, disc_base) %>%
  full_join(df_AEFI2 %>% 
              dplyr::filter(!scenario == "runA", 
                            !delay_time==54321) %>% 
              ungroup, 
            by=c("scenario", "ScenA", 
                 "run", "parm", "indicator", 
                 "lambda_value", "epi_rate", "costs_AEFI", "QALYs_AEFI")) %>%
  dplyr::mutate(ratio  = disc/disc_base,
                ratio2 = ifelse(ratio>1, -(ratio-1), 1-ratio),
                disc   = ifelse(disc<0, -(disc_base - disc), disc_base - disc)) %>%
  dplyr::select(-disc_base) %>%
  dplyr::filter(scenario != "runA") %>%
  mutate(epi_rate    = factor(epi_rate, levels = order_scen),
         delay_time  = factor(delay_time, levels = as.numeric(unique(df_AEFI2$delay_time)) %>% sort),
         QALYs_AEFI2 = factor(label_AEFI[paste0(QALYs_AEFI)], levels = label_AEFI)) %>%
  # indicator of likely AEFIs for COVID vaccines (as observed in phase 3 trials and reported in March 2021)
  mutate(disc = ifelse(disc==0, NA, disc),
         ratio3 = ratio2, #ifelse(ratio2 < (-0.03), -0.03, ratio2),
         rel2 = cut(ratio3, 
                    breaks = c(-Inf, 
                               seq(-signif(round(max(abs(range(ratio3))), digits = 3)*1000, 1)/1000,
                                   signif(round(max(abs(range(ratio3))), digits = 3)*1000, 1)/1000,
                                   0.01), 
                               Inf), 
                    labels = as.factor(c(
                      paste0("<", seq(-signif(round(max(abs(range(ratio3))), digits = 3)*1000, 1)/1000,
                                      0, 0.01)*100, "%"),
                      "<1%",
                      paste0(">", seq(0.01, signif(round(max(abs(range(ratio3))), digits = 3)*1000, 1)/1000,
                                      0.01)*100, "%")
                    )), 
                    include.lowest=TRUE, right=FALSE),
         disc2 = cut(disc, breaks = c(-Inf, -1e+9, -1e+8, -1e+7, -1e+6, -1e+5, -1e+4, -1e+3,
                                      0, 
                                      1e+3, 1e+4, 1e+5, 1e+6, 1e+7, 1e+8, 1e+9, 1e+12, Inf), 
                     labels = as.factor(c("<-100 mln.", "<-100 mln.", "<-10 mln.", "<-1 mln.", "<-100,000", "<-10,000", "<-1,000", 
                                          "<0 & >-1,000", ">0 & <1,000", 
                                          ">1,000", ">10,000", ">100,000", ">1 mln.", ">10 mln.", ">100 mln.", ">1 bln.", ">1 trln.")), 
                     include.lowest=TRUE, right=FALSE)
  ) %>%
  mutate(observed = ifelse((QALYs_AEFI %in% c(0.0015, 0.0027) & ScenA %in% paste0(seq(0.5, 0.9, 0.1))) |
                             (QALYs_AEFI == 20 & ScenA %in% c(1e-06, 1e-05)), TRUE, FALSE)) %>%
  mutate(alpha_val = ifelse(disc2 %in% c("<0 & >-1,000", ">0 & <1,000"), 0.7, 1),
         alpha_rel = ifelse(rel2 %in% c("<0%", "<1%"), 0.7, 1))


# plot
p1 = df_AEFI_incr_delay %>%
  dplyr::filter(!ScenA == 0, 
                indicator %in% c("health benefit", "gross health benefit"), 
                delay_time == 56,
                lambda_value == 20000) %>%
  ggplot(aes(x = factor(ScenA), y= QALYs_AEFI2, fill = rel2)) + 
  geom_tile(aes(#alpha = alpha_rel, 
    color = disc>=0), size=1.5) + #rel>=6), size=1.5) +
  scale_fill_manual(name = "health value\nof delay\n(compared to\nno delay)\nin terms of\nQALYs",
                    values=c("#DC143C", colorRampPalette(c("red","#FFCCCC"))(4), #"#FF9999"))(6),
                             colorRampPalette(c("#cdeefd","#03a9f4"))(7) ), #(12) ),
                    na.value="#EEEEEE") +
  scale_colour_manual(name = "", 
                      values = c("FALSE" = "darkred",  "TRUE" = "#0080FF"),
                      labels = c("FALSE" = "negative", "TRUE" = "positive")) +
  ggh4x::facet_nested(epi_rate +delay_time~parm + scenario, scales = "free",
                      labeller = as_labeller(c(label_scen)) ) +
  guides(fill = guide_legend(order = 1), alpha = FALSE, colour = guide_legend(order = 2, override.aes = list(fill = NA))) +
  theme_bw() +
  geom_point(aes(size=ifelse(observed, "dot", "no_dot")), shape = 4, show.legend = FALSE) +
  scale_size_manual(values=c(dot=4, no_dot=NA), guide="none") +
  labs(x = "rate of adverse event", y = "severity of adverse event") +
  fun_theme() +
  theme(axis.title.x   = element_text(size=24),
        axis.title.y   = element_text(size=24, vjust=3.5),
        strip.text.x = element_text(size=24),
        strip.text.y = element_text(size=24),
        legend.title = element_text(size=24),
        legend.text  = element_text(size=24),
        plot.caption = element_text(size=20),
        axis.text.x  = element_text(size=18, angle = 45, vjust = 1, hjust = 1),
        axis.text.y  = element_text(size=24))

# add footnote
title <- ggdraw() + 
  draw_label(paste0("Benefit-risk ratios are positive when the benefit (and net benefit; dQALY_loss - dCosts/£20,000) of the vaccine exceeds the (net) harm of AEFIs, and negative when the (net) harm of AEFIs exceeds the (net) benefit of the vaccine.
                     V1: vaccination with 62% VE against symptomatic disease, 4% against asymptomatic infection (AstraZeneca-Oxford Uni. vaccine, low estimate); assumed vaccine-induced protection of 45-weeks
                     V2: vaccination with 91% VE against symptomatic disease, 56% against asymptomatic infection (AstraZeneca-Oxford Uni. vaccine, high estimate); assumed vaccine-induced protection of 3-years
                     X indicate the observed rates and severities of AEFIs of COVID-19 vaccines."),
             size=18, x = 1, y=0.78, hjust = 1) 

cowplot::plot_grid(p1, title, 
                   ncol = 1,
                   rel_heights = c(1, 0.1) )



ggsave(paste0("FigS3-CoVAEFI_determ_DELAYratios_", format(Sys.time(), "%d%b%Y"), "_v2.png"), path = econ_path, width = 90, height = 40, units = "cm")

