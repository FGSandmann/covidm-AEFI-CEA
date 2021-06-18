# Value of COVID-19 vaccine
# original code of "covidm" by ND, CEA code adapted by FGS

# main functions for lockdown, epi, and econ

################################
## set-up
################################

# data wrangling for the econ analysis
library(dplyr)

# Build simpler lockdown function
library(stringr)

# lockdown function
control_lockdown = function(
  threshold_lockdown1,  # Daily new cases above which lockdown is implemented. Vector: recycled to number of populations if needed.
  contact_lockdown1,    # Contact during lockdown (e.g. c(0.8, 0.1, 0.1, 0.1) ). Same for all populations.
  threshold_lockdown2,  # Daily new cases above which lockdown is implemented. Vector: recycled to number of populations if needed.
  contact_lockdown2,    # Contact during lockdown (e.g. c(0.8, 0.1, 0.1, 0.1) ). Same for all populations.
  threshold_lockdown3,  # Daily new cases above which lockdown is implemented. Vector: recycled to number of populations if needed.
  contact_lockdown3,    # Contact during lockdown (e.g. c(0.8, 0.1, 0.1, 0.1) ). Same for all populations.
  threshold_release,    # Daily new cases below which lockdown is released. Recycled if necessary.
  contact_release,      # Contact during lockdown release (e.g. c(1.0, 0.8, 0.8, 0.8))
  contact_holidays,     # Contact rates during school holidays
  t_noFurtherLockdowns, # no further lockdowns after this day
  t_1,   t_2,  t_3,  t_4,  t_5,  t_6,  t_7,  t_8,  t_9, t_10,  # Time of ignoring/imposing lockdowns (school holidays)
  t_11, t_12, t_13, t_14, t_15, t_16, t_17, t_18, t_19, t_20,  # hard coded for 10 years as observer not working well with schedule
  t_21, t_22, t_23, t_24, t_25, t_26, t_27, t_28, t_29, t_30, 
  t_31, t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, t_40
) 
{
  as_init_list = function(x) paste0("{ ", paste0(x, collapse = ", "), " }");
  str_interp('
    const int obsStage = 0;
    const int obsR0 = 1;
    const int obsRt = 2;
    const vector<double> threshold_lockdown1 = { ${ paste0(threshold_lockdown1, collapse = ", ") } };
    const vector<double> threshold_lockdown2 = { ${ paste0(threshold_lockdown2, collapse = ", ") } };
    const vector<double> threshold_lockdown3 = { ${ paste0(threshold_lockdown3, collapse = ", ") } };
    const vector<double> threshold_release   = { ${ paste0(threshold_release,   collapse = ", ") } };
    const vector<double> t_noFurtherLockdowns  = { ${ paste0(t_noFurtherLockdowns,  collapse = ", ") } };
    const vector<double> t_1  = { ${ paste0(t_1,  collapse = ", ") } };
    const vector<double> t_2  = { ${ paste0(t_2,  collapse = ", ") } };
    const vector<double> t_3  = { ${ paste0(t_3,  collapse = ", ") } };
    const vector<double> t_4  = { ${ paste0(t_4,  collapse = ", ") } };
    const vector<double> t_5  = { ${ paste0(t_5,  collapse = ", ") } };
    const vector<double> t_6  = { ${ paste0(t_6,  collapse = ", ") } };
    const vector<double> t_7  = { ${ paste0(t_7,  collapse = ", ") } };
    const vector<double> t_8  = { ${ paste0(t_8,  collapse = ", ") } };
    const vector<double> t_9  = { ${ paste0(t_9,  collapse = ", ") } };
    const vector<double> t_10 = { ${ paste0(t_10, collapse = ", ") } };
    const vector<double> t_11 = { ${ paste0(t_11, collapse = ", ") } };
    const vector<double> t_12 = { ${ paste0(t_12, collapse = ", ") } };
    const vector<double> t_13 = { ${ paste0(t_13, collapse = ", ") } };
    const vector<double> t_14 = { ${ paste0(t_14, collapse = ", ") } };
    const vector<double> t_15 = { ${ paste0(t_15, collapse = ", ") } };
    const vector<double> t_16 = { ${ paste0(t_16, collapse = ", ") } };
    const vector<double> t_17 = { ${ paste0(t_17, collapse = ", ") } };
    const vector<double> t_18 = { ${ paste0(t_18, collapse = ", ") } };
    const vector<double> t_19 = { ${ paste0(t_19, collapse = ", ") } };
    const vector<double> t_20 = { ${ paste0(t_20, collapse = ", ") } };
    const vector<double> t_21 = { ${ paste0(t_21, collapse = ", ") } };
    const vector<double> t_22 = { ${ paste0(t_22, collapse = ", ") } };
    const vector<double> t_23 = { ${ paste0(t_23, collapse = ", ") } };
    const vector<double> t_24 = { ${ paste0(t_24, collapse = ", ") } };
    const vector<double> t_25 = { ${ paste0(t_25, collapse = ", ") } };
    const vector<double> t_26 = { ${ paste0(t_26, collapse = ", ") } };
    const vector<double> t_27 = { ${ paste0(t_27, collapse = ", ") } };
    const vector<double> t_28 = { ${ paste0(t_28, collapse = ", ") } };
    const vector<double> t_29 = { ${ paste0(t_29, collapse = ", ") } };
    const vector<double> t_30 = { ${ paste0(t_30, collapse = ", ") } };
    const vector<double> t_31 = { ${ paste0(t_31, collapse = ", ") } };
    const vector<double> t_32 = { ${ paste0(t_32, collapse = ", ") } };
    const vector<double> t_33 = { ${ paste0(t_33, collapse = ", ") } };
    const vector<double> t_34 = { ${ paste0(t_34, collapse = ", ") } };
    const vector<double> t_35 = { ${ paste0(t_35, collapse = ", ") } };
    const vector<double> t_36 = { ${ paste0(t_36, collapse = ", ") } };
    const vector<double> t_37 = { ${ paste0(t_37, collapse = ", ") } };
    const vector<double> t_38 = { ${ paste0(t_38, collapse = ", ") } };
    const vector<double> t_39 = { ${ paste0(t_39, collapse = ", ") } };
    const vector<double> t_40 = { ${ paste0(t_40, collapse = ", ") } };

    auto set_contact = [&](unsigned int p, vector<double> con, int stage) {
        P.pop[p].contact = con;
        P.pop[p].needs_recalc = true;
        P.pop[p].Recalculate();
        dyn.Obs(t, p, obsStage, 0) = stage;
    };

    for (unsigned int p = 0; p < P.pop.size(); ++p)
    {
    
    double sum_pop = std::accumulate(P.pop[p].size.begin(), P.pop[p].size.end(), 0.0);
    double curr_incidence  = (dyn("cases", t, { p }, {})     / sum_pop *100000);

        if (
        ((t >= t_1[p % t_1.size()]   && t < t_2[p % t_2.size()])   && (curr_incidence < threshold_lockdown3[p % threshold_lockdown3.size()])) ||
        ((t >= t_3[p % t_3.size()]   && t < t_4[p % t_4.size()])   && (curr_incidence < threshold_lockdown3[p % threshold_lockdown3.size()])) ||
        ((t >= t_5[p % t_5.size()]   && t < t_6[p % t_6.size()])   && (curr_incidence < threshold_lockdown3[p % threshold_lockdown3.size()])) ||
        ((t >= t_7[p % t_7.size()]   && t < t_8[p % t_8.size()])   && (curr_incidence < threshold_lockdown3[p % threshold_lockdown3.size()])) ||
        ((t >= t_9[p % t_9.size()]   && t < t_10[p % t_10.size()]) && (curr_incidence < threshold_lockdown3[p % threshold_lockdown3.size()])) ||
        ((t >= t_11[p % t_11.size()] && t < t_12[p % t_12.size()]) && (curr_incidence < threshold_lockdown3[p % threshold_lockdown3.size()])) ||
        ((t >= t_13[p % t_13.size()] && t < t_14[p % t_14.size()]) && (curr_incidence < threshold_lockdown3[p % threshold_lockdown3.size()])) ||
        ((t >= t_15[p % t_15.size()] && t < t_16[p % t_16.size()]) && (curr_incidence < threshold_lockdown3[p % threshold_lockdown3.size()])) ||
        ((t >= t_17[p % t_17.size()] && t < t_18[p % t_18.size()]) && (curr_incidence < threshold_lockdown3[p % threshold_lockdown3.size()])) ||
        ((t >= t_19[p % t_19.size()] && t < t_20[p % t_20.size()]) && (curr_incidence < threshold_lockdown3[p % threshold_lockdown3.size()])) ||
        ((t >= t_21[p % t_21.size()] && t < t_22[p % t_22.size()]) && (curr_incidence < threshold_lockdown3[p % threshold_lockdown3.size()])) ||
        ((t >= t_23[p % t_23.size()] && t < t_24[p % t_24.size()]) && (curr_incidence < threshold_lockdown3[p % threshold_lockdown3.size()])) ||
        ((t >= t_25[p % t_25.size()] && t < t_26[p % t_26.size()]) && (curr_incidence < threshold_lockdown3[p % threshold_lockdown3.size()])) ||
        ((t >= t_27[p % t_27.size()] && t < t_28[p % t_28.size()]) && (curr_incidence < threshold_lockdown3[p % threshold_lockdown3.size()])) ||
        ((t >= t_29[p % t_29.size()] && t < t_30[p % t_30.size()]) && (curr_incidence < threshold_lockdown3[p % threshold_lockdown3.size()])) ||
        ((t >= t_31[p % t_31.size()] && t < t_32[p % t_32.size()]) && (curr_incidence < threshold_lockdown3[p % threshold_lockdown3.size()])) ||
        ((t >= t_33[p % t_33.size()] && t < t_34[p % t_34.size()]) && (curr_incidence < threshold_lockdown3[p % threshold_lockdown3.size()])) ||
        ((t >= t_35[p % t_35.size()] && t < t_36[p % t_36.size()]) && (curr_incidence < threshold_lockdown3[p % threshold_lockdown3.size()])) ||
        ((t >= t_37[p % t_37.size()] && t < t_38[p % t_38.size()]) && (curr_incidence < threshold_lockdown3[p % threshold_lockdown3.size()])) ||
        ((t >= t_39[p % t_39.size()] && t < t_40[p % t_40.size()]) && (curr_incidence < threshold_lockdown3[p % threshold_lockdown3.size()]))
        )
            set_contact(p, ${ as_init_list(contact_holidays) }, 5);
            
            
       else if (t > t_noFurtherLockdowns[p % t_noFurtherLockdowns.size()])
       {
       set_contact(p, ${ as_init_list(contact_release)  }, 4);
       }
    
        else if (t > 0)
        {
            int prev_stage = dyn.Obs(t - 1, p, obsStage, 0);
            
        
            double sum_pop = std::accumulate(P.pop[p].size.begin(), P.pop[p].size.end(), 0.0);
            double curr_incidence  = (dyn("cases", t, { p }, {})  / sum_pop *100000);
            
            if (prev_stage==0 && curr_incidence >= threshold_lockdown1[p % threshold_lockdown1.size()])
                set_contact(p, ${ as_init_list(contact_lockdown1) }, 1);
            else if (prev_stage!=0 && prev_stage!=1 && curr_incidence < threshold_release[p % threshold_release.size()])
                set_contact(p, ${ as_init_list(contact_release)  }, 4);
            else if (prev_stage!=0 && prev_stage!=1 && (curr_incidence >= threshold_lockdown3[p % threshold_lockdown3.size()]))
                set_contact(p, ${ as_init_list(contact_lockdown3) }, 3);
            else if (prev_stage!=0 && prev_stage!=1 && curr_incidence >= threshold_lockdown2[p % threshold_lockdown2.size()] )
                set_contact(p, ${ as_init_list(contact_lockdown2) }, 2);
            else if (prev_stage == 5)
                set_contact(p, ${ as_init_list(contact_release) }, 4);
            else
                dyn.Obs(t, p, obsStage, 0) = prev_stage;
        }
        
        dyn.Obs(t, p, obsR0, 0) = estimate_R0(P, t, p, 20);
        dyn.Obs(t, p, obsRt, 0) = estimate_Rt(P, dyn, t, p, 20);
    }');
}


##############################
## B) Data for econ model
##############################

# probabilistic parameters (method of modes)
beta_dis_alpha <- function(mean, sd){ mean*(((mean*(1-mean))/(sd^2))-1) }
beta_dis_beta  <- function(mean, sd){ (1-mean)*(((mean*(1-mean))/(sd^2))-1) }

gamma_dis_shape <- function(mean, low, up){(mean/sqrt((up-low)/(2*1.96)) )^2 }
gamma_dis_scale <- function(mean, low, up){(sqrt((up-low)/(2*1.96))^2)/mean }

lognorm_location <- function(mean, sd){ log(mean^2 / sqrt(sd^2 + mean^2)) }
lognorm_shape    <- function(mean, sd){ sqrt(log(1 + (sd^2 / mean^2))) }


##############################
## C) Set-up epi model (SEEIIR)
##############################

#scenario A: No vaccination

#scenario B: vaccinate ages 15+ years, pessimistic VE and 45 wks vaccine protection

#scenario C: vaccinate ages 15+ years, optimistic VE and 3 yrs vaccine protection

#In addition, all scenarios considered natural waning of 45 wks; no contacts during school holiday
#Model run with 5 new infections each day for 7 days every month in first 6 months (for seeding)

# model set up for baseline A plus scenarios B and C (allowing differences per vaccine that are mostly equal currently)
epi_mod <- function(para                = params,
                    waning_nat          = 45*7,          # natural waning of common coronaviruses in days
                    includeC            = TRUE,          # 2 vaccination alternatives, B and C
                    vacc_evenly         = FALSE,         # if false, trying to resemble current JCVI priority groups (as of Nov 2020)
                    n_init_vaccinees_B  = 100000,        # 4000000/31
                    n_vaccinees_B       = 200000,        # based on 40 mio doses by end March; https://www.bbc.co.uk/news/uk-55096434
                    n_revac_vaccinees_B = 100000,        # revaccinate
                    n_init_vaccinees_C  = 100000,        # 4000000/31
                    n_vaccinees_C       = 200000,        # based on 40 mio doses by end March; https://www.bbc.co.uk/news/uk-55096434
                    n_revac_vaccinees_C = 100000,        # revaccinate
                    start_vacc_B        = as.Date("2020-12-08")+28, # start on 08-12-2020, but protection after 28 days
                    start_vacc_C        = as.Date("2020-12-08")+28, # start on 08-12-2020, but protection after 28 days
                    vacc_minAge_B       = 15,            # vaccinating from age 15
                    vacc_maxAge_B       = 100,           # vaccinating all 15+
                    vacc_minAge_C       = 15,            # vaccinating from age 15
                    vacc_maxAge_C       = 100,           # vaccinating all 15+
                    wane_vacc_B         = 45*7,          # sub-optimal: like natural waning
                    wane_vacc_C         = 365.25*3,      # optimal: 3 years 
                    #VE_B                = 0.50,          # sub-optimal: 30% # don't use any longer; changed to "take"/probability of moving to V
                    #VE_C                = 0.95,          # optimal: 95%     # don't use any longer; changed to "take"/probability of moving to V
                    VE_B_infection      = rep(0, 16),
                    VE_B_disease        = rep(0.50, 16), 
                    VE_C_infection      = rep(0.95, 16), # based on Pfizer vacc; https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/940565/Information_for_Healthcare_Professionals_on_Pfizer_BioNTech_COVID-19_vaccine.pdf
                    VE_C_disease        = rep(0,    16), # set to 0 as otherwise VE additive and not reflecting phase 3 trials; making it more complicated although no massive difference for C
                    cov1_vacc_O50_B     = 0.90,          # assume 75% initially
                    cov1_vacc_U50_B     = 0.75,          # assume 75% initially
                    cov2_vacc_B         = c(0, rep(0.5, 12), rep(0.75, 3)), # assume 50% and 75% revaccinated similar to flu coverage
                    cov1_vacc_O50_C     = 0.90,          # assume 75% initially
                    cov1_vacc_U50_C     = 0.75,          # assume 75% initially
                    cov2_vacc_C         = c(0, rep(0.5, 12), rep(0.75, 3)), # assume 50% and 75% revaccinated similar to flu coverage
                    vacc_delay_prop     = NULL,          # NEW: proportion reduction in available doses in ages <50
                    vacc_delay_time     = NULL,          # NEW: extra time to vaccinate ages <50
                    n_iter              = 5,             # number of iterations run the model
                    seed_val            = 0,             # specify random seed
                    run_determ          = TRUE,          # set stochastic or deterministic
                    storeMinParms       = TRUE){         # keep only compartments needed for epi/econ analysis

  # if deterministic, set n_iter to 1
  if(run_determ==TRUE){
    n_iter <- 1
    para$deterministic = TRUE
  } else {
    para$deterministic = FALSE
  }
  
  
  ##############################
  # (A) no vaccination
  ##############################

  # There are 3 parameters for controlling vaccination.
  # v -- number of vaccines administered per day, for each age group.
  # ev -- vaccine effectiveness, assuming all-or-nothing protection. (So effectively, v[i] * ev[i] people get vaccine protection each day in age group i.)
  # wv -- rate at which vaccine protection wanes (per day), for each age group.

  n_age_groups     = length(para$pop[[1]]$size)
  para$pop[[1]]$v  = rep(0, n_age_groups) # no vaccines administered
  para$pop[[1]]$wv = rep(0, n_age_groups) # no vaccine protection
  para$pop[[1]]$wn = rep((1/waning_nat), n_age_groups) # # natural waning of 45 weeks on average
  
  runA <- cm_simulate(para, model_seed = seed_val, n_iter)
  
  # nr of vaccinated individuals (keep for econ analysis later)
  df_vaccpts <- do.call(rbind,
                        rep(list(rep(0, n_age_groups)), 
                            length(unique(runA$dynamics$t)) ))
  
  df_vaccpts <-  as.data.table(df_vaccpts)
  
  # replicate by number of iterations
  df_vaccpts <- do.call("rbind", rep(list(df_vaccpts), n_iter))
  
  df_vaccpts$t <- rep( c(0:(max(unique(runA$dynamics$t))) ), n_iter)
  df_vaccpts <- df_vaccpts %>%
    setNames(c(para$pop[[1]]$group_names, "t")) %>%
    tidyr::gather(group, value, -t)
  df_vaccpts$run <- rep(1:n_iter, each=length(unique(runA$dynamics$t)) )
  df_vaccpts$compartment="vaccinees"
  df_vaccpts$population <- para$pop[[1]]$name
  
  # merge
  runA$dynamics <- rbind(runA$dynamics, df_vaccpts)
  
  
  ####################################################
  # (B) pessimistic vaccination
  ####################################################
  
  # running epi model separately for B and C works well as seed re-set to same value each time
  n_age_gr_vacc_B   <- cm_age_coefficients(vacc_minAge_B, vacc_maxAge_B, 5 * (0:length(para$pop[[1]]$size)))
  n_age_gr_ttl_B    <- sum(para$pop[[1]]$size[which(cm_age_coefficients(vacc_minAge_B, vacc_maxAge_B, 5 * (0:length(para$pop[[1]]$size)))==1)])
  
  # revisions:
  #ev = probability that a person who has received a V vaccine actually moves into V, otherwise they stay in S, so it's like a "take" vaccine effectivemess.
  #ei_v = effectiveness of V vaccine against infection, given that a person is in compartment V.
  #ed_vi = effectiveness of V vaccine against disease, given that a person in compartment V got infected.
  #wv and wv2 are the waning rate of V and V2 back to S.
  para$pop[[1]]$ev    = rep(1, 16)
  para$pop[[1]]$ei_v  = VE_B_infection #rep(0, 16) #
  para$pop[[1]]$ed_vi = VE_B_disease #rep(0.5, 16) #
  para$pop[[1]]$wv    = rep((1/wane_vacc_B), n_age_groups) # vaccine protection for x years on average
  
  # numbers of vaccinees
  # uniformly vaccinating; no longer used other than sensitivity analyses
  if(vacc_evenly==TRUE){
    
    vacc_vals_B <- list(rep(0, n_age_groups), # initially no vaccine available
                        # total nr of vaccinees distributed evenly (assumed)
                        c(n_vaccinees_B / sum(n_age_gr_vacc_B) * n_age_gr_vacc_B), # could redistribute more to certain ages here
                        # total nr of vaccinees after having vaccinated all initially (accounting for coverage)
                        c(n_age_gr_ttl_B *cov2_vacc_B /365.25 /sum(n_age_gr_vacc_B) * n_age_gr_vacc_B) ) #re-vaccinate
    
    # timing of initial vaccination and revaccination
    vacc_times_B <- c(0, 
                      # start vaccinating
                      as.numeric(as.Date(start_vacc_B) - as.Date(para$date0) ),
                      # how long initially when vaccinating x a day (all reached, accounting for coverage)
                      ceiling(as.numeric(as.Date(start_vacc_B) - as.Date(para$date0) ) + 
                                n_age_gr_ttl_B *cov1_vacc_B /n_vaccinees_B) ) # if want to use this by different VE, would need to adjust over/under 50..
    
  } else {
  # targeted vaccinating in line with provisional vaccination priorities
      
    # moderate to high risk (to do: move this out of function to be customisable; or wait for revision of new compartments)
    #high risk adults under 65 years of age; Clarke et al, females and males in the UK
    n_highRisk  = c(40, 570, 1038, 13140, 27014, 74912, 108766, 144361, 144723, 201420, 270600, 371662, 368724, 413316, 479443, 1178643)
    n_highRisk   = n_highRisk * cm_age_coefficients(0, 65, 5 * (0:length(para$pop[[1]]$size)))
    #moderate/increased risk adults under 65 years of age; Clarke et al, females and males in the UK
    n_modRisk  = c(73604, 118254, 131285, 206125, 366435, 551319, 725170, 827328, 965984, 1177495, 1563243, 1930802, 2011177, 2045102, 2332117, 4746535)
    n_modRisk   = n_modRisk * cm_age_coefficients(0, 65, 5 * (0:length(para$pop[[1]]$size)))
    
    # work out proportions of who received the daily number of vaccines; priority order of JCVI
    #https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/939119/Greenbook_chapter_14a___provisional_guidance_subject_to_MHRA_approval_of_vaccine_supply_.pdf
    vacc_vals_B <-list()
    # initially no vaccine available
    vacc_vals_B[[1]] = rep(0, n_age_groups)
    
    # start with prop in care homes and social care first and be done in 4 weeks, then move on to 75+ (as there are no)
    #291,000 care home residents (260k 75+, 31k 65-74) 
    #1.52M social care  (20-64 y)
    #1.31M million hospital staff (20-64 y)
    n_HCW = 1520000+1310000 # 1.52M social care staff, 1.31M million hospital staff (20-64 y)
    # vaccinate care home pop
    ch_pop = c((31000/2), (31000/2), (260000))
    # vaccinate HCW of working age and CH residents
    pop_v1 = n_HCW / 
      sum(cm_age_coefficients(20, 65, 5*(0:length(para$pop[[1]]$size)))) * 
      cm_age_coefficients(20, 65, 5 * (0:length(para$pop[[1]]$size)))
    pop_v1[14:16] = ch_pop
    # check if aligning with vaccination age groups
    pop_v1 = pop_v1 * cm_age_coefficients(vacc_minAge_B, vacc_maxAge_B, 5 * (0:length(para$pop[[1]]$size)))
    vacc_vals_B[[2]] = round(n_init_vaccinees_B / sum(pop_v1) * pop_v1)
    
    #75+, (once care workers exhausted move to remaining 75+ only)
    vacc_vals_B[[3]] = n_vaccinees_B * cm_age_coefficients(75, 100, 5*(0:length(para$pop[[1]]$size)))
    # check if aligning with vaccination age groups
    vacc_vals_B[[3]] = vacc_vals_B[[3]] * cm_age_coefficients(vacc_minAge_B, vacc_maxAge_B, 5 * (0:length(para$pop[[1]]$size)))
    
    #70+ and high-risk adults under 65 years of age; Clarke et al, females and males in the UK
    # ages 70-74
    n_highRisk[which(cm_age_coefficients(70, 75, 5 * (0:length(para$pop[[1]]$size)))==1)] = 
      para$pop[[1]]$size[which(cm_age_coefficients(70, 75, 5 * (0:length(para$pop[[1]]$size)))==1)]-
      31000/2 # minus care home residents
    # check if aligning with vaccination age groups
    n_highRisk70 = n_highRisk * cm_age_coefficients(vacc_minAge_B, vacc_maxAge_B, 5 * (0:length(para$pop[[1]]$size)))
    vacc_vals_B[[4]] = round(n_vaccinees_B / sum(n_highRisk70) * n_highRisk70)
    if(any(is.na(vacc_vals_B[[4]] )==TRUE)) vacc_vals_B[[4]] = rep(0, n_age_groups)
    
    #65+
    vacc_vals_B[[5]] = n_vaccinees_B * cm_age_coefficients(65, 70, 5*(0:length(para$pop[[1]]$size)))
    # check if aligning with vaccination age groups
    vacc_vals_B[[5]] = vacc_vals_B[[5]] * cm_age_coefficients(vacc_minAge_B, vacc_maxAge_B, 5 * (0:length(para$pop[[1]]$size)))
    
    #moderate/increased risk adults under 65 years of age; Clarke et al, females and males in the UK
    # check if aligning with vaccination age groups
    n_modRisk = n_modRisk * cm_age_coefficients(vacc_minAge_B, vacc_maxAge_B, 5 * (0:length(para$pop[[1]]$size)))
    vacc_vals_B[[6]] = round(n_vaccinees_B / sum(n_modRisk) * n_modRisk)
    if(any(is.na(vacc_vals_B[[6]] )==TRUE)) vacc_vals_B[[6]] = rep(0, n_age_groups)
    
    #60+,
    vacc_vals_B[[7]] = n_vaccinees_B * cm_age_coefficients(60, 65, 5*(0:length(para$pop[[1]]$size)))
    # check if aligning with vaccination age groups
    vacc_vals_B[[7]] = vacc_vals_B[[7]] * cm_age_coefficients(vacc_minAge_B, vacc_maxAge_B, 5 * (0:length(para$pop[[1]]$size)))
    
    #55+, 
    vacc_vals_B[[8]] = n_vaccinees_B * cm_age_coefficients(55, 60, 5*(0:length(para$pop[[1]]$size)))
    # check if aligning with vaccination age groups
    vacc_vals_B[[8]] = vacc_vals_B[[8]] * cm_age_coefficients(vacc_minAge_B, vacc_maxAge_B, 5 * (0:length(para$pop[[1]]$size)))
    
    #50+, 
    vacc_vals_B[[9]] = n_vaccinees_B * cm_age_coefficients(50, 55, 5*(0:length(para$pop[[1]]$size)))
    # check if aligning with vaccination age groups
    vacc_vals_B[[9]] = vacc_vals_B[[9]] * cm_age_coefficients(vacc_minAge_B, vacc_maxAge_B, 5 * (0:length(para$pop[[1]]$size)))
    
    #rest (vacc_minAge_B-49)
    vacc_vals_B[[10]] = n_vaccinees_B / 
      sum(cm_age_coefficients(vacc_minAge_B, 50, 5*(0:length(para$pop[[1]]$size)))) * 
      cm_age_coefficients(vacc_minAge_B, 50, 5 * (0:length(para$pop[[1]]$size)))
    # check if aligning with vaccination age groups
    vacc_vals_B[[10]] = vacc_vals_B[[10]] * cm_age_coefficients(vacc_minAge_B, vacc_maxAge_B, 5 * (0:length(para$pop[[1]]$size)))
    if(any(is.na(vacc_vals_B[[10]] )==TRUE)) vacc_vals_B[[10]] = rep(0, n_age_groups)
    
    # revaccinate; total nr of vaccinees after having vaccinated all initially (accounting for coverage)
    # only need to revaccinate as many as want to (based on coverage..), or at assumed daily capacity
    vacc_vals_B[[11]] = min(n_revac_vaccinees_B, c(sum(para$pop[[1]]$size *cov2_vacc_B) /365.25)) / 
      sum(n_age_gr_vacc_B) * n_age_gr_vacc_B
    # check if aligning with vaccination age groups
    vacc_vals_B[[11]] = vacc_vals_B[[11]] * cm_age_coefficients(vacc_minAge_B, vacc_maxAge_B, 5 * (0:length(para$pop[[1]]$size)))
    
    
    # start simulation
    vacc_t0 = 0
    
    # vaccination introduction
    vacc_t1 = ceiling(vacc_t0 + as.numeric(as.Date(start_vacc_B) - as.Date(para$date0) ))
    
    # care homes and HCW vaccinated in first ~4 weeks
    vacc_t2 = pop_v1* cm_age_coefficients(vacc_minAge_B, vacc_maxAge_B, 5 * (0:length(para$pop[[1]]$size)))
    vacc_t2 = ceiling(vacc_t1 + sum(vacc_t2*cov1_vacc_O50_B)/n_init_vaccinees_B)
    
    #75+ minus care home residents
    n_vacc = rep(0, n_age_groups)
    # assign 75+
    n_vacc[which(cm_age_coefficients(75, 100, 5 * (0:length(para$pop[[1]]$size)))==1)] = 
      para$pop[[1]]$size[which(cm_age_coefficients(75, 100, 5 * (0:length(para$pop[[1]]$size)))==1)]-260000 # minus healt care residents
    # check if aligning with vaccination age groups
    n_vacc = n_vacc * cm_age_coefficients(vacc_minAge_B, vacc_maxAge_B, 5 * (0:length(para$pop[[1]]$size)))
    # calculate vaccination rate
    vacc_t3 = ceiling(vacc_t2 + sum(n_vacc*cov1_vacc_O50_B)/n_vaccinees_B)
    
    #70+ and high-risk <65 minus care home residents (checked and subtracted within n_highRisk70 already)
    vacc_t4 = ceiling(vacc_t3 + sum(n_highRisk70*cov1_vacc_O50_B) / n_vaccinees_B)
    
    #65+ minus care home residents
    n_vacc = rep(0, n_age_groups)
    n_vacc[which(cm_age_coefficients(65, 70, 5 * (0:length(para$pop[[1]]$size)))==1)] = 
      para$pop[[1]]$size[which(cm_age_coefficients(65, 70, 5 * (0:length(para$pop[[1]]$size)))==1)]-31000/2 # minus care home residents
    # check if aligning with vaccination age groups
    n_vacc = n_vacc * cm_age_coefficients(vacc_minAge_B, vacc_maxAge_B, 5 * (0:length(para$pop[[1]]$size)))
    
    # calculate vaccination rate
    vacc_t5 = ceiling(vacc_t4 + sum(n_vacc*cov1_vacc_O50_B)/n_vaccinees_B)
    
    #moderate risk <65 (checked within n_modRisk already)
    vacc_t6 = ceiling(vacc_t5 + sum(n_modRisk*cov1_vacc_O50_B) / n_vaccinees_B)
    
    #60+ minus moderate-high risk
    n_vacc = rep(0, n_age_groups)
    n_vacc[which(cm_age_coefficients(60, 65, 5 * (0:length(para$pop[[1]]$size)))==1)] = 
      para$pop[[1]]$size[which(cm_age_coefficients(60, 65, 5 * (0:length(para$pop[[1]]$size)))==1)]-
      # minus high- and moderate risk
      n_highRisk[which(cm_age_coefficients(60, 65, 5 * (0:length(para$pop[[1]]$size)))==1)]-
      n_modRisk[which(cm_age_coefficients(60, 65, 5 * (0:length(para$pop[[1]]$size)))==1)]- 
      # minus HCW
      n_HCW/sum(cm_age_coefficients(20, 65, 5 * (0:length(para$pop[[1]]$size)))) 
    # check if aligning with vaccination age groups
    n_vacc = n_vacc * cm_age_coefficients(vacc_minAge_B, vacc_maxAge_B, 5 * (0:length(para$pop[[1]]$size)))
    # calculate vaccination rate
    vacc_t7 = ceiling(vacc_t6 + sum(n_vacc*cov1_vacc_O50_B)/n_vaccinees_B)
    
    #55+ minus moderate-high risk 
    n_vacc = rep(0, n_age_groups)
    n_vacc[which(cm_age_coefficients(55, 60, 5 * (0:length(para$pop[[1]]$size)))==1)] = 
      para$pop[[1]]$size[which(cm_age_coefficients(55, 60, 5 * (0:length(para$pop[[1]]$size)))==1)]-
      # minus high- and moderate risk
      n_highRisk[which(cm_age_coefficients(55, 60, 5 * (0:length(para$pop[[1]]$size)))==1)]-
      n_modRisk[which(cm_age_coefficients(55, 60, 5 * (0:length(para$pop[[1]]$size)))==1)]- 
      # minus HCW
      n_HCW/sum(cm_age_coefficients(20, 65, 5 * (0:length(para$pop[[1]]$size))))
    # check if aligning with vaccination age groups
    n_vacc = n_vacc * cm_age_coefficients(vacc_minAge_B, vacc_maxAge_B, 5 * (0:length(para$pop[[1]]$size)))
    # calculate vaccination rate
    vacc_t8 = ceiling(vacc_t7 + sum(n_vacc*cov1_vacc_O50_B)/n_vaccinees_B)
    
    #50+ minus moderate-high risk 
    n_vacc = rep(0, n_age_groups)
    n_vacc[which(cm_age_coefficients(50, 55, 5 * (0:length(para$pop[[1]]$size)))==1)] = 
      para$pop[[1]]$size[which(cm_age_coefficients(50, 55, 5 * (0:length(para$pop[[1]]$size)))==1)]-
      # minus high- and moderate risk
      n_highRisk[which(cm_age_coefficients(50, 55, 5 * (0:length(para$pop[[1]]$size)))==1)]-
      n_modRisk[which(cm_age_coefficients(50, 55, 5 * (0:length(para$pop[[1]]$size)))==1)]- 
      # minus HCW
      n_HCW/sum(cm_age_coefficients(20, 65, 5 * (0:length(para$pop[[1]]$size))))
    # check if aligning with vaccination age groups
    n_vacc = n_vacc * cm_age_coefficients(vacc_minAge_B, vacc_maxAge_B, 5 * (0:length(para$pop[[1]]$size)))
    # calculate vaccination rate
    vacc_t9 = ceiling(vacc_t8 + sum(n_vacc*cov1_vacc_O50_B)/n_vaccinees_B)
    
    #rest (vacc_minAge_B-49) minus moderate-high risk; only 50% #that's the time after which we start re-vaccinating
    n_vacc = rep(0, n_age_groups)
    n_vacc[which(cm_age_coefficients(vacc_minAge_B, 50, 5 * (0:length(para$pop[[1]]$size)))==1)] = 
      para$pop[[1]]$size[which(cm_age_coefficients(vacc_minAge_B, 50, 5 * (0:length(para$pop[[1]]$size)))==1)]-
      # minus high- and moderate risk
      n_highRisk[which(cm_age_coefficients(vacc_minAge_B, 50, 5 * (0:length(para$pop[[1]]$size)))==1)]-
      n_modRisk[which(cm_age_coefficients(vacc_minAge_B, 50, 5 * (0:length(para$pop[[1]]$size)))==1)]- 
      # minus HCW
      n_HCW/
       sum(cm_age_coefficients(20, 65, 5 * (0:length(para$pop[[1]]$size))))*
       sum(cm_age_coefficients(20, 50, 5 * (0:length(para$pop[[1]]$size))))
    # check if aligning with vaccination age groups
    n_vacc = n_vacc * cm_age_coefficients(vacc_minAge_B, vacc_maxAge_B, 5 * (0:length(para$pop[[1]]$size)))
    
    # calculate vaccination rate (accounting for possible delays of AZ vaccine)
    
    # implementation of delay in programme by x days
    #"This safety update report is based on detailed analysis of data up to 28 April 2021. 
    #At this date, an estimated 11.4 million first doses of the Pfizer/BioNTech vaccine and 22.6 million first doses
    #of the COVID-19 Vaccine AstraZeneca vaccine had been administered, and around 8.1 million and 5.9 million 
    #second doses of the Pfizer/BioNTech vaccine and COVID-19 Vaccine AstraZeneca vaccine respectively. 
    #An approximate 0.1 million first doses of the Moderna vaccine have also now been administered."
    #    https://www.gov.uk/government/publications/coronavirus-covid-19-vaccine-adverse-reactions/coronavirus-vaccine-summary-of-yellow-card-reporting
    #    # only 40% of usual vaccines available if AZ stopped; (22.6+5.9) / (11.4+22.6+8.1+5.9)
    
    # implementation either via delay time or via doses
    # drop of all AZ vaccine would be equivalent to roughly 2 months delay
    #sum(n_vacc*cov1_vacc_U50_B)/n_vaccinees_B_delay
    #((sum(n_vacc*cov1_vacc_U50_B)/n_vaccinees_B) +60)
    
    # adjust schedule and timing due to delays
    # original function (if nothing delayed)
    if(is.null(vacc_delay_prop) & is.null(vacc_delay_time)){
      vacc_t99 = ceiling(vacc_t9 + sum(n_vacc*cov1_vacc_U50_B)/n_vaccinees_B )
    }
    
    # if prop drop in vacc doses given
    if(!is.null(vacc_delay_prop) & is.null(vacc_delay_time)){
      #new time
      vacc_t99 = ceiling(vacc_t9 + sum(n_vacc*cov1_vacc_U50_B)/ (n_vaccinees_B*vacc_delay_prop) )
      #new dosing
      vacc_vals_B[[10]][ which(vacc_vals_B[[10]]>0) ] = 
        vacc_vals_B[[10]][ which(vacc_vals_B[[10]]>0) ] /
        vacc_vals_B[[10]][ which(vacc_vals_B[[10]]>0) ] * 
        (n_vaccinees_B*vacc_delay_prop)/length(which(vacc_vals_B[[10]]>0) )
      
    }
    
    # if delay time given
    if(is.null(vacc_delay_prop) & !is.null(vacc_delay_time)){
      
      # new nr of doses
      n_vaccinees_B_delay = 
        n_vaccinees_B*(sum(n_vacc*cov1_vacc_U50_B)/n_vaccinees_B) / 
        ((sum(n_vacc*cov1_vacc_U50_B)/n_vaccinees_B) +vacc_delay_time)
      
      #new time
      vacc_t99 = ceiling(vacc_t9 + sum(n_vacc*cov1_vacc_U50_B)/ n_vaccinees_B_delay )
      #new dosing
      vacc_vals_B[[10]][ which(vacc_vals_B[[10]]>0) ] = 
        vacc_vals_B[[10]][ which(vacc_vals_B[[10]]>0) ] /
        vacc_vals_B[[10]][ which(vacc_vals_B[[10]]>0) ] * 
        (n_vaccinees_B_delay/length(which(vacc_vals_B[[10]]>0) ))
      
    }
    
    # stop if both given.. could optimise min(vacc_delay_prop) and max(vacc_delay_time)
    if(!is.null(vacc_delay_prop) & !is.null(vacc_delay_time)){
      stop("Both delay time and proportion reduction given; not implemented currently")
    }

    # get all times of vaccination changes
    vacc_times_B <- c(vacc_t0,vacc_t1,vacc_t2,vacc_t3,vacc_t4,vacc_t5,vacc_t6,vacc_t7,vacc_t8,vacc_t9,vacc_t99)
  }

  # save to parms file
  para$schedule[["vaccination"]] = list(
    parameter = "v",
    pops = numeric(),
    mode = "assign",
    values = vacc_vals_B,
    times = vacc_times_B);
  
  #v = number of people going into compartment S -> V each day (as before).
  #v2 = number of people going from S -> V2 each day.
  #v12 = number of people going from V -> V2 each day.
  #Having these three options is supposed to allow modelling V and V2 as two different vaccines, or as a primary and booster dose of the same vaccine.
  #ev2, ei_v, and ed_vi2 are analogous for vaccine V2.
  
  
  runB = cm_simulate(para, model_seed = seed_val,  n_iter)
  
  # nr of vaccinated individuals (keep for econ analysis later)
  df_vaccpts <- sapply(1:length(para$schedule[["vaccination"]]$values), function(x) {
    
    do.call(rbind,
            rep(para$schedule[["vaccination"]]$values[x], 
                c(c(para$schedule[["vaccination"]]$times, length(unique(runB$dynamics$t)))-
                    lag(c(para$schedule[["vaccination"]]$times, length(unique(runB$dynamics$t)))) )[-1][x] )
    )
  } )
  
  
  # save and return results (for further analysis)
  df_vaccpts <-  as.data.table(do.call(rbind, df_vaccpts))
  
  # replicate by number of iterations
  df_vaccpts <- do.call("rbind", rep(list(df_vaccpts), n_iter))
  
  df_vaccpts$t <- rep( c(0:(max(unique(runB$dynamics$t))) ), n_iter)
  df_vaccpts <- df_vaccpts %>%
    setNames(c(para$pop[[1]]$group_names, "t")) %>%
    tidyr::gather(group, value, -t)
  df_vaccpts$run <- rep(1:n_iter, each=length(unique(runB$dynamics$t)) )
  df_vaccpts$compartment="vaccinees"
  df_vaccpts$population <- para$pop[[1]]$name
  
  # merge
  runB$dynamics <- rbind(runB$dynamics, df_vaccpts)
  
  
  
  ####################################################
  # (C) optimistic vaccination
  ####################################################
  
  if(includeC == TRUE){
    
    n_age_gr_vacc_C   <- cm_age_coefficients(vacc_minAge_C, vacc_maxAge_C, 5 * (0:length(para$pop[[1]]$size)))
    n_age_gr_ttl_C    <- sum(para$pop[[1]]$size[which(cm_age_coefficients(vacc_minAge_C, vacc_maxAge_C, 5 * (0:length(para$pop[[1]]$size)))==1)])
    
    #ev = probability that a person who has received a V vaccine actually moves into V, otherwise they stay in S, so it's like a "take" vaccine effectivemess.
    #ei_v = effectiveness of V vaccine against infection, given that a person is in compartment V.
    #ed_vi = effectiveness of V vaccine against disease, given that a person in compartment V got infected.
    #wv and wv2 are the waning rate of V and V2 back to S.
    para$pop[[1]]$ev    = rep(1, 16) #VE_C * n_age_gr_vacc_C
    para$pop[[1]]$ei_v  = VE_C_infection
    para$pop[[1]]$ed_vi = VE_C_disease
    para$pop[[1]]$wv    = rep((1/wane_vacc_C), n_age_groups) # vaccine protection for x years on average
    
    # numbers of vaccinees
    # uniformly vaccinating; no longer used other than sensitivity analyses
    if(vacc_evenly==TRUE){
      
      vacc_vals_C <- list(rep(0, n_age_groups), # initially no vaccine available
                          # total nr of vaccinees distributed evenly (assumed)
                          c(n_vaccinees_C / sum(n_age_gr_vacc_C) * n_age_gr_vacc_C), # could redistribute more to certain ages here
                          # total nr of vaccinees after having vaccinated all initially (accounting for coverage)
                          c(n_age_gr_ttl_C *cov2_vacc_C /365.25 /sum(n_age_gr_vacc_C) * n_age_gr_vacc_C) ) #re-vaccinate
      
      # timing of initial vaccination and revaccination
      vacc_times_C <- c(0, 
                        # start vaccinating
                        as.numeric(as.Date(start_vacc_C) - as.Date(para$date0) ),
                        # how long initially when vaccinating x a day (all reached, accounting for coverage)
                        ceiling(as.numeric(as.Date(start_vacc_C) - as.Date(para$date0) ) + 
                                  n_age_gr_ttl_C *cov1_vacc_C /n_vaccinees_C) ) # if want to use this by different VE, would need to adjust over/under 50..
    } else {
    # targeted vaccinating in line with provisional vaccination priorities
      
      # moderate to high risk (to do: move this out of function to be customisable; or wait for revision of new compartments)
      #high risk adults under 65 years of age; Clarke et al, females and males in the UK
      n_highRisk  = c(40, 570, 1038, 13140, 27014, 74912, 108766, 144361, 144723, 201420, 270600, 371662, 368724, 413316, 479443, 1178643)
      n_highRisk   = n_highRisk * cm_age_coefficients(0, 65, 5 * (0:length(para$pop[[1]]$size)))
      #moderate/increased risk adults under 65 years of age; Clarke et al, females and males in the UK
      n_modRisk  = c(73604, 118254, 131285, 206125, 366435, 551319, 725170, 827328, 965984, 1177495, 1563243, 1930802, 2011177, 2045102, 2332117, 4746535)
      n_modRisk   = n_modRisk * cm_age_coefficients(0, 65, 5 * (0:length(para$pop[[1]]$size)))
      
      # work out proportions of who received the daily number of vaccines; priority order of JCVI
      #https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/939119/Greenbook_chapter_14a___provisional_guidance_subject_to_MHRA_approval_of_vaccine_supply_.pdf
      vacc_vals_C <-list()
      # initially no vaccine available
      vacc_vals_C[[1]] = rep(0, n_age_groups)
      
      # start with prop in care homes and social care first and be done in 4 weeks, then move on to 75+ (as there are no)
      #291,000 care home residents (260k 75+, 31k 65-74) 
      #1.52M social care  (20-64 y)
      #1.31M million hospital staff (20-64 y)
      n_HCW = 1520000+1310000 # 1.52M social care staff, 1.31M million hospital staff (20-64 y)
      # vaccinate care home pop
      ch_pop = c((31000/2), (31000/2), (260000))
      # vaccinate HCW and CH residents
      pop_v1 = n_HCW / 
        sum(cm_age_coefficients(20, 65, 5*(0:length(para$pop[[1]]$size)))) * 
        cm_age_coefficients(20, 65, 5 * (0:length(para$pop[[1]]$size)))
      pop_v1[14:16] = ch_pop
      # check if aligning with vaccination age groups
      pop_v1 = pop_v1 * cm_age_coefficients(vacc_minAge_C, vacc_maxAge_C, 5 * (0:length(para$pop[[1]]$size)))
      vacc_vals_C[[2]] = round(n_init_vaccinees_C / sum(pop_v1) * pop_v1)
      
      #75+, (once care workers exhausted move to remaining 75+ only)
      vacc_vals_C[[3]] = n_vaccinees_C * cm_age_coefficients(75, 100, 5*(0:length(para$pop[[1]]$size)))
      # check if aligning with vaccination age groups
      vacc_vals_C[[3]] = vacc_vals_C[[3]] * cm_age_coefficients(vacc_minAge_C, vacc_maxAge_C, 5 * (0:length(para$pop[[1]]$size)))
      
      #70+ and high-risk adults under 65 years of age; Clarke et al, females and males in the UK
      # ages 70-74
      n_highRisk[which(cm_age_coefficients(70, 75, 5 * (0:length(para$pop[[1]]$size)))==1)] = 
        para$pop[[1]]$size[which(cm_age_coefficients(70, 75, 5 * (0:length(para$pop[[1]]$size)))==1)]-
        31000/2 # minus healt care residents
      # check if aligning with vaccination age groups
      n_highRisk70 = n_highRisk * cm_age_coefficients(vacc_minAge_C, vacc_maxAge_C, 5 * (0:length(para$pop[[1]]$size)))
      vacc_vals_C[[4]] = round(n_vaccinees_C / sum(n_highRisk70) * n_highRisk70)
      if(any(is.na(vacc_vals_C[[4]] )==TRUE)) vacc_vals_C[[4]] = rep(0, n_age_groups)
      
      #65+
      vacc_vals_C[[5]] = n_vaccinees_C * cm_age_coefficients(65, 70, 5*(0:length(para$pop[[1]]$size)))
      # check if aligning with vaccination age groups
      vacc_vals_C[[5]] = vacc_vals_C[[5]] * cm_age_coefficients(vacc_minAge_C, vacc_maxAge_C, 5 * (0:length(para$pop[[1]]$size)))
      
      #moderate/increased risk adults under 65 years of age; Clarke et al, females and males in the UK
      # check if aligning with vaccination age groups
      n_modRisk = n_modRisk * cm_age_coefficients(vacc_minAge_C, vacc_maxAge_C, 5 * (0:length(para$pop[[1]]$size)))
      vacc_vals_C[[6]] = round(n_vaccinees_C / sum(n_modRisk) * n_modRisk)
      if(any(is.na(vacc_vals_C[[6]] )==TRUE)) vacc_vals_C[[6]] = rep(0, n_age_groups)
      
      #60+,
      vacc_vals_C[[7]] = n_vaccinees_C * cm_age_coefficients(60, 65, 5*(0:length(para$pop[[1]]$size)))
      # check if aligning with vaccination age groups
      vacc_vals_C[[7]] = vacc_vals_C[[7]] * cm_age_coefficients(vacc_minAge_C, vacc_maxAge_C, 5 * (0:length(para$pop[[1]]$size)))
      
      #55+, 
      vacc_vals_C[[8]] = n_vaccinees_C * cm_age_coefficients(55, 60, 5*(0:length(para$pop[[1]]$size)))
      # check if aligning with vaccination age groups
      vacc_vals_C[[8]] = vacc_vals_C[[8]] * cm_age_coefficients(vacc_minAge_C, vacc_maxAge_C, 5 * (0:length(para$pop[[1]]$size)))
      
      #50+, 
      vacc_vals_C[[9]] = n_vaccinees_C * cm_age_coefficients(50, 55, 5*(0:length(para$pop[[1]]$size)))
      # check if aligning with vaccination age groups
      vacc_vals_C[[9]] = vacc_vals_C[[9]] * cm_age_coefficients(vacc_minAge_C, vacc_maxAge_C, 5 * (0:length(para$pop[[1]]$size)))
      
      #rest (vacc_minAge_C-49)
      vacc_vals_C[[10]] = n_vaccinees_C / 
        sum(cm_age_coefficients(vacc_minAge_C, 50, 5*(0:length(para$pop[[1]]$size)))) * 
        cm_age_coefficients(vacc_minAge_C, 50, 5 * (0:length(para$pop[[1]]$size)))
      # check if aligning with vaccination age groups
      vacc_vals_C[[10]] = vacc_vals_C[[10]] * cm_age_coefficients(vacc_minAge_C, vacc_maxAge_C, 5 * (0:length(para$pop[[1]]$size)))
      if(any(is.na(vacc_vals_C[[10]] )==TRUE)) vacc_vals_C[[10]] = rep(0, n_age_groups)
      
      # revaccinate; total nr of vaccinees after having vaccinated all initially (accounting for coverage)
      # only need to revaccinate as many as want to (based on coverage..), or at assumed daily capacity
      vacc_vals_C[[11]] = min(n_revac_vaccinees_C, c(sum(para$pop[[1]]$size *cov2_vacc_C) /365.25)) / 
        sum(n_age_gr_vacc_C) * n_age_gr_vacc_C
      # check if aligning with vaccination age groups
      vacc_vals_C[[11]] = vacc_vals_C[[11]] * cm_age_coefficients(vacc_minAge_C, vacc_maxAge_C, 5 * (0:length(para$pop[[1]]$size)))
      
      
      # start simulation
      vacc_t0 = 0
      
      # vaccination introduction
      vacc_t1 = ceiling(vacc_t0 + as.numeric(as.Date(start_vacc_C) - as.Date(para$date0) ))
      
      # care homes and HCW vaccinated in first ~4 weeks
      vacc_t2 = pop_v1* cm_age_coefficients(vacc_minAge_C, vacc_maxAge_C, 5 * (0:length(para$pop[[1]]$size)))
      vacc_t2 = ceiling(vacc_t1 + sum(vacc_t2*cov1_vacc_O50_C)/n_init_vaccinees_C)
      
      #75+ minus care home residents
      n_vacc = rep(0, n_age_groups)
      # assign 75+
      n_vacc[which(cm_age_coefficients(75, 100, 5 * (0:length(para$pop[[1]]$size)))==1)] = 
        para$pop[[1]]$size[which(cm_age_coefficients(75, 100, 5 * (0:length(para$pop[[1]]$size)))==1)]-260000 # minus healt care residents
      # check if aligning with vaccination age groups
      n_vacc = n_vacc * cm_age_coefficients(vacc_minAge_C, vacc_maxAge_C, 5 * (0:length(para$pop[[1]]$size)))
      # calculate vaccination rate
      vacc_t3 = ceiling(vacc_t2 + sum(n_vacc*cov1_vacc_O50_C)/n_vaccinees_C)
      
      #70+ and high-risk <65 minus care home residents (checked and subtracted within n_highRisk70 already)
      vacc_t4 = ceiling(vacc_t3 + sum(n_highRisk70*cov1_vacc_O50_C) / n_vaccinees_C)
      
      #65+ minus care home residents
      n_vacc = rep(0, n_age_groups)
      n_vacc[which(cm_age_coefficients(65, 70, 5 * (0:length(para$pop[[1]]$size)))==1)] = 
        para$pop[[1]]$size[which(cm_age_coefficients(65, 70, 5 * (0:length(para$pop[[1]]$size)))==1)]-31000/2 # minus care home residents
      # check if aligning with vaccination age groups
      n_vacc = n_vacc * cm_age_coefficients(vacc_minAge_C, vacc_maxAge_C, 5 * (0:length(para$pop[[1]]$size)))
      
      # calculate vaccination rate
      vacc_t5 = ceiling(vacc_t4 + sum(n_vacc*cov1_vacc_O50_C)/n_vaccinees_C)
      
      #moderate risk <65 (checked within n_modRisk already)
      vacc_t6 = ceiling(vacc_t5 + sum(n_modRisk*cov1_vacc_O50_C) / n_vaccinees_C)
      
      #60+ minus moderate-high risk
      n_vacc = rep(0, n_age_groups)
      n_vacc[which(cm_age_coefficients(60, 65, 5 * (0:length(para$pop[[1]]$size)))==1)] = 
        para$pop[[1]]$size[which(cm_age_coefficients(60, 65, 5 * (0:length(para$pop[[1]]$size)))==1)]-
        # minus high- and moderate risk
        n_highRisk[which(cm_age_coefficients(60, 65, 5 * (0:length(para$pop[[1]]$size)))==1)]-
        n_modRisk[which(cm_age_coefficients(60, 65, 5 * (0:length(para$pop[[1]]$size)))==1)]- 
        # minus HCW
        n_HCW/sum(cm_age_coefficients(20, 65, 5 * (0:length(para$pop[[1]]$size)))) 
      # check if aligning with vaccination age groups
      n_vacc = n_vacc * cm_age_coefficients(vacc_minAge_C, vacc_maxAge_C, 5 * (0:length(para$pop[[1]]$size)))
      # calculate vaccination rate
      vacc_t7 = ceiling(vacc_t6 + sum(n_vacc*cov1_vacc_O50_C)/n_vaccinees_C)
      
      #55+ minus moderate-high risk 
      n_vacc = rep(0, n_age_groups)
      n_vacc[which(cm_age_coefficients(55, 60, 5 * (0:length(para$pop[[1]]$size)))==1)] = 
        para$pop[[1]]$size[which(cm_age_coefficients(55, 60, 5 * (0:length(para$pop[[1]]$size)))==1)]-
        # minus high- and moderate risk
        n_highRisk[which(cm_age_coefficients(55, 60, 5 * (0:length(para$pop[[1]]$size)))==1)]-
        n_modRisk[which(cm_age_coefficients(55, 60, 5 * (0:length(para$pop[[1]]$size)))==1)]- 
        # minus HCW
        n_HCW/sum(cm_age_coefficients(20, 65, 5 * (0:length(para$pop[[1]]$size))))
      # check if aligning with vaccination age groups
      n_vacc = n_vacc * cm_age_coefficients(vacc_minAge_C, vacc_maxAge_C, 5 * (0:length(para$pop[[1]]$size)))
      # calculate vaccination rate
      vacc_t8 = ceiling(vacc_t7 + sum(n_vacc*cov1_vacc_O50_C)/n_vaccinees_C)
      
      #50+ minus moderate-high risk 
      n_vacc = rep(0, n_age_groups)
      n_vacc[which(cm_age_coefficients(50, 55, 5 * (0:length(para$pop[[1]]$size)))==1)] = 
        para$pop[[1]]$size[which(cm_age_coefficients(50, 55, 5 * (0:length(para$pop[[1]]$size)))==1)]-
        # minus high- and moderate risk
        n_highRisk[which(cm_age_coefficients(50, 55, 5 * (0:length(para$pop[[1]]$size)))==1)]-
        n_modRisk[which(cm_age_coefficients(50, 55, 5 * (0:length(para$pop[[1]]$size)))==1)]- 
        # minus HCW
        n_HCW/sum(cm_age_coefficients(20, 65, 5 * (0:length(para$pop[[1]]$size))))
      # check if aligning with vaccination age groups
      n_vacc = n_vacc * cm_age_coefficients(vacc_minAge_C, vacc_maxAge_C, 5 * (0:length(para$pop[[1]]$size)))
      # calculate vaccination rate
      vacc_t9 = ceiling(vacc_t8 + sum(n_vacc*cov1_vacc_O50_C)/n_vaccinees_C)
      
      #rest (vacc_minAge_C-49) minus moderate-high risk; only 50% #that's the time after which we start re-vaccinating
      n_vacc = rep(0, n_age_groups)
      n_vacc[which(cm_age_coefficients(vacc_minAge_C, 50, 5 * (0:length(para$pop[[1]]$size)))==1)] = 
        para$pop[[1]]$size[which(cm_age_coefficients(vacc_minAge_C, 50, 5 * (0:length(para$pop[[1]]$size)))==1)]-
        # minus high- and moderate risk
        n_highRisk[which(cm_age_coefficients(vacc_minAge_C, 50, 5 * (0:length(para$pop[[1]]$size)))==1)]-
        n_modRisk[which(cm_age_coefficients(vacc_minAge_C, 50, 5 * (0:length(para$pop[[1]]$size)))==1)]- 
        # minus HCW
        n_HCW/
        sum(cm_age_coefficients(20, 65, 5 * (0:length(para$pop[[1]]$size))))*
        sum(cm_age_coefficients(20, 50, 5 * (0:length(para$pop[[1]]$size))))
      # check if aligning with vaccination age groups
      n_vacc = n_vacc * cm_age_coefficients(vacc_minAge_C, vacc_maxAge_C, 5 * (0:length(para$pop[[1]]$size)))
      
      # implementation of delay in programme by x days
      # implementation either via delay time or via doses
      # drop of all AZ vaccine would be equivalent to roughly 2 months delay
      #sum(n_vacc*cov1_vacc_U50_C)/n_vaccinees_B_delay
      #((sum(n_vacc*cov1_vacc_U50_C)/n_vaccinees_C) +60)
      
      # adjust schedule and timing due to delays
      # original function (if nothing delayed)
      if(is.null(vacc_delay_prop) & is.null(vacc_delay_time)){
        vacc_t99 = ceiling(vacc_t9 + sum(n_vacc*cov1_vacc_U50_C)/n_vaccinees_C )
      }
      
      # if prop drop in vacc doses given
      if(!is.null(vacc_delay_prop) & is.null(vacc_delay_time)){
        #new time
        vacc_t99 = ceiling(vacc_t9 + sum(n_vacc*cov1_vacc_U50_C)/ (n_vaccinees_C*vacc_delay_prop) )
        #new dosing
        vacc_vals_C[[10]][ which(vacc_vals_C[[10]]>0) ] = 
          vacc_vals_C[[10]][ which(vacc_vals_C[[10]]>0) ] /
          vacc_vals_C[[10]][ which(vacc_vals_C[[10]]>0) ] * 
          (n_vaccinees_C*vacc_delay_prop)/length(which(vacc_vals_C[[10]]>0) )
        
      }
      
      # if delay time given
      if(is.null(vacc_delay_prop) & !is.null(vacc_delay_time)){
        
        # new nr of doses
        n_vaccinees_C_delay = 
          n_vaccinees_C*(sum(n_vacc*cov1_vacc_U50_C)/n_vaccinees_C) / 
          ((sum(n_vacc*cov1_vacc_U50_C)/n_vaccinees_C) +vacc_delay_time)
        
        #new time
        vacc_t99 = ceiling(vacc_t9 + sum(n_vacc*cov1_vacc_U50_C)/ n_vaccinees_C_delay )
        #new dosing
        vacc_vals_C[[10]][ which(vacc_vals_C[[10]]>0) ] = 
          vacc_vals_C[[10]][ which(vacc_vals_C[[10]]>0) ] /
          vacc_vals_C[[10]][ which(vacc_vals_C[[10]]>0) ] * 
          (n_vaccinees_C_delay/length(which(vacc_vals_C[[10]]>0) ))
        
      }
      
      # stop if both given.. could optimise min(vacc_delay_prop) and max(vacc_delay_time)
      if(!is.null(vacc_delay_prop) & !is.null(vacc_delay_time)){
        stop("Both delay time and proportion reduction given; not implemented currently")
      }
      
      # get all times of vaccination changes
      vacc_times_C <- c(vacc_t0,vacc_t1,vacc_t2,vacc_t3,vacc_t4,vacc_t5,vacc_t6,vacc_t7,vacc_t8,vacc_t9,vacc_t99)
    }
    
    # save to parms file
    para$schedule[["vaccination"]] = list(
      parameter = "v",
      pops = numeric(),
      mode = "assign",
      values = vacc_vals_C,
      times = vacc_times_C);
    
    runC = cm_simulate(para, model_seed = seed_val, n_iter)
    
    # nr of vaccinated individuals (keep for econ analysis later)
    df_vaccpts <- sapply(1:length(para$schedule[["vaccination"]]$values), function(x) {
      
      do.call(rbind,
              rep(para$schedule[["vaccination"]]$values[x], 
                  c(c(para$schedule[["vaccination"]]$times, length(unique(runC$dynamics$t)))-
                      lag(c(para$schedule[["vaccination"]]$times, length(unique(runC$dynamics$t)))) )[-1][x] )
      )
    } )
    
    
    # save and return results (for further analysis)
    df_vaccpts <-  as.data.table(do.call(rbind, df_vaccpts))
    
    # replicate by number of iterations
    df_vaccpts <- do.call("rbind", rep(list(df_vaccpts), n_iter))
    
    df_vaccpts$t <- rep( c(0:(max(unique(runC$dynamics$t))) ), n_iter)
    df_vaccpts <- df_vaccpts %>%
      setNames(c(para$pop[[1]]$group_names, "t")) %>%
      tidyr::gather(group, value, -t)
    df_vaccpts$run <- rep(1:n_iter, each=length(unique(runC$dynamics$t)) )
    df_vaccpts$compartment <- "vaccinees"
    df_vaccpts$population <- para$pop[[1]]$name
    
    # merge
    runC$dynamics <- rbind(runC$dynamics, df_vaccpts)
  }
  
  
  # save and return results (for further analysis)
  df_res <- data.table::rbindlist(lapply(mget(ls(pattern="^run[A-Z]")), "[[", "dynamics"), idcol = "scenario") 
  
  # get year
  df_res$year <- ifelse(df_res$t == 0, 1, ceiling(df_res$t/365.25))
  
  # keep only compartments needed for epi/econ analysis to make code run more efficiently
  if(storeMinParms==TRUE){
    
    df1 <- df_res[compartment %in% c("Is", "death_o", "nonicu_i", "icu_i", "vaccinees"), .(value = sum(value)), 
                  by = .(scenario, run, group, compartment, year)][order(scenario, run, group, compartment,  year)]
    
    df2 <- df_res[compartment %in% c("icu_p","nonicu_p"), .(value = sum(value)), 
                  by = .(scenario, run, t, group, year)][order(scenario, run, t, group, year)]
    df2$compartment <- "hosp_p"
    
    df3 <- df_res[compartment == "cases", .(value = sum(value)), 
                  by = .(scenario, run, t, 
                         group, 
                         compartment, year)][order(scenario, run, t, 
                                                   group, 
                                                   compartment, year)]
    
    df4 <- df_res[compartment == "obs0" & group %in% c("0-4", "10-14")]# & value %in% 1:3]
    
    df_res <- data.table::rbindlist(list(df1, df2, df3, df4), fill=TRUE)[order(scenario, run, group, compartment,  year)]
  }
  
  return(df_res)
}



##############################
## D) Set-up econ model
##############################

# model set up for baseline A plus 2 scenarios B and C (allowing differences per scenario that are mostly not used currently)
econ_mod <- function(epi_dataset,             # input based on epi_mod
                     n_agegrs        = NULL,  # faster if specified
                     n_years         = NULL,  # faster if specified
                     n_scens         = NULL,  # faster if specified
                     
                     cost_vacc_B     = 30,    # based on 2 doses (2*10); similar to long-run influenza vaccine list prices (BNF; <=£10 even for quadrivalent)
                     cost_vacc_C     = 30,    # based on 2 doses (2*20); similar to Pfizer price for COVID vaccine and accounting for freezers, other opportunity costs of more time potentially needed
                     costs_GPvisit   = TRUE, 
                     costs_calls     = TRUE, 
                     costs_ICU       = TRUE, 
                     costs_ICUbedday = FALSE, 
                     costs_hosp      = TRUE, 
                     costs_vaccAEFI  = TRUE, 
                     costs_vaccAdmin = TRUE, 
                     costs_vacc      = TRUE,  
                     costs_lockdown  = TRUE, 
                     costs_PPE       = TRUE,
                     
                     cost_vaccRD     = 250000000, # £250m lumpsum in total (later divided by age)
                     cost_vaccFreezers = 3300000, # £3.3M supply of ultra-low temperature (ULT) freezers required to store vaccines for the Covid 19 Vaccination Programme; https://bidstats.uk/tenders/2020/W48/739770817
                     
                     c_PD_min          = 115146593,  #2% GDP
                     c_PD_min_trigger  = 1000,
                     c_PD_lockdown     = 287866484,  #5% GDP
                     
                     QALYs_cases    = TRUE, 
                     QALYs_hosp     = TRUE, 
                     QALYs_ICU      = TRUE, 
                     QALYs_vaccAEFI = TRUE, 
                     QALYs_mort     = TRUE, 
                     LEs_mort       = TRUE,  
                     
                     comorb_reduc    = 0.90,      # assume 90% comorbidity reduction
                     discount_rate_c = 0.035,     # discount rate costs
                     discount_rate_b = 0.035,     # discount rate benefits
                     addRD_costs     = TRUE,      # public RD vaccination costs
                     add_Freezercost = TRUE,      # public RD vaccination costs
                     save_by_agegr   = FALSE,     # save by age and year
                     n_iter          = 5,         # number of iterations run the model
                     PSA_var         = 0.25,      # variance used for uncertainty (+/-25%)
                     run_determ      = TRUE       # set stochastic or deterministic
                     ){
  
  # get parameters if not supplied
  if(is.null(n_scens)) { n_scens  <- length(unique(epi_dataset$scenario)) }
  if(is.null(n_agegrs)){ n_agegrs <- length(unique(epi_dataset$group)[!is.na(unique(epi_dataset$group))]) }
  if(is.null(n_years)) { n_years  <- length(unique(epi_dataset$year)[!is.na(unique(epi_dataset$year))]) }
  
  # check order of age groups for faster vector combination
  if(any(unlist(epi_dataset[1:length(params$pop[[1]]$group_names), "group"]) == params$pop[[1]]$group_names)==FALSE) warning("age groups not aligned")
  
  # if deterministic, set n_iter to 1
  if(run_determ==TRUE){
    n_iter <- 1
    params$deterministic = TRUE
  } else {
    params$deterministic = FALSE
    
    # random latin hypercube sampling for PSA later (more efficient form of MC sampling)
    set.seed(0)
    LHS_samples <- lhs::randomLHS(n = n_iter, 
                                  k = sum(costs_GPvisit, costs_calls, costs_ICU, costs_ICUbedday, 
                                          costs_hosp, costs_vaccAEFI, costs_vaccAdmin, costs_vacc,
                                          costs_PPE, costs_lockdown, 
                                          QALYs_cases, QALYs_hosp, QALYs_ICU, QALYs_vaccAEFI, 
                                          QALYs_mort, LEs_mort) )
    n_LHS <- 1
  }
  
  
  # collect results per element
  pr_res <- list()
  
  # dataframe if epi iterations different to econ iterations
  output_res           <- epi_dataset[compartment == "nonicu_i"][order(scenario, run,  year, group, compartment)]
  output_res2          <- data.table(run = rep(1:n_iter, each = length(output_res$run)))
  output_res2$scenario <- rep(output_res$scenario, n_iter)
  output_res2$group    <- rep(output_res$group, n_iter)
  output_res2$year     <- rep(output_res$year, n_iter)
  rm(output_res)
  
  
  ##############################
  ##############################
  # costs NHS perspective (hospitalisation, vaccination, primary care)
  ##############################
  ##############################

  ##############################
  # costs hospitalisation
  ##############################
  
  if(costs_hosp==TRUE){ 
    
    if(run_determ==TRUE){
      
      input_val <- data.frame("group" = names(cost_hosp_m), "mean" = cost_hosp_m)
      
    } else {
      
      # more efficient replication if identical value
      if(length(unique(cost_hosp_m))==TRUE){
        
        input_val <- unique(cost_hosp_m) %>%
          setNames("mean") %>% #need named vector for map_dfc
          purrr::map_dfc(function(x) qlnorm(LHS_samples[,n_LHS], 
                                            meanlog = lognorm_location(x, x*PSA_var), 
                                            sdlog   = lognorm_shape(   x, x*PSA_var)) ) %>%
          mutate(run = 1:n_iter)
        
      } else{
        
        # input data by agegroups
        input_val <- cost_hosp_m %>%
          purrr::map_dfc(function(x) qlnorm(LHS_samples[,n_LHS], 
                                            meanlog= lognorm_location(x, x*PSA_var), 
                                            sdlog  = lognorm_shape(   x, x*PSA_var)) ) %>%
          tidyr::gather(group, mean) %>%
          group_by(group) %>%
          mutate(run = 1:n_iter)
        
      }
      n_LHS <- n_LHS +1
    }

    # subset outcomes needed (re-order in required format; not working well by grouping alone)
    output_val <- epi_dataset[compartment == "nonicu_i"][order(scenario, run,  year, group, compartment)]
    
    # calculate results
    if(identical(length(output_val$value), length(input_val$mean))){
      output_val             <- suppressMessages(full_join(output_val, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val$undisc      <- output_val$value * output_val$mean
      output_val$mean        <- NULL
      output_val$disc        <- output_val$undisc * (1+discount_rate_c)^-(output_val$year-1)
      output_val$compartment <- "costs_hosp"
      
    } else { # different approach if epi results different from econ iterations
      
      output_val2             <- output_res2 # take dataframe by econ iterations
      output_val2$compartment <- "costs_hosp"
      output_val2$undisc      <- rep(output_val$value, n_iter)
      output_val2             <- suppressMessages(full_join(output_val2, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val2$undisc      <- output_val2$undisc * output_val2$mean
      output_val2$mean        <- NULL
      output_val2$disc        <- output_val2$undisc * (1+discount_rate_c)^-c(output_val$year-1)
      output_val              <- output_val2[, lapply(.SD, sum, na.rm=TRUE), by = .(scenario, run, compartment), 
                                             .SDcols=c("undisc", "disc") ]
    }
    
    # save
    if(save_by_agegr==TRUE){
      pr_res[["costs_hosp"]] <- output_val2 
    } else {
      pr_res[["costs_hosp"]] <- output_val
    }
    
  }
  
  
  ##############################
  # costs ICU
  ##############################
  
  if(costs_ICU==TRUE){
    
    if(costs_ICUbedday==TRUE){  
      icu_comp <- "icu_p" # icu_p given bed-day costs
    } else {
      icu_comp <- "icu_i" # icu_i given costs per stay
    }
    
    if(run_determ==TRUE){
      
      input_val <- data.frame("group" = names(cost_icu_m), "mean" = cost_icu_m)
      
    } else {
      
      # more efficient replication if identical value
      if(length(unique(cost_icu_m))==TRUE){
        
        input_val <- unique(cost_icu_m) %>%
          setNames("mean") %>% #need named vector for map_dfc
          purrr::map_dfc(function(x) qlnorm(LHS_samples[,n_LHS], 
                                            meanlog = lognorm_location(x, x*PSA_var), 
                                            sdlog   = lognorm_shape(   x, x*PSA_var)) ) %>%
          mutate(run = 1:n_iter)
        
      } else{
        
        # input data by agegroups
        input_val <- cost_icu_m %>%
          purrr::map_dfc(function(x) qlnorm(LHS_samples[,n_LHS], 
                                            meanlog= lognorm_location(x, x*PSA_var), 
                                            sdlog  = lognorm_shape(   x, x*PSA_var)) ) %>%
          tidyr::gather(group, mean) %>%
          group_by(group) %>%
          mutate(run = 1:n_iter)
        
      }
      n_LHS <- n_LHS +1
    }
    
    # subset outcomes needed (re-order in required format; not working well by grouping alone)
    output_val <- epi_dataset[compartment == icu_comp][order(scenario, run,  year, group, compartment)]
    
    # calculate results if epi iterations equal econ iterations
    if(identical(length(output_val$value), length(input_val$mean))){
      output_val             <- suppressMessages(full_join(output_val, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val$undisc      <- output_val$value * output_val$mean
      output_val$mean        <- NULL
      output_val$disc        <- output_val$undisc * (1+discount_rate_c)^-(output_val$year-1)
      output_val$compartment <- "costs_icu"
      
    } else { # different approach if epi results different from econ iterations
      
      output_val2             <- output_res2 # take dataframe by econ iterations
      output_val2$compartment <- "costs_icu"
      output_val2$undisc      <- rep(output_val$value, n_iter)
      output_val2             <- suppressMessages(full_join(output_val2, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val2$undisc      <- output_val2$undisc * output_val2$mean
      output_val2$mean        <- NULL
      output_val2$disc        <- output_val2$undisc * (1+discount_rate_c)^-c(output_val$year-1)
      output_val              <- output_val2[, lapply(.SD, sum, na.rm=TRUE), by = .(scenario, run, compartment), 
                                             .SDcols=c("undisc", "disc") ]
    }
    
    # save
    if(save_by_agegr==TRUE){
      pr_res[["costs_icu"]] <- output_val2 
    } else {
      pr_res[["costs_icu"]] <- output_val
    }

  }
  
  
  ##############################
  # costs PPE in hospital
  ##############################
  
  if(costs_PPE==TRUE){ 
    
    if(run_determ==TRUE){
      
      input_val <- data.frame("group" = names(cost_hosp_m), "mean" = cost_PPE_m, "run" = 1)
      
    } else {
      
      # more efficient replication if identical value
      if(length(unique(cost_PPE_m))==TRUE){
        
        input_val <- unique(cost_PPE_m) %>%
          setNames("mean") %>% #need named vector for map_dfc
          purrr::map_dfc(function(x) qlnorm(LHS_samples[,n_LHS], 
                                            meanlog = lognorm_location(x, x*PSA_var), 
                                            sdlog   = lognorm_shape(   x, x*PSA_var)) ) %>%
          mutate(run = 1:n_iter)
        
      } else{
        
        # input data by agegroups
        input_val <- cost_PPE_m %>%
          purrr::map_dfc(function(x) qlnorm(LHS_samples[,n_LHS], 
                                            meanlog= lognorm_location(x, x*PSA_var), 
                                            sdlog  = lognorm_shape(   x, x*PSA_var)) ) %>%
          tidyr::gather(group, mean) %>%
          group_by(group) %>%
          mutate(run = 1:n_iter)
        
      }
      n_LHS <- n_LHS +1
    }
    
    # subset outcomes needed (re-order in required format; not working well by grouping alone)
    output_val <- epi_dataset[compartment == "hosp_p"][order(scenario, run, t, group, year)]
    
    # calculate results if epi iterations equal econ iterations
    if(identical(length(output_val$value), length(input_val$mean))){
      output_val             <- output_val[, .(value = sum(value)), by = .(scenario, run, t, year)][order(scenario, run, t, year)]
      output_val             <- suppressMessages(full_join(output_val, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val$undisc2     <- output_val$value * output_val$mean
      output_val$mean        <- NULL
      output_val$undisc      <- ifelse(output_val$value<1, 0, output_val$undisc2)
      output_val$disc        <- output_val$undisc * (1+discount_rate_c)^-(output_val$year-1)
      output_val$compartment <- "costs_PPE"
      output_val$undisc2     <- output_val$value <- NULL
      
    } else { # different approach if epi results different from econ iterations
      
      output_val2             <- data.table(run = rep(1:n_iter, each = length(output_val$run)))
      output_val2$scenario    <- rep(output_val$scenario, n_iter)
      output_val2$t           <- rep(output_val$t, n_iter)
      output_val2$year        <- rep(output_val$year, n_iter)
      output_val2$value       <- rep(output_val$value, n_iter)
      output_val2             <- output_val2[, .(value = sum(value)), by = .(scenario, run, t, year)][order(scenario, run, t, year)]
      output_val2$compartment <- "costs_PPE"
      output_val2             <- suppressMessages(full_join(output_val2, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val2$undisc2     <- output_val2$value * output_val2$mean
      output_val2$mean        <- NULL
      output_val2$undisc      <- ifelse(output_val2$value<1, 0, output_val2$undisc2)
      output_val2$disc        <- output_val2$undisc * (1+discount_rate_c)^-c(output_val2$year-1)
      output_val              <- output_val2[, lapply(.SD, sum, na.rm=TRUE), by = .(scenario, run, compartment), 
                                             .SDcols=c("undisc", "disc") ]
    }
    
    # save
    if(save_by_agegr==TRUE){
      pr_res[["costs_PPE"]] <- output_val2 
    } else {
      pr_res[["costs_PPE"]] <- output_val
    }
    
  }
  
  
  ##############################
  # costs GP visits
  ##############################
  
  if(costs_GPvisit==TRUE){
    
    # replicate values for each year and all scenarios
    if(run_determ==TRUE){
      
      input_val <- data.frame("group" = names(cost_GP_m), "mean" = cost_GP_m)
      
    } else {
      
      # more efficient replication if identical input value
      if(length(unique(cost_GP_m))==TRUE){
        
        input_val <- unique(cost_GP_m) %>%
          setNames("mean") %>% #need named vector for map_dfc
          purrr::map_dfc(function(x) qlnorm(LHS_samples[,n_LHS], 
                                            meanlog = lognorm_location(x, x*PSA_var), 
                                            sdlog   = lognorm_shape(   x, x*PSA_var)) ) %>%
          mutate(run = 1:n_iter)
        
      } else{
        
        # input data by agegroups
        input_val <- cost_GP_m %>%
          purrr::map_dfc(function(x) qlnorm(LHS_samples[,n_LHS], 
                                            meanlog= lognorm_location(x, x*PSA_var), 
                                            sdlog  = lognorm_shape(   x, x*PSA_var)) ) %>%
          tidyr::gather(group, mean) %>%
          group_by(group) %>%
          mutate(run = 1:n_iter)
      }
      
      n_LHS <- n_LHS +1
    }
    
    # subset outcomes needed (re-order in required format; not working well by grouping alone)
    output_val <- epi_dataset[compartment == "Is"][order(scenario, run,  year, group, compartment)]
    
    # calculate results if epi iterations equal econ iterations
    if(identical(length(output_val$value), length(input_val$mean))){
      output_val             <- suppressMessages(full_join(output_val, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val$undisc      <- output_val$value * output_val$mean
      output_val$mean        <- NULL
      output_val$disc        <- output_val$undisc * (1+discount_rate_c)^-(output_val$year-1)
      output_val$compartment <- "costs_GPvisit"
      
    } else { # different approach if epi iterations different from econ iterations
      
      output_val2             <- output_res2 # take dataframe by econ iterations
      output_val2$compartment <- "costs_GPvisit"
      output_val2$undisc      <- rep(output_val$value, n_iter)
      output_val2             <- suppressMessages(full_join(output_val2, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val2$undisc      <- output_val2$undisc * output_val2$mean
      output_val2$mean        <- NULL
      output_val2$disc        <- output_val2$undisc * (1+discount_rate_c)^-c(output_val$year-1)
      output_val              <- output_val2[, lapply(.SD, sum, na.rm=TRUE), by = .(scenario, run, compartment), 
                                             .SDcols=c("undisc", "disc") ]
    }
    
    # save
    if(save_by_agegr==TRUE){
      pr_res[["costs_GPvisit"]] <- output_val2 
    } else {
      pr_res[["costs_GPvisit"]] <- output_val
    }
  }
  
  
  ##############################
  # costs phone calls
  ##############################
  
  if(costs_calls==TRUE){
    
    if(run_determ==TRUE){
      
      input_val <- data.frame("group" = names(cost_calls_m), "mean" = cost_calls_m)
      
    } else {
      
      # more efficient replication if identical value
      if(length(unique(cost_calls_m))==TRUE){
        
        input_val <- unique(cost_calls_m) %>%
          setNames("mean") %>% #need named vector for map_dfc
          purrr::map_dfc(function(x) qlnorm(LHS_samples[,n_LHS], 
                                            meanlog = lognorm_location(x, x*PSA_var), 
                                            sdlog   = lognorm_shape(   x, x*PSA_var)) ) %>%
          mutate(run = 1:n_iter)
        
      } else{
        
        # input data by agegroups
        input_val <- cost_calls_m %>%
          purrr::map_dfc(function(x) qlnorm(LHS_samples[,n_LHS], 
                                            meanlog= lognorm_location(x, x*PSA_var), 
                                            sdlog  = lognorm_shape(   x, x*PSA_var)) ) %>%
          tidyr::gather(group, mean) %>%
          group_by(group) %>%
          mutate(run = 1:n_iter)
        
      }
      n_LHS <- n_LHS +1
    }
    
    # subset outcomes needed (re-order in required format; not working well by grouping alone)
    output_val <- epi_dataset[compartment == "Is"][order(scenario, run,  year, group, compartment)]
    
    # calculate results if epi iterations equal econ iterations
    if(identical(length(output_val$value), length(input_val$mean))){
      output_val             <- suppressMessages(full_join(output_val, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val$undisc      <- output_val$value * output_val$mean
      output_val$mean        <- NULL
      output_val$disc        <- output_val$undisc * (1+discount_rate_c)^-(output_val$year-1)
      output_val$compartment <- "costs_calls"
      
    } else { # different approach if epi results different from econ iterations
      
      output_val2             <- output_res2 # take dataframe by econ iterations
      output_val2$compartment <- "costs_calls"
      output_val2$undisc      <- rep(output_val$value, n_iter)
      output_val2             <- suppressMessages(full_join(output_val2, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val2$undisc      <- output_val2$undisc * output_val2$mean
      output_val2$mean        <- NULL
      output_val2$disc        <- output_val2$undisc * (1+discount_rate_c)^-c(output_val$year-1)
      output_val              <- output_val2[, lapply(.SD, sum, na.rm=TRUE), by = .(scenario, run, compartment), 
                                             .SDcols=c("undisc", "disc") ]
    }
    
    # save
    if(save_by_agegr==TRUE){
      pr_res[["costs_calls"]] <- output_val2 
    } else {
      pr_res[["costs_calls"]] <- output_val
    }
    
  }
  

  ##############################
  ## costs adverse events following immunisation (AEFI)
  ##############################
  
  if(costs_vaccAEFI==TRUE){ 
    
    if(run_determ==TRUE){
      
      input_val <- data.frame("group" = names(cost_vaccAEFI_m), "mean" = cost_vaccAEFI_m)
      
    } else {
      
      # more efficient replication if identical value
      if(length(unique(cost_vaccAEFI_m))==TRUE){
        
        input_val <- unique(cost_vaccAEFI_m) %>%
          setNames("mean") %>% #need named vector for map_dfc
          purrr::map_dfc(function(x) qlnorm(LHS_samples[,n_LHS], 
                                            meanlog = lognorm_location(x, x*PSA_var), 
                                            sdlog   = lognorm_shape(   x, x*PSA_var)) ) %>%
          mutate(run = 1:n_iter)
        
      } else{
        
        # input data by agegroups
        input_val <- cost_vaccAEFI_m %>%
          purrr::map_dfc(function(x) qlnorm(LHS_samples[,n_LHS], 
                                            meanlog= lognorm_location(x, x*PSA_var), 
                                            sdlog  = lognorm_shape(   x, x*PSA_var)) ) %>%
          tidyr::gather(group, mean) %>%
          group_by(group) %>%
          mutate(run = 1:n_iter)
        
      }
      n_LHS <- n_LHS +1
    }
    
    # subset outcomes needed (re-order in required format; not working well by grouping alone)
    output_val <- epi_dataset[compartment == "vaccinees"][order(scenario, run,  year, group, compartment)]
    
    # calculate results if epi iterations equal econ iterations
    if(identical(length(output_val$value), length(input_val$mean))){
      output_val             <- suppressMessages(full_join(output_val, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val$undisc      <- output_val$value * output_val$mean
      output_val$mean        <- NULL
      output_val$disc        <- output_val$undisc * (1+discount_rate_c)^-(output_val$year-1)
      output_val$compartment <- "costs_vaccAEFI"
      
    } else { # different approach if epi results different from econ iterations
      
      output_val2             <- output_res2 # take dataframe by econ iterations
      output_val2$compartment <- "costs_vaccAEFI"
      output_val2$undisc      <- rep(output_val$value, n_iter)
      output_val2             <- suppressMessages(full_join(output_val2, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val2$undisc      <- output_val2$undisc * output_val2$mean
      output_val2$mean        <- NULL
      output_val2$disc        <- output_val2$undisc * (1+discount_rate_c)^-c(output_val$year-1)
      output_val              <- output_val2[, lapply(.SD, sum, na.rm=TRUE), by = .(scenario, run, compartment), 
                                             .SDcols=c("undisc", "disc") ]
    }
    
    # save
    if(save_by_agegr==TRUE){
      pr_res[["costs_vaccAEFI"]] <- output_val2 
    } else {
      pr_res[["costs_vaccAEFI"]] <- output_val
    }
  }
  
  
  ##############################
  ## costs of vaccine administration
  ##############################
  
  if(costs_vaccAdmin==TRUE){
    
    if(run_determ==TRUE){
      
      input_val <- data.frame("group" = names(cost_vaccAdmin_m), "mean" = cost_vaccAdmin_m)
      
    } else {
      
      # more efficient replication if identical value
      if(length(unique(cost_vaccAdmin_m))==TRUE){
        
        input_val <- unique(cost_vaccAdmin_m) %>%
          setNames("mean") %>% #need named vector for map_dfc
          purrr::map_dfc(function(x) qlnorm(LHS_samples[,n_LHS], 
                                            meanlog = lognorm_location(x, x*PSA_var), 
                                            sdlog   = lognorm_shape(   x, x*PSA_var)) ) %>%
          mutate(run = 1:n_iter)
        
      } else{
        
        # input data by agegroups
        input_val <- cost_vaccAdmin_m %>%
          purrr::map_dfc(function(x) qlnorm(LHS_samples[,n_LHS], 
                                            meanlog= lognorm_location(x, x*PSA_var), 
                                            sdlog  = lognorm_shape(   x, x*PSA_var)) ) %>%
          tidyr::gather(group, mean) %>%
          group_by(group) %>%
          mutate(run = 1:n_iter)
        
      }
      n_LHS <- n_LHS +1
    }
    
    # subset outcomes needed (re-order in required format; not working well by grouping alone)
    output_val <- epi_dataset[compartment == "vaccinees"][order(scenario, run,  year, group, compartment)]
    
    # calculate results if epi iterations equal econ iterations
    if(identical(length(output_val$value), length(input_val$mean))){
      output_val             <- suppressMessages(full_join(output_val, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val$undisc      <- output_val$value * output_val$mean
      output_val$mean        <- NULL
      output_val$disc        <- output_val$undisc * (1+discount_rate_c)^-(output_val$year-1)
      output_val$compartment <- "costs_vaccAdmin"
      
    } else { # different approach if epi results different from econ iterations
      
      output_val2             <- output_res2 # take dataframe by econ iterations
      output_val2$compartment <- "costs_vaccAdmin"
      output_val2$undisc      <- rep(output_val$value, n_iter)
      output_val2             <- suppressMessages(full_join(output_val2, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val2$undisc      <- output_val2$undisc * output_val2$mean
      output_val2$mean        <- NULL
      output_val2$disc        <- output_val2$undisc * (1+discount_rate_c)^-c(output_val$year-1)
      output_val              <- output_val2[, lapply(.SD, sum, na.rm=TRUE), by = .(scenario, run, compartment), 
                                             .SDcols=c("undisc", "disc") ]
    }
    
    # save
    if(save_by_agegr==TRUE){
      pr_res[["costs_vaccAdmin"]] <- output_val2 
    } else {
      pr_res[["costs_vaccAdmin"]] <- output_val
    }
  }
  
  
  ##############################
  ## costs/price of vaccination, per vaccinee (non-probabilistic)
  ##############################
  
  if(costs_vacc==TRUE){ 
    
    # calculate results (non-probabilistic prices)
    input_val <- data.frame("group" = names(cost_vaccAdmin_m), 
                            "runA"  = 0,
                            "runB"  = cost_vacc_B,
                            "runC"  = cost_vacc_C) %>%
      tidyr::gather(scenario, mean, -group)

    # subset outcomes needed (re-order in required format; not working well by grouping alone)
    output_val <- epi_dataset[compartment == "vaccinees"][order(scenario, run,  year, group, compartment)]
    
    # calculate results if epi iterations equal econ iterations
    if(identical(length(output_val$value), length(input_val$mean))){
      output_val             <- suppressMessages(full_join(output_val, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val$undisc      <- output_val$value * output_val$mean
      output_val$mean        <- NULL
      output_val$disc        <- output_val$undisc * (1+discount_rate_c)^-(output_val$year-1)
      output_val$compartment <- "costs_vacc"
      
    } else { # different approach if epi results different from econ iterations
      
      output_val2             <- output_res2 # take dataframe by econ iterations
      output_val2$compartment <- "costs_vacc"
      output_val2$undisc      <- rep(output_val$value, n_iter)
      output_val2             <- suppressMessages(full_join(output_val2, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val2$undisc      <- output_val2$undisc * output_val2$mean
      output_val2$mean        <- NULL
      output_val2$disc        <- output_val2$undisc * (1+discount_rate_c)^-c(output_val$year-1)
      output_val              <- output_val2[, lapply(.SD, sum, na.rm=TRUE), by = .(scenario, run, compartment), 
                                             .SDcols=c("undisc", "disc") ]
    }
    
    # save
    if(save_by_agegr==TRUE){
      pr_res[["costs_vacc"]] <- output_val2 
    } else {
      pr_res[["costs_vacc"]] <- output_val
    }
  }
  
  
  ##############################
  ## GDP impact lockdown
  ##############################
  
  if(costs_lockdown==TRUE){ 
    
    # replicate values for each year and all scenarios
    if(run_determ==TRUE){
      
      input_val <- data.frame("group" = names(cost_GP_m), "mean" = 1)
      
    } else {
      
      input_val <- 1 %>%
        setNames("mean") %>% #need named vector for map_dfc
        purrr::map_dfc(function(x) qlnorm(LHS_samples[,n_LHS], 
                                          meanlog = lognorm_location(x, x*PSA_var), 
                                          sdlog   = lognorm_shape(   x, x*PSA_var)) ) %>%
        mutate(run = 1:n_iter)
      
      n_LHS <- n_LHS +1
    }
    
    # subset outcomes needed (re-order in required format; not working well by grouping alone)
    output_val <- epi_dataset[compartment == "cases"][order(scenario, run, t, year, compartment)]
    
    # subset periods of PD/lockdown
    output_val2 <- epi_dataset[compartment == "obs0" & group %in% c("0-4") & value %in% 1:3]
    output_val2$lockdown <- output_val2$value
    
    # merge
    output_val <- full_join(output_val, 
                            output_val2[,c("scenario", "run", "t","year", "lockdown")], 
                            by = c("scenario", "run", "t", "year"))
    
    # calculate results; indicate whether low PD or lockdown
    # NEW for AEFI: get cases per day to enact economic lockdown losses
    output_val = output_val %>%
      group_by(scenario, run, compartment, year, t, population) %>%
      mutate(cases_threshold = sum(value, na.rm=TRUE))
    
    output_val$undisc <- ifelse(is.na(output_val$lockdown) & output_val$cases_threshold > c_PD_min_trigger, 1, #c_PD_min,
                                ifelse(!is.na(output_val$lockdown), 2, #c_PD_lockdown,
                                       0))

    # calculate results if epi iterations equal econ iterations
    if(identical(length(output_val$value), length(input_val$mean))){
      output_val             <- suppressMessages(full_join(output_val, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val$undisc      <- ifelse(output_val$undisc == 1, output_val$mean*c_PD_min,
                                ifelse(output_val$undisc == 2, output_val$mean*c_PD_lockdown,
                                       0))
      #output_val$undisc      <- output_val$value * output_val$mean
      output_val$mean        <- NULL
      output_val$disc        <- output_val$undisc * (1+discount_rate_c)^-(output_val$year-1)
      output_val$compartment <- "costs_lockdown"
      output_val             <- output_val[, lapply(.SD, sum, na.rm=TRUE), by = .(scenario, run, compartment), 
                                           .SDcols=c("undisc", "disc") ]
      
    } else { # different approach if epi iterations different from econ iterations
    
      # calculate results
      output_val2             <- data.table(run = rep(1:n_iter, each = length(output_val$run)))
      output_val2$scenario    <- rep(output_val$scenario, n_iter)
      output_val2$group       <- rep(output_val$group, n_iter)
      
      output_val2$t           <- rep(output_val$t, n_iter)
      output_val2$year        <- rep(output_val$year, n_iter)
      
      output_val2$compartment <- "costs_lockdown"
      output_val2$undisc      <- rep(output_val$undisc, n_iter)
      output_val2             <- suppressMessages(full_join(output_val2, input_val)) # suppress as either join by run, or run and AgeGroup
      
      output_val2$undisc        <- ifelse(output_val2$undisc == 1, output_val2$mean*c_PD_min,
                                   ifelse(output_val2$undisc == 2, output_val2$mean*c_PD_lockdown,
                                          0))
      #output_val2$undisc      <- output_val2$undisc * output_val2$mean
      output_val2$mean        <- NULL
      output_val2$disc        <- output_val2$undisc * (1+discount_rate_c)^-c(output_val$year-1)
      output_val              <- output_val2[, lapply(.SD, sum, na.rm=TRUE), by = .(scenario, run, compartment), 
                                             .SDcols=c("undisc", "disc") ]
    }
    
    # save
    if(save_by_agegr==TRUE){
      pr_res[["costs_lockdown"]] <- output_val2 
    } else {
      pr_res[["costs_lockdown"]] <- output_val
    }
  }
  
  
  
  ##############################
  ##############################
  # QALYs lost (NHS perspective)
  ##############################
  ##############################
  
  ##############################
  # QALYs lost per symptomatic episode
  ##############################
  
  if(QALYs_cases==TRUE){
    
    if(run_determ==TRUE){
      
      input_val <- data.frame("group" = names(QALY_ILI_m), "mean" = QALY_ILI_m)
      
    } else {
      
      # more efficient replication if identical value
      if(length(unique(QALY_ILI_m))==TRUE){
        
        input_val <- unique(QALY_ILI_m) %>%
          setNames("mean") %>% #need named vector for map_dfc
          purrr::map_dfc(function(x) qbeta(LHS_samples[,n_LHS], 
                                           shape1 = beta_dis_alpha(x, x*PSA_var), 
                                           shape2 = beta_dis_beta(   x, x*PSA_var)) ) %>%
          mutate(run = 1:n_iter)
        
      } else{
        
        # input data by agegroups
        input_val <- QALY_ILI_m %>%
          purrr::map_dfc(function(x) qbeta(LHS_samples[,n_LHS], 
                                           shape1=beta_dis_alpha(x, x*PSA_var), 
                                           shape2=beta_dis_beta(   x, x*PSA_var)) ) %>%
          tidyr::gather(group, mean) %>%
          group_by(group) %>%
          mutate(run = 1:n_iter)
        
      }
      n_LHS <- n_LHS +1
    }
    
    # subset outcomes needed (re-order in required format; not working well by grouping alone)
    output_val <- epi_dataset[compartment == "Is"][order(scenario, run,  year, group, compartment)]
    
    # calculate results if epi iterations equal econ iterations
    if(identical(length(output_val$value), length(input_val$mean))){
      output_val             <- suppressMessages(full_join(output_val, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val$undisc      <- output_val$value * output_val$mean
      output_val$mean        <- NULL
      output_val$disc        <- output_val$undisc * (1+discount_rate_b)^-(output_val$year-1)
      output_val$compartment <- "QALYs_cases"
      
    } else { # different approach if epi results different from econ iterations
      
      output_val2             <- output_res2 # take dataframe by econ iterations
      output_val2$compartment <- "QALYs_cases"
      output_val2$undisc      <- rep(output_val$value, n_iter)
      output_val2             <- suppressMessages(full_join(output_val2, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val2$undisc      <- output_val2$undisc * output_val2$mean
      output_val2$mean        <- NULL
      output_val2$disc        <- output_val2$undisc * (1+discount_rate_b)^-c(output_val$year-1)
      output_val              <- output_val2[, lapply(.SD, sum, na.rm=TRUE), by = .(scenario, run, compartment), 
                                             .SDcols=c("undisc", "disc") ]
    }
    
    # save
    if(save_by_agegr==TRUE){
      pr_res[["QALYs_cases"]] <- output_val2 
    } else {
      pr_res[["QALYs_cases"]] <- output_val
    }
  }
  
  ##############################
  # QALYs lost per hospitalisation
  ##############################
  
  if(QALYs_hosp==TRUE){ 
    
    if(run_determ==TRUE){
      
      input_val <- data.frame("group" = names(QALY_hosp_m), "mean" = QALY_hosp_m)
      
    } else {
      
      # more efficient replication if identical value
      if(length(unique(QALY_hosp_m))==TRUE){
        
        input_val <- unique(QALY_hosp_m) %>%
          setNames("mean") %>% #need named vector for map_dfc
          purrr::map_dfc(function(x) qbeta(LHS_samples[,n_LHS], 
                                           shape1 = beta_dis_alpha(x, x*PSA_var), 
                                           shape2 = beta_dis_beta(   x, x*PSA_var)) ) %>%
          mutate(run = 1:n_iter)
        
      } else{
        
        # input data by agegroups
        input_val <- QALY_hosp_m %>%
          purrr::map_dfc(function(x) qbeta(LHS_samples[,n_LHS], 
                                           shape1=beta_dis_alpha(x, x*PSA_var), 
                                           shape2=beta_dis_beta(   x, x*PSA_var)) ) %>%
          tidyr::gather(group, mean) %>%
          group_by(group) %>%
          mutate(run = 1:n_iter)
        
      }
      n_LHS <- n_LHS +1
    }
    
    # subset outcomes needed (re-order in required format; not working well by grouping alone)
    output_val <- epi_dataset[compartment == "nonicu_i"][order(scenario, run,  year, group, compartment)]
    
    # calculate results if epi iterations equal econ iterations
    if(identical(length(output_val$value), length(input_val$mean))){
      output_val             <- suppressMessages(full_join(output_val, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val$undisc      <- output_val$value * output_val$mean
      output_val$mean        <- NULL
      output_val$disc        <- output_val$undisc * (1+discount_rate_b)^-(output_val$year-1)
      output_val$compartment <- "QALYs_hosp"
      
    } else { # different approach if epi results different from econ iterations
      
      output_val2             <- output_res2 # take dataframe by econ iterations
      output_val2$compartment <- "QALYs_hosp"
      output_val2$undisc      <- rep(output_val$value, n_iter)
      output_val2             <- suppressMessages(full_join(output_val2, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val2$undisc      <- output_val2$undisc * output_val2$mean
      output_val2$mean        <- NULL
      output_val2$disc        <- output_val2$undisc * (1+discount_rate_b)^-c(output_val$year-1)
      output_val              <- output_val2[, lapply(.SD, sum, na.rm=TRUE), by = .(scenario, run, compartment), 
                                             .SDcols=c("undisc", "disc") ]
    }
    
    # save
    if(save_by_agegr==TRUE){
      pr_res[["QALYs_hosp"]] <- output_val2 
    } else {
      pr_res[["QALYs_hosp"]] <- output_val
    }
  }
  
  ##############################
  # QALYs lost per ICU bed-day/stay
  ##############################
  
  if(QALYs_ICU==TRUE){
    
    if(run_determ==TRUE){
      
      input_val <- data.frame("group" = names(QALY_ICU_m), "mean" = QALY_ICU_m)
      
    } else {
      
      # more efficient replication if identical value
      if(length(unique(QALY_ICU_m))==TRUE){
        
        input_val <- unique(QALY_ICU_m) %>%
          setNames("mean") %>% #need named vector for map_dfc
          purrr::map_dfc(function(x) qbeta(LHS_samples[,n_LHS], 
                                           shape1 = beta_dis_alpha(x, x*PSA_var), 
                                           shape2 = beta_dis_beta(   x, x*PSA_var)) ) %>%
          mutate(run = 1:n_iter)
        
      } else{
        
        # input data by agegroups
        input_val <- QALY_ICU_m %>%
          purrr::map_dfc(function(x) qbeta(LHS_samples[,n_LHS], 
                                           shape1=beta_dis_alpha(x, x*PSA_var), 
                                           shape2=beta_dis_beta(   x, x*PSA_var)) ) %>%
          tidyr::gather(group, mean) %>%
          group_by(group) %>%
          mutate(run = 1:n_iter)
        
      }
      n_LHS <- n_LHS +1
    }
    
    # subset outcomes needed (re-order in required format; not working well by grouping alone)
    output_val <- epi_dataset[compartment == "icu_i"][order(scenario, run,  year, group, compartment)]
    
    # calculate results if epi iterations equal econ iterations
    if(identical(length(output_val$value), length(input_val$mean))){
      output_val             <- suppressMessages(full_join(output_val, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val$undisc      <- output_val$value * output_val$mean
      output_val$mean        <- NULL
      output_val$disc        <- output_val$undisc * (1+discount_rate_b)^-(output_val$year-1)
      output_val$compartment <- "QALYs_ICU"
      
    } else { # different approach if epi results different from econ iterations
      
      output_val2             <- output_res2 # take dataframe by econ iterations
      output_val2$compartment <- "QALYs_ICU"
      output_val2$undisc      <- rep(output_val$value, n_iter)
      output_val2             <- suppressMessages(full_join(output_val2, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val2$undisc      <- output_val2$undisc * output_val2$mean
      output_val2$mean        <- NULL
      output_val2$disc        <- output_val2$undisc * (1+discount_rate_b)^-c(output_val$year-1)
      output_val              <- output_val2[, lapply(.SD, sum, na.rm=TRUE), by = .(scenario, run, compartment), 
                                             .SDcols=c("undisc", "disc") ]
    }
    
    # save
    if(save_by_agegr==TRUE){
      pr_res[["QALYs_ICU"]] <- output_val2 
    } else {
      pr_res[["QALYs_ICU"]] <- output_val
    }
  }
  
  
  ##############################
  # QALYs lost per vaccination (adverse events)
  ##############################
  
  if(QALYs_vaccAEFI==TRUE){
    
    if(run_determ==TRUE){
      
      input_val <- data.frame("group" = names(QALY_vaccAEFI_m), "mean" = QALY_vaccAEFI_m)
      
    } else {
      
      # more efficient replication if identical value
      if(length(unique(QALY_vaccAEFI_m))==TRUE){
        
        input_val <- unique(QALY_vaccAEFI_m) %>%
          setNames("mean") %>% #need named vector for map_dfc
          purrr::map_dfc(function(x) qbeta(LHS_samples[,n_LHS], 
                                           shape1 = beta_dis_alpha(x, x*PSA_var), 
                                           shape2 = beta_dis_beta(   x, x*PSA_var)) ) %>%
          mutate(run = 1:n_iter)
        
      } else{
        
        # input data by agegroups
        input_val <- QALY_vaccAEFI_m %>%
          purrr::map_dfc(function(x) qbeta(LHS_samples[,n_LHS], 
                                           shape1=beta_dis_alpha(x, x*PSA_var), 
                                           shape2=beta_dis_beta(   x, x*PSA_var)) ) %>%
          tidyr::gather(group, mean) %>%
          group_by(group) %>%
          mutate(run = 1:n_iter)
        
      }
      n_LHS <- n_LHS +1
    }
    
    # subset outcomes needed (re-order in required format; not working well by grouping alone)
    output_val <- epi_dataset[compartment == "vaccinees"][order(scenario, run,  year, group, compartment)]
    
    # calculate results if epi iterations equal econ iterations
    if(identical(length(output_val$value), length(input_val$mean))){
      output_val             <- suppressMessages(full_join(output_val, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val$undisc      <- output_val$value * output_val$mean
      output_val$mean        <- NULL
      output_val$disc        <- output_val$undisc * (1+discount_rate_b)^-(output_val$year-1)
      output_val$compartment <- "QALYs_vaccAEFI"
      
    } else { # different approach if epi results different from econ iterations
      
      output_val2             <- output_res2 # take dataframe by econ iterations
      output_val2$compartment <- "QALYs_vaccAEFI"
      output_val2$undisc      <- rep(output_val$value, n_iter)
      output_val2             <- suppressMessages(full_join(output_val2, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val2$undisc      <- output_val2$undisc * output_val2$mean
      output_val2$mean        <- NULL
      output_val2$disc        <- output_val2$undisc * (1+discount_rate_b)^-c(output_val$year-1)
      output_val              <- output_val2[, lapply(.SD, sum, na.rm=TRUE), by = .(scenario, run, compartment), 
                                             .SDcols=c("undisc", "disc") ]
    }
    
    # save
    if(save_by_agegr==TRUE){
      pr_res[["QALYs_vaccAEFI"]] <- output_val2 
    } else {
      pr_res[["QALYs_vaccAEFI"]] <- output_val
    }
  }
  
  
  ##############################
  # QALYs lost per premature death
  ##############################
  
  if(QALYs_mort==TRUE){
    
    # discount quality-adjusted life-expectancy (QALE); different utilities by sex
    for(j in c("female","male")){
      if(j == "female") {fem=TRUE} else {fem=FALSE}
      
      # discounting
      LE_UK[[paste0(j, "_QALE_disc", sep="")]] <- 
        LE_UK[, j] %>% unlist() %>% 
        purrr::map_dbl(function(x) discount_LY(x, x, comorb_reduc, discount_rate_b, female=fem) )
    }
    
    # mean QALYs total by age
    LE_UK <- LE_UK %>%
      transmute(group  = age,
                undisc = (female                 +male)    /2, 
                disc   = (female_QALE_disc +male_QALE_disc)/2 )

    
    if(run_determ==TRUE){
      
      input_val <- data.frame("group" = LE_UK$group, "mean" = LE_UK$disc)

    } else {

      # input data by agegroups
      input_val <- LE_UK$disc %>%
        setNames(LE_UK$group) %>%
        purrr::map_dfc(function(x) qnorm(LHS_samples[, n_LHS], x*1, x*0.1) ) %>%
        tidyr::gather(group, mean) %>%
        group_by(group) %>%
        mutate(run = 1:n_iter)

      n_LHS <- n_LHS +1
    }
    
    # subset outcomes needed
    output_val <- epi_dataset[compartment == "death_o"][order(scenario, run,  year, group, compartment)]
    
    # calculate results if epi iterations equal econ iterations
    if(identical(length(output_val$value), length(input_val$mean))){
      output_val             <- suppressMessages(full_join(output_val, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val$undisc      <- output_val$value * output_val$mean
      output_val$mean        <- NULL
      output_val$disc        <- output_val$undisc * (1+discount_rate_b)^-(output_val$year-1)
      output_val$compartment <- "QALYs_mort"
      
    } else { # different approach if epi results different from econ iterations
      
      output_val2             <- output_res2 # take dataframe by econ iterations
      output_val2$compartment <- "QALYs_mort"
      output_val2$undisc      <- rep(output_val$value, n_iter)
      output_val2             <- suppressMessages(full_join(output_val2, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val2$undisc      <- output_val2$undisc * output_val2$mean
      output_val2$mean        <- NULL
      output_val2$disc        <- output_val2$undisc * (1+discount_rate_b)^-c(output_val$year-1)
      output_val              <- output_val2[, lapply(.SD, sum, na.rm=TRUE), by = .(scenario, run, compartment), 
                                             .SDcols=c("undisc", "disc") ]
    }
    
    # save
    if(save_by_agegr==TRUE){
      pr_res[["QALYs_mort"]] <- output_val2 
    } else {
      pr_res[["QALYs_mort"]] <- output_val
    }
  }
  
  
  ##############################
  # LEs lost per premature death
  ##############################
  
  if(LEs_mort==TRUE){
    
    # discount life-expectancy (LE); different utilities by sex
    for(j in c("female","male")){
      if(j == "female") {fem=TRUE} else {fem=FALSE}
      
      # discounting
      LE_UK[[paste0(j, "_QALE_disc", sep="")]] <- 
        LE_UK[, j] %>% unlist() %>% 
        purrr::map_dbl(function(x) discount_LY(x, x, comorb_reduc, discount_rate_b, female=fem) )
    }
    
    # mean LEs total by age
    LE_UK <- LE_UK %>%
      transmute(group  = age,
                undisc = (female                 +male)    /2, 
                disc   = (female_QALE_disc +male_QALE_disc)/2 )
    
    
    if(run_determ==TRUE){
      
      input_val <- data.frame("group" = LE_UK$group, "mean" = LE_UK$undisc)
      
    } else {
      
      # input data by agegroups
      input_val <- LE_UK$undisc %>%
        setNames(LE_UK$group) %>%
        purrr::map_dfc(function(x) qnorm(LHS_samples[, n_LHS], x*1, x*0.1) ) %>%
        tidyr::gather(group, mean) %>%
        group_by(group) %>%
        mutate(run = 1:n_iter)
      
      n_LHS <- n_LHS +1
    }
    
    # subset outcomes needed
    output_val <- epi_dataset[compartment == "death_o"][order(scenario, run,  year, group, compartment)]
    
    # calculate results if epi iterations equal econ iterations
    if(identical(length(output_val$value), length(input_val$mean))){
      output_val             <- suppressMessages(full_join(output_val, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val$undisc      <- output_val$value * output_val$mean
      output_val$mean        <- NULL
      output_val$disc        <- output_val$undisc * (1+discount_rate_b)^-(output_val$year-1)
      output_val$compartment <- "LEs_mort"
      
    } else { # different approach if epi results different from econ iterations
      
      output_val2             <- output_res2 # take dataframe by econ iterations
      output_val2$compartment <- "LEs_mort"
      output_val2$undisc      <- rep(output_val$value, n_iter)
      output_val2             <- suppressMessages(full_join(output_val2, input_val)) # suppress as either join by run, or run and AgeGroup
      output_val2$undisc      <- output_val2$undisc * output_val2$mean
      output_val2$mean        <- NULL
      output_val2$disc        <- output_val2$undisc * (1+discount_rate_b)^-c(output_val$year-1)
      output_val              <- output_val2[, lapply(.SD, sum, na.rm=TRUE), by = .(scenario, run, compartment), 
                                             .SDcols=c("undisc", "disc") ]
    }
    
    # save
    if(save_by_agegr==TRUE){
      pr_res[["LEs_mort"]] <- output_val2 
    } else {
      pr_res[["LEs_mort"]] <- output_val
    }
  }
  
  ##############################
  # print costs and QALYs
  ##############################
  
#  # save and return results (for further analysis)
  pr_res <- data.table::rbindlist(pr_res, fill=TRUE)
  
  # add total costs/price of public/government vaccination R&D (added to the total, and at NPV)
  if(addRD_costs==TRUE){
    
    df_vaccRD <- data.table(group       = rep(rep(params$pop[[1]]$group_names, n_scens), n_iter),
                            run         = rep(1:n_iter, each=length(rep(unique(epi_dataset$scenario), each=n_agegrs))),
                            compartment = c("costs_vaccRD"),
                            scenario    = rep(rep(unique(epi_dataset$scenario), each=n_agegrs), n_iter),
                            undisc      = rep( c(rep(0, n_agegrs),  rep(rep(cost_vaccRD/n_agegrs, n_agegrs), n_scens-1)), n_iter),
                            year        = c(1),
                            disc        = rep( c(rep(0, n_agegrs),  rep(rep(cost_vaccRD/n_agegrs, n_agegrs), n_scens-1)), n_iter) )
    
    pr_res <- df_vaccRD %>%
      merge(., pr_res, sort=FALSE, all=TRUE) %>%
      ungroup()
  }
  
  # add total costs/price of public/government vaccination R&D (added to the total, and at NPV)
  if(add_Freezercost==TRUE){
    
    df_vaccFreezer <- data.table(group       = rep(rep(params$pop[[1]]$group_names, n_scens), n_iter),
                                 run         = rep(1:n_iter, each=length(rep(unique(epi_dataset$scenario), each=n_agegrs))),
                                 compartment = c("costs_vaccFreezers"),
                                 scenario    = rep(rep(unique(epi_dataset$scenario), each=n_agegrs), n_iter),
                                 undisc      = rep( c(rep(rep(0, n_agegrs), n_scens-1),  rep(rep(cost_vaccFreezers/n_agegrs, n_agegrs), n_scens-2)), n_iter),
                                 year        = c(1),
                                 disc        = rep( c(rep(rep(0, n_agegrs), n_scens-1),  rep(rep(cost_vaccFreezers/n_agegrs, n_agegrs), n_scens-2)), n_iter) )
    
    pr_res <- df_vaccFreezer %>%
      merge(., pr_res, sort=FALSE, all=TRUE) %>%
      ungroup()
  }
  
  return( as.data.table(pr_res) )
}
