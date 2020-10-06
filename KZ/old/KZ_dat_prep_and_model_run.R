## ----render, eval=FALSE, message=FALSE, warning=FALSE, include=FALSE, results='hide'----------------------------------
## rmarkdown::render(here::here('KZ','KZ_dat_prep_and_model_run.Rmd'),
##                   output_file = "README.md")
## 
## knitr::purl(input=here::here("KZ",'KZ_dat_prep_and_model_run.Rmd'),
##             output=here::here("KZ",'KZ_dat_prep_and_model_run.R'))


## ----Load packages and data, results='hide', message=FALSE, warning=FALSE---------------------------------------------
library(ggmcmc)
library(jagsUI)
library(knitr)
library(gt)
library(ggallin)
library(hrbrthemes)
library(tidyverse)
library(RColorBrewer)
library(here)

##define fn
gm_mean = function(a){prod(a)^(1/length(a))}


# Read in raw data
adult_female_survival <- read.csv(here::here("data", "KZ", "adult_female_survival_KZ.csv"))
adult_female_recruit <- read.csv(here::here("data", "KZ", "adult_female_recruit_KZ.csv"))
adult_sex_ratio <- read.csv(here::here("data", "KZ", "adult_sexratio_KZ_withimm.csv"))

adult_female_mincount <- read.csv(here::here("data", "KZ", "adult_female_mincount_KZ_withimm.csv"))
calf_mincount <- read.csv(here::here("data", "KZ", "calf_mincount_KZ.csv"))
sa1_mincount <- read.csv(here::here("data", "KZ", "sa1_mincount_KZ.csv"))

adult_female_pencount_in <-  read.csv(here::here("data", "KZ", "adult_female_pencount_in_KZ.csv"))
subadult_female_pencount_in <-  read.csv(here::here("data", "KZ", "sa1_pencount_in_KZ.csv"))
calf_female_pencount_recruit <-  read.csv(here::here("data", "KZ", "calf_female_pencount_recruit_KZ.csv"))


## ----Prep for IPM, results='hide', message=FALSE, warning=FALSE-------------------------------------------------------
#  Years of study
yrs <-  seq(from = 1996, to = 2020, by = 1)
nyr <- length(yrs)
yr_idx <- seq(from = 1, to = nyr, by = 1)
yr_df <- as.data.frame(cbind(yrs, yr_idx))

#  Years to simulate data
sim_yr_df <- yr_df %>%
	dplyr::filter(yrs >= 2014)
sim_yrs <- sim_yr_df$yr_idx
nsimyr <- length(sim_yrs)

sim_yrs2 <- yr_df %>%
	dplyr::filter(yrs >= 2013)%>%
  pull(yr_idx)

	
#  Wolf and pen indicator per year
wolfpen_ind <- adult_female_mincount$Wolf_control_pen


#  Create vectors of known numbers of adult females and calves going into pen
nPenInF <- adult_female_pencount_in$Mean

nPenInSaF <- subadult_female_pencount_in$Mean

nPenC <- calf_female_pencount_recruit$Mean
adult_sex_ratio <- adult_sex_ratio$Mean

#  Adjust minimum count of calves from control (not penned) population unit to 50% female:male ratio
female_calf_mincount <- calf_mincount %>%
	mutate(Mean = Mean * 0.5) # partial animals okay
#  Mean minimum coun of calves from not penned population
mean_calf_mincount <- mean(female_calf_mincount$Mean, na.rm = TRUE)

#  Combined minimum count of female calves across penned and control (not penned) population units
combined_female_calf_mincount <- female_calf_mincount %>%
    rbind(calf_female_pencount_recruit) %>%
    group_by(Year) %>%
    mutate(Mean = sum(Mean)) %>%
    slice(1) %>%
    as.data.frame() %>%
    mutate(Pop = NA, Pop_desc = NA)
mean_combined_female_calf_mincount <- mean(combined_female_calf_mincount$Mean, na.rm = TRUE)


#  Mean minimum count of adult females combined across penned 
#   and control (not penned) population units
mean_combined_female_adult_mincount <- mean(adult_female_mincount$Mean, na.rm = TRUE)


#  Prepare values for priors for survival and recruitment from mean data values

#  Recruitment - only need recruitment of control (not penned) population unit
adult_female_recruit <- adult_female_recruit%>%
  mutate(
	  Mean=Mean/2,
	  SD=SD/2)

meanr <- array(NA, c(3,1,2)) 
meanr[3,1,1] <- mean(adult_female_recruit$Mean, na.rm = TRUE)
meanr[3,1,2] <- 100

# Survival - need estimates for both penned and control (not penned) population units
#  Survival rate is assumed the same across subadult 1, subadult 2 and adult
means <- array(NA, c(3,2,2))

adult_female_survival_control <- adult_female_survival %>%
	dplyr::filter(Pop == 1) 
adult_female_survival_pen <-adult_female_survival %>%
	dplyr::filter(Pop == 2) 

means[2:3,1,1] <- mean(adult_female_survival_control$Mean, na.rm = TRUE)
means[2:3,2,1] <- mean(adult_female_survival_pen$Mean, na.rm = TRUE)
means[2:3,1,2] <- 100
means[2:3,2,2] <- 100


# Starting population size of control (not penned) population unit (vector of values for each age class)
n1 <- numeric()
n1[1] <- combined_female_calf_mincount$Mean[1]
n1[2] <- adult_female_mincount$Mean[1]*.15  # using year 2 bcause yar 1 is NA; arbitrary proportion of sub adults 1
n1[3] <- adult_female_mincount$Mean[1]*.85  # using year 2 bcause yar 1 is NA; arbitrary proportion of adults 


#  Format data for use in JAGS model

#  Minium counts of all individuals (only females in this set up) combined over penned and control (not penned) 
#   population units 
mindat <- adult_female_mincount %>%
	rbind(combined_female_calf_mincount) %>%
	group_by(Year) %>%
	mutate(sum = sum(Mean, na.rm = TRUE)) %>%
	slice(1) %>%
	as.data.frame() %>%
	mutate(dau = 1, 
		age = NA,
		pop = NA,
		mu = round(sum),
		tau = NA) %>%
	dplyr::filter(mu > 4) %>%
	left_join(yr_df, by = c("Year" = "yrs")) %>%
	dplyr::select(dau, yr = yr_idx, age, pop, mu, tau) 
nmin <- nrow(mindat)

#  Survival estimates by year of penned and control (not penned) population units separately
sdat <- rbind(
	adult_female_survival %>%
	mutate(dau = 1, 
		age = 3,
		pop = Pop,
		mu = Mean,
		tau = 1/ (SD * SD)) %>%
	dplyr::filter(!is.na(mu)) %>%
	left_join(yr_df, by = c("Year" = "yrs")) %>%
	dplyr::select(dau, yr = yr_idx, age, pop, mu, tau),
	adult_female_survival %>%
	mutate(dau = 1, 
		age = 2,
		pop = Pop,
		mu = Mean,
		tau = 1/ (SD * SD)) %>%
	dplyr::filter(!is.na(mu)) %>%
	left_join(yr_df, by = c("Year" = "yrs")) %>%
	dplyr::select(dau, yr = yr_idx, age, pop, mu, tau))
sdat$tau[sdat$tau == "Inf"] <- 100
ns <- nrow(sdat)	
	
#  Recruitment estimates by year of control (not penned) population unit
rdat <- adult_female_recruit  %>%
	mutate(
	  dau = 1, 
		age = NA,
		pop = Pop,
		mu = Mean,
		tau = 1/ (SD * SD)) %>%
	dplyr::filter(!is.na(mu)) %>%
	left_join(yr_df, by = c("Year" = "yrs")) %>%
	dplyr::select(dau, yr = yr_idx, age, pop, mu, tau) 
nr <- nrow(rdat)	


## ----Run model in JAGS, results='hide', message=FALSE, warning=FALSE--------------------------------------------------
# MCMC params
nt <- 1
nb <- 8000
nc <- 3
nad <- 1000
ni <- 20000


#  Gather data inputs in a list
ipm_dat <- list(nyr = nyr,
	nsimyr = nsimyr,
	nmin = nmin,
	ns = ns,
	nr = nr,
	wolfpen_ind = wolfpen_ind,
	means = means,
	meanr = meanr,
	n1 = n1, 
	nPenC = nPenC,
	nPenInF = nPenInF,
	nPenInSaF = nPenInSaF,
	adult_sex_ratio = adult_sex_ratio,
	sim_yrs = sim_yrs,
	sim_yrs2 = sim_yrs2,
	mindat = mindat, 
	sdat = sdat, 
	rdat = rdat)


#  Initial values for N to avoid parent node erros
Nst <- array(10, c(nyr,3,2))
Nst[,,2] <- NA
ipm_inits <- function(){ 
		list(N = Nst)}
	
	
#  Model parameters to monitor
model_parms <- c("lambda","logla", "lambdaC", 
	"totN", "totC", "totP", 
	"N", "S", "R", "r",
	"wolfpen_eff_r", "wolfpen_eff_s",
	"sa_yr", "ssa_yr", "r_yr", 
	"totCalves", "totCalvesC", "totCalvesP",
	"totAdults", "totAdultsC", "totAdultsP",
	"totAdultsMF", "totCalvesMF", "totNMF",
	"geom_mean_lambda",
	"geom_mean_lambda_prepen", "geom_mean_lambda_postpen","geom_mean_lambda_postwolf",
	"diff_geom_mean_lambda_post_to_pre", "diff_geom_mean_lambda_post_to_pre_iwolf",
	"geom_mean_lambda_SimPen","geom_mean_lambda_SimC", "geom_mean_lambda_SimC2",
	"geom_mean_lambda_postpen_iWolf", "geom_mean_lambda_SimPen_iWolf", "geom_mean_lambda_SimC_iWolf",
	"simN", "simTotC", "simTotP", "simTotBAU", "simTotCMF", "simTotPMF", "simTotBAUMF",
	"mean_surv_pre","mean_surv_pen","mean_surv_wolf",
	"mean_surv_presa1", "mean_surv_presa2",
	"mean_r_pre","mean_r_pen","mean_r_wolf",
	"mean_r3_pre","mean_r3_pen","mean_r3_wolf",
	"mean_r3_C", "mean_r3_P",
	"wolf_eff_proportion","pen_eff_proportion",
	"R.sim", "lambdaSimC", "lambdaSimC2","lambdaSimP")



#  Run model with survival and recruitment varying per year
#  Adjust n.cores for your computer
out_rnd_eff <- jagsUI::jags(ipm_dat, 
		inits = ipm_inits,
		model_parms,
		model.file = here::here("jags","KZ_IPM_rnd_eff_CL.txt"),
		n.chains = nc, 
		n.iter = ni,
		n.burnin = nb,
		n.thin = nt,
		n.adapt = nad,
		parallel = TRUE,
		n.cores = 7)	

#save(out_rnd_eff, file = "out_rnd_eff.Rdata")


## ----Plot results-abundance, results='hide', message=FALSE, warning=FALSE---------------------------------------------
res_df <- data.frame(rbind(yr_df,yr_df),
                     est=c(out_rnd_eff$mean$totN, out_rnd_eff$mean$lambda),
                     q2.5=c(out_rnd_eff$q2.5$totN, out_rnd_eff$q2.5$lambda),
                     q97.5=c(out_rnd_eff$q97.5$totN, out_rnd_eff$q97.5$lambda),
                     param=rep(c("MinCount_F", "Lambda"), each=nyr))

res_df <- res_df[-(nyr+1),]

ggplot(res_df,aes(x = yrs, y = est, ymin=q2.5, ymax=q97.5, fill=param)) +
  geom_cloud(steps=20, max_alpha = 1,se_mult=1.96)+
  geom_line() +
  geom_point() +
  theme_ipsum()+
  theme(legend.position = "none")+
  ylab("Estimate")+
  xlab("Year")+
  facet_wrap(vars(param), scales="free_y")+
  geom_vline(xintercept = 2012.5)+
  geom_vline(xintercept = 2013.5,linetype="dashed")+
  labs(x="Year",title="Population Growth Trajectory")



res_df%>%rbind(data.frame(rbind(yr_df),
                          est=c(out_rnd_eff$mean$totNMF),
                          q2.5=c(out_rnd_eff$q2.5$totNMF),
                          q97.5=c(out_rnd_eff$q97.5$totNMF),
                          param=rep(c("MinCount_All"), each=nyr)))%>%
  filter(param!="Lambda")%>%
  ggplot(aes(x = yrs, y = est, ymin=q2.5, ymax=q97.5, fill=param)) +
  geom_cloud(steps=20, max_alpha = 1,se_mult=1.96)+
  geom_line() +
  geom_point() +
  theme_ipsum()+
  theme(legend.position = "none")+
  ylab("Estimate")+
  xlab("Year")+
  facet_wrap(vars(param), scales="free_y")+
  geom_vline(xintercept = 2013.5)+
  geom_vline(xintercept = 2014.5,linetype="dashed")+
    labs(x="Year",title="Population Abundance and Trajectory")+
  expand_limits(y=0)

#ggsave(here::here("plots", "abundance_MF.png"), width=5, height=8)


## ----Plot results-vital rates, results='hide', message=FALSE, warning=FALSE-------------------------------------------
#R
pop_df_r <- data.frame(estimate=c(out_rnd_eff$mean$totCalvesC/out_rnd_eff$mean$totAdultsC,
                                  out_rnd_eff$mean$totCalvesP/out_rnd_eff$mean$totAdultsP),
                       lower=c(out_rnd_eff$q2.5$totCalvesC/out_rnd_eff$q2.5$totAdultsC,
                                  out_rnd_eff$q2.5$totCalvesP/out_rnd_eff$q2.5$totAdultsP),
                       upper=c(out_rnd_eff$q97.5$totCalvesC/out_rnd_eff$q97.5$totAdultsC,
                                  out_rnd_eff$q97.5$totCalvesP/out_rnd_eff$q97.5$totAdultsP),
                       pop=rep(c("Free", "Pen"), each=nyr),
                       param="Recruitment",
                       yrs=rep(yrs, times=2))%>%
  mutate(lower=case_when(lower<0~0,
                         TRUE~lower))

#S
pop_df_s <- data.frame(estimate=c(out_rnd_eff$mean$S[,3,1],
                                  out_rnd_eff$mean$S[,3,2]),
                       lower=c(out_rnd_eff$q2.5$S[,3,1],
                               out_rnd_eff$q2.5$S[,3,2]),
                       upper=c(out_rnd_eff$q97.5$S[,3,1],
                               out_rnd_eff$q97.5$S[,3,2]),
                       pop=rep(c("Free", "Pen"), each=nyr),
                       param="Survival",
                       yrs=rep(yrs, times=2))



# pop_df <- data.frame(estimate=pop_df_s$estimate/(1-pop_df_r$estimate),
#                      lower=pop_df_s$lower/(1-pop_df_r$lower),
#                      upper=pop_df_s$upper/(1-pop_df_r$upper),
#                      pop=rep(c("Free", "Pen"), each=nyr),
#                      yrs=rep(yrs, times=2))

pop_df <- data.frame(estimate=pop_df_s$estimate+pop_df_r$estimate,
                     lower=pop_df_s$lower+pop_df_r$lower,
                     upper=pop_df_s$upper+pop_df_r$upper,
                     pop=rep(c("Free", "Pen"), each=nyr),
                     yrs=rep(yrs, times=2))

pop_df2 <- rbind(pop_df_r, pop_df_s)

ggplot(pop_df2%>%
         filter(yrs>1996)%>%
         mutate(estimate=case_when(param=="Recruitment" & yrs<2015 & pop=="Pen"~NA_real_, 
                                   param=="Survival" & yrs<2014 & pop=="Pen"~NA_real_, 
                                   TRUE~estimate)),
       aes(x = yrs, y = estimate, fill=pop)) +
  geom_line(aes(color=pop)) +
  theme_ipsum()+
  theme_ipsum()+
  labs(x="Year", y="Rate", title="Annual vital rates for penned and free animals")+
  facet_wrap(vars(param), scales="free_y")+
    geom_vline(xintercept = 2013)+
  geom_vline(xintercept = 2014,linetype="dashed")

#ggsave(here::here("plots", "vital_rates.png"), width=8, height=5)


## ----Plot results-vital rates vs raw data, results='hide', message=FALSE, warning=FALSE-------------------------------
calc.vr <- pop_df2%>%
         filter(yrs>1996)%>%
         mutate(estimate=case_when(param=="Recruitment" & yrs<2015 & pop=="Pen"~NA_real_, 
                                   param=="Survival" & yrs<2014 & pop=="Pen"~NA_real_, 
                                   TRUE~estimate))

raw.vr <- rbind(adult_female_survival%>%mutate(param="Survival",
                                         pop=case_when(Pop_desc%in%"Control"~"Free",TRUE~as.character(Pop_desc)),
                                         estimate=Mean,
                                         lower=NA,
                                         upper=NA,
                                         yrs=Year)%>%
  select(estimate,lower,upper,pop,param,yrs),  ###ADD AFS
  
  adult_female_recruit%>%mutate(param="Recruitment",
                                         pop=case_when(Pop_desc%in%"Control"~"Free",TRUE~as.character(Pop_desc)),
                                         estimate=Mean,
                                         lower=NA,
                                         upper=NA,
                                         yrs=Year)%>%
  select(estimate,lower,upper,pop,param,yrs),  ###ADD FREE RECRUIT
  
  calf_female_pencount_recruit%>%
    left_join(adult_female_pencount_endofyr%>%select(Year,F.cnt=Mean))%>%
    mutate(estimate=Mean/F.cnt)%>%
    filter(Year>2014)%>% 
                                mutate(param="Recruitment",
                                         lower=NA,
                                         upper=NA,
                                         yrs=Year,
                                       pop=Pop_desc)%>%
  select(estimate,lower,upper,pop,param,yrs))  ###ADD PEN RECRUIT
                                 
  
calc.vr%>%mutate(class="modelled")%>%
  rbind(raw.vr%>%mutate(class="raw"))%>%
ggplot(aes(x = yrs, y = estimate, fill=class)) +
  geom_line(aes(color=class)) +
  theme_ipsum()+
  theme_ipsum()+
  labs(x="Year", y="Rate", title="Comparing modelled vs raw vital rate")+
  facet_wrap(vars(param,pop), scales="free_y")+
    geom_vline(xintercept = 2013.5)+
  geom_vline(xintercept = 2014.5,linetype="dashed")

#ggsave(here::here("plots", "vital_rates.png"), width=8, height=5)


raw.vr%>%
  filter(pop=="Pen" & yrs>2013)%>%
  arrange(yrs,param)%>%
  mutate(l=case_when(param=="Recruitment"~lag(estimate)/(1-estimate),
                     TRUE~NA_real_))%>%
  drop_na(l)%>%
  pull(l)%>%
  gm_mean()


calc.vr%>%
  filter(pop=="Pen" & yrs>2013)%>%
  arrange(yrs,param)%>%
  mutate(l=case_when(param=="Recruitment"~lag(estimate)/(1-estimate),
                     TRUE~NA_real_))%>%
  drop_na(l)%>%
  pull(l)%>%
  gm_mean()


## ----Plot results-counts vs raw data, fig.height=6, fig.width=6, message=FALSE, warning=FALSE, results='hide'---------

calc.abund <- res_df%>%rbind(data.frame(rbind(yr_df),
                          est=c(out_rnd_eff$mean$totNMF),
                          q2.5=c(out_rnd_eff$q2.5$totNMF),
                          q97.5=c(out_rnd_eff$q97.5$totNMF),
                          param=rep(c("MinCount_All"), each=nyr)))%>%
  filter(param=="MinCount_F")


raw.abund <-rbind(adult_female_mincount%>%mutate(param="MinCount_F",
                                         est=Mean,
                                         q2.5=NA,
                                         q97.5=NA,
                                         yrs=Year,
                                         yr_idx=NA)%>%
  select(colnames(calc.abund)),  ###ADD AF Count
  
        calf_mincount%>%mutate(param="MinCount_Calf",
                                         est=Mean/2,
                                         q2.5=NA,
                                         q97.5=NA,
                                         yrs=Year,
                                         yr_idx=NA)%>%
  select(colnames(calc.abund)),  ###ADD WILD CALVES
  
          calf_female_pencount_recruit%>%mutate(param="MinCount_Calf",
                                         est=Mean,
                                         q2.5=NA,
                                         q97.5=NA,
                                         yrs=Year,
                                         yr_idx=NA)%>%
  select(colnames(calc.abund)) ###ADD PEN CALVES
  )%>%
  group_by(yrs,yr_idx)%>%
  summarise(est=sum(est),
            q2.5=sum(q2.5),
            q97.5=sum(q97.5))%>%
mutate(param="MinCount_F")
         




calc.abund%>%mutate(class="modelled")%>%
  rbind(raw.abund%>%as.data.frame()%>%mutate(class="observed"))%>%
ggplot(aes(x = yrs, y = est, fill=class)) +
  geom_cloud(aes( ymin=q2.5, ymax=q97.5),steps=20, max_alpha = 1,se_mult=1.96)+
  geom_line(aes(color=class)) +
  geom_point(aes(color=class)) +
  theme_ipsum()+
  theme_ipsum()+
  labs(x="Year", y="Abundance", title="Comparing modelled vs observed female counts")+
    geom_vline(xintercept = 2013.5)+
  geom_vline(xintercept = 2014.5,linetype="dashed")+
  expand_limits(y=0)




## ----Popsim, fig.height=6, fig.width=6, results='asis'----------------------------------------------------------------
pop.sim <-data.frame(yrs=rep(c(2014:2020), times=3),
                     est=c(out_rnd_eff$mean$simTotC, out_rnd_eff$mean$simTotP, out_rnd_eff$mean$simTotBAU),
                     q2.5=c(out_rnd_eff$q2.5$simTotC, out_rnd_eff$q2.5$simTotP, out_rnd_eff$q2.5$simTotBAU),
                     q97.5=c(out_rnd_eff$q97.5$simTotC, out_rnd_eff$q97.5$simTotP, out_rnd_eff$q97.5$simTotBAU),
                     pop=rep(c("Free", "Pen", "Control"), each=length(out_rnd_eff$q97.5$simTotC)))



ggplot(pop.sim%>%mutate(pop=fct_relevel(pop,"Pen","Free","Control")),
       aes(x = yrs, y = est, ymin=q2.5, ymax=q97.5, fill=pop, linetype=pop)) +
  geom_cloud(steps=20, max_alpha = 0.8,se_mult=1.96)+
  geom_line() +
  theme_ipsum()+
  labs(x="Year", y="Abundance", title="Simulated Female Abundance", subtitle="Same starting abundance and treatment applied to all females")+
  expand_limits(y=0)


## ----Age structure, fig.height=6, fig.width=6, results='asis'---------------------------------------------------------

as.data.frame(out_rnd_eff$mean$N[,,1])%>%
  mutate(pop="Free",
         year=yrs)%>%
  pivot_longer(-c("pop","year"))%>%
  mutate(name=str_sub(name,2,-1))%>%
rbind(
  as.data.frame(out_rnd_eff$mean$N[,,2])%>%
  mutate(pop="Pen",
         year=yrs)%>%
  pivot_longer(-c("pop","year"))%>%
  mutate(name=str_sub(name,2,-1))
)%>%
  ggplot(aes(x=year, y=value, fill = name)) +
  geom_area()+
  ylab("Abundance")+
  xlab("Year")+
  facet_wrap(vars(pop))+
    theme_ipsum()


as.data.frame(out_rnd_eff$mean$N[,,1])%>%
  mutate(pop="Free",
         year=yrs)%>%
  pivot_longer(-c("pop","year"))%>%
  mutate(name=str_sub(name,2,-1))%>%
rbind(
  as.data.frame(out_rnd_eff$mean$N[,,2])%>%
  mutate(pop="Pen",
         year=yrs)%>%
  pivot_longer(-c("pop","year"))%>%
  mutate(name=str_sub(name,2,-1))
)%>%
  group_by(pop,year)%>%
  mutate(sum=sum(value),
         value2=value/sum)%>%
  
  ggplot(aes(x=year, y=value2, fill = name)) +
  geom_area()+
  ylab("Proportion in each ageclass")+
  xlab("Year")+
  facet_wrap(vars(pop), scales="free_x")+
    theme_ipsum()


as.data.frame(out_rnd_eff$mean$N[,,1])%>%
  mutate(pop="Free",
         year=yrs)%>%
  pivot_longer(-c("pop","year"))%>%
  mutate(name=str_sub(name,2,-1))%>%
rbind(
  as.data.frame(out_rnd_eff$mean$N[,,2])%>%
  mutate(pop="Pen",
         year=yrs)%>%
  pivot_longer(-c("pop","year"))%>%
  mutate(name=str_sub(name,2,-1))
)%>%
  filter(year>2013)%>%
  ggplot(aes(x=year, y=value, fill = name)) +
  geom_area()+
  ylab("Abundance")+
  xlab("Year")+
  facet_wrap(vars(pop))+
    theme_ipsum()


as.data.frame(out_rnd_eff$mean$N[,,1])%>%
  mutate(pop="Free",
         year=yrs)%>%
rbind(
  as.data.frame(out_rnd_eff$mean$N[,,2])%>%
  mutate(pop="Pen",
         year=yrs)
)%>%
  write_csv(here::here("ages_forScott.csv"))



## ----Plot results-summarise growth rates, results='asis'--------------------------------------------------------------

summary.l <- tribble(
  ~pop,~l, ~l.lower, ~l.upper,
"pre-mgmt",out_rnd_eff$mean$geom_mean_lambda_prepen, out_rnd_eff$q2.5$geom_mean_lambda_prepen,out_rnd_eff$q97.5$geom_mean_lambda_prepen,
"post-mgmt", out_rnd_eff$mean$geom_mean_lambda_postpen, out_rnd_eff$q2.5$geom_mean_lambda_postpen, out_rnd_eff$q97.5$geom_mean_lambda_postpen,
"post-mgmt_intensewolf", out_rnd_eff$mean$geom_mean_lambda_postpen_iWolf, out_rnd_eff$q2.5$geom_mean_lambda_postpen_iWolf, out_rnd_eff$q97.5$geom_mean_lambda_postpen_iWolf,
"Free_intensewolf", out_rnd_eff$mean$geom_mean_lambda_SimC_iWolf, out_rnd_eff$q2.5$geom_mean_lambda_SimC_iWolf,out_rnd_eff$q97.5$geom_mean_lambda_SimC_iWolf,
"Pen_intensewolf", out_rnd_eff$mean$geom_mean_lambda_SimPen_iWolf, out_rnd_eff$q2.5$geom_mean_lambda_SimPen_iWolf, out_rnd_eff$q97.5$geom_mean_lambda_SimPen_iWolf,
"Free",out_rnd_eff$mean$geom_mean_lambda_SimC2, out_rnd_eff$q2.5$geom_mean_lambda_SimC2, out_rnd_eff$q97.5$geom_mean_lambda_SimC2,
"Pen",out_rnd_eff$mean$geom_mean_lambda_SimPen, out_rnd_eff$q2.5$geom_mean_lambda_SimPen, out_rnd_eff$q97.5$geom_mean_lambda_SimPen
)%>%
  mutate_if(is.numeric,function(x) round(x,2))


summary.l$Years <- c("1996-2013", "2015-2020","2017-2020","2017-2020", "2017-2020", "2014-2020","2015-2020")
colnames(summary.l) <- c("Group", "Lambda", "Lamba.Lower", "Lambda.Upper", "Years")


summary.l <- summary.l%>%
  mutate(`95% CI`=paste(Lamba.Lower,Lambda.Upper, sep="-"))%>%
  select(Group, Years, Lambda,`95% CI`)

gt(summary.l)%>%
    tab_header(
    title = md("Population Growth Rates")
  ) 


data.frame(rbind(yr_df),
           est=c(out_rnd_eff$mean$totNMF),
           q2.5=c(out_rnd_eff$q2.5$totNMF),
           q97.5=c(out_rnd_eff$q97.5$totNMF),
           param=rep(c("MinCount_All"), each=nyr))%>%
  mutate(l = est/lag(est))%>%
  filter(yrs>2013)%>%
  drop_na(l)%>%
  summarise(l.all=gm_mean(l)%>%round(2))


mean_r4_C <- c()
mean_r4_P <- c()
q2.5_r4_C <- c()
q2.5_r4_P <- c()
q97.5_r4_C <- c()
q97.5_r4_P <- c()
mean_r4_C[1] <- 1
mean_r4_P[1] <- 1
q2.5_r4_C[1] <- 1
q2.5_r4_P[1] <- 1
q97.5_r4_C[1] <- 1
q97.5_r4_P[1] <- 1
for(y in 2:nyr){
# mean_r4_C[y] <- out_rnd_eff$mean$N[y,1,1]/((out_rnd_eff$mean$N[y-1,4,1]*out_rnd_eff$mean$S[y-1,4,1])+(out_rnd_eff$mean$N[y-1,3,1]*out_rnd_eff$mean$S[y-1,3,1]))
# mean_r4_P[y] <- out_rnd_eff$mean$N[y,1,2]/((out_rnd_eff$mean$N[y-1,4,2]*out_rnd_eff$mean$S[y-1,4,2])+(out_rnd_eff$mean$N[y-1,3,2]*out_rnd_eff$mean$S[y-1,3,2]))
mean_r4_C[y] <- out_rnd_eff$mean$N[y,1,1]/(out_rnd_eff$mean$N[y-1,3,1]*out_rnd_eff$mean$S[y-1,3,1])
mean_r4_P[y] <- out_rnd_eff$mean$N[y,1,2]/(out_rnd_eff$mean$N[y-1,3,2]*out_rnd_eff$mean$S[y-1,3,2])

q2.5_r4_C[y] <- out_rnd_eff$q2.5$N[y,1,1]/(out_rnd_eff$mean$N[y-1,3,1]*out_rnd_eff$mean$S[y-1,3,1])
q2.5_r4_P[y] <- out_rnd_eff$q2.5$N[y,1,2]/(out_rnd_eff$mean$N[y-1,3,2]*out_rnd_eff$mean$S[y-1,3,2])

q97.5_r4_C[y] <- out_rnd_eff$q97.5$N[y,1,1]/(out_rnd_eff$mean$N[y-1,3,1]*out_rnd_eff$mean$S[y-1,3,1])
q97.5_r4_P[y] <- out_rnd_eff$q97.5$N[y,1,2]/(out_rnd_eff$mean$N[y-1,3,2]*out_rnd_eff$mean$S[y-1,3,2])
}

mean(mean_r4_C[2:18])
mean(mean_r4_P[20:25])
mean(mean_r4_C[19:25])


out_rnd_eff$mean$wolf_eff_proportion

out_rnd_eff$mean$pen_eff_proportion


out_rnd_eff$mean$mean_r3_pen


## ----Plot results-summarise vital rates, results='asis'---------------------------------------------------------------

summary.s <- tribble(
  ~pop,~s, ~s.lower, ~s.upper,
"pre-mgmt",out_rnd_eff$mean$mean_surv_pre, out_rnd_eff$q2.5$mean_surv_pre,out_rnd_eff$q97.5$geom_mean_lambda_pre,
"pen", out_rnd_eff$mean$mean_surv_pen, out_rnd_eff$q2.5$mean_surv_pen, out_rnd_eff$q97.5$mean_surv_pen,
"wolf", out_rnd_eff$mean$mean_surv_wolf, out_rnd_eff$q2.5$mean_surv_wolf, out_rnd_eff$q97.5$mean_surv_wolf)

summary.r <- tribble(
  ~pop,~r, ~r.lower, ~r.upper,
"pre-mgmt",out_rnd_eff$mean$mean_r_pre, out_rnd_eff$q2.5$mean_r_pre,out_rnd_eff$q97.5$geom_mean_lambda_pre,
"pen", out_rnd_eff$mean$mean_r_pen, out_rnd_eff$q2.5$mean_r_pen, out_rnd_eff$q97.5$mean_r_pen,
"wolf", out_rnd_eff$mean$mean_r_wolf, out_rnd_eff$q2.5$mean_r_wolf, out_rnd_eff$q97.5$mean_r_wolf)

summary.r3 <- tribble(
  ~pop,~r3, ~r3.lower, ~r3.upper,
"pre-mgmt",mean(mean_r4_C[2:18]), mean(q2.5_r4_C[2:18]), mean(q97.5_r4_C[2:18]),
"pen", mean(mean_r4_P[20:25]), mean(q2.5_r4_P[20:25]), mean(q97.5_r4_P[20:25]),
"wolf", mean(mean_r4_C[19:25]), mean(q2.5_r4_C[19:25]), mean(q97.5_r4_C[19:25]))




summary.vr <- summary.s%>%
  left_join(summary.r)%>%
  left_join(summary.r3)%>%
  mutate_if(is.numeric,function(x) round(x,2))


summary.vr$Years <- c("1996-2013", "2015-2020", "2014-2020")

summary.vr <- summary.vr%>%
  mutate(`s 95% CI`=paste(s.lower,s.upper, sep="-"),
         `r 95% CI`=paste(r.lower,r.upper, sep="-"),
         `r3 95% CI`=paste(r3.lower,r3.upper, sep="-"))%>%
  select(pop, Years, s,`s 95% CI`, r,`r 95% CI`, r3, `r3 95% CI`)



colnames(summary.vr) <- c("Group", "Years", "AF Survival","95% CI", "Recruitment","r95% CI", "Adult F Recruitment", "r3_95% CI")


gt(summary.vr)%>%
  tab_header(
    title = md("Vital Rates")
  ) 


## ----effects, results='asis'------------------------------------------------------------------------------------------

summary.effect <- tribble(
  ~pop,~lambda.dif, ~lower, ~upper,
"pre vs post",out_rnd_eff$mean$diff_geom_mean_lambda_post_to_pre, out_rnd_eff$q2.5$diff_geom_mean_lambda_post_to_pre,out_rnd_eff$q97.5$diff_geom_mean_lambda_post_to_pre,
"pre vs post_iwolf",out_rnd_eff$mean$diff_geom_mean_lambda_post_to_pre_iwolf, out_rnd_eff$q2.5$diff_geom_mean_lambda_post_to_pre_iwolf,out_rnd_eff$q97.5$diff_geom_mean_lambda_post_to_pre_iwolf)%>%
    mutate_if(is.numeric,function(x) round(x,2))



gt(summary.effect)%>%
  tab_header(
    title = md("Treatment Effect")
  ) 

