## ----render, eval=FALSE, message=FALSE, warning=FALSE, include=FALSE, results='hide'---------------------------------------------------------
## rmarkdown::render(here::here("QT",'QT_data_prep_and_model_run.Rmd'),
##                   output_file = "README.md")
## 
## knitr::purl(input=here::here("QT",'QT_data_prep_and_model_run.Rmd'),
##             output=here::here("QT",'QT_data_prep_and_model_run.R'))


## ----Load packages and data, results='hide', message=FALSE, warning=FALSE--------------------------------------------------------------------
library(ggmcmc)
library(jagsUI)
library(knitr)
library(gt)
library(ggallin)
library(hrbrthemes)
library(tidyverse)
library(rjags)
library(MCMCvis)
library(here)

##define fn
gm_mean = function(a){prod(a)^(1/length(a))}


# Read in raw data
adult_female_survival <- read.csv(here::here("data", "QT", "adult_female_survival_QT.csv"))
adult_female_recruit <- read.csv(here::here("data", "QT", "adult_female_recruit_QT.csv"))

adult_sex_ratio <- read.csv(here::here("data", "QT", "adult_sexratio_QT.csv"))


count_dat <- read_xlsx(here::here("data", "QT", "Count_summary_QT.xlsx"))%>%
  mutate_at(c("SurveryCount_ADULTMF","SurveryCount_CALFMF", "Estimate_ADULTMF","Estimate_CALFMF"), funs(round(., 0)))

##remove 2019 est
##significant evidence that this is a poor estimate of population abundance
##min count was 34% > than the sightability corrected survey count.
count_dat[19,"Estimate_ADULTMF"] <- NA
count_dat[19,"SD_Estimate_ADULTMF"] <- NA
count_dat[19,"Estimate_CALFMF"] <- NA
count_dat[19,"SD_Estimate_CALFMF"] <- NA

count_dat[19,"SurveryCount_ADULTMF"] <- NA
count_dat[19,"SurveryCount_CALFMF"] <- NA

sightability <- read.csv(here::here("data", "QT", "sightability_QT.csv"))


## ----sight bootstrap, eval=FALSE, message=FALSE, warning=FALSE, include=FALSE, results='hide'------------------------------------------------
## 
## for(i in 1:nrow(sightability)){
##   if(!is.na(sightability[i,"seen"])){
##     missed <- sightability[i,"out"]-sightability[i,"seen"]
##     if(missed==0){missed<-1}
##     a <-data.frame(outcome=c(rep(1,times=sightability[i,"seen"]),
##                          rep(0,times=missed)))
## 
##     boot.summary <- c()
##     for(boot in 1:1000){
##   b <- a%>%sample_frac(1, replace=TRUE)%>%summarise(mean=mean(outcome))
##      boot.summary[boot]<- b[[1]]
##     }
## 
##     mean(boot.summary)
##     sightability[i,"SD"] <-sd(boot.summary)
##   }
## }
## 
## sightability$mean <- sightability$seen/sightability$out
## 
## ggplot(sightability, aes(x = Year, y = mean)) +
##     geom_errorbar(aes( ymin=mean-SD, ymax=mean+SD),alpha=0.5)+
##   geom_point() +
##   theme_ipsum()+
##   theme_ipsum()+
##   geom_hline(yintercept = 1)+
##   labs(x="Year", y="Sightability", title="Sightability during census")+
##   expand_limits(y=0)


## ----Prep for IPM, results='hide', message=FALSE, warning=FALSE------------------------------------------------------------------------------
#  Years of study
yrs <-  seq(from = 2001, to = 2020, by = 1)
nyr <- length(yrs)
yr_idx <- seq(from = 1, to = nyr, by = 1)
yr_df <- as.data.frame(cbind(yrs, yr_idx))


#  Create vectors of sex ratio
adult_sex_ratio <- adult_sex_ratio$Mean


#  Recruitment 
adult_female_recruit <- adult_female_recruit%>%	mutate(
	  Mean=Mean/2,
	  SD=SD/2)

meanr <- array(NA, c(3,1,2)) 
meanr[3,1,1] <- mean(adult_female_recruit$Mean, na.rm = TRUE)
meanr[3,1,2] <- 100


# Survival - need estimates for both penned and control (not penned) population units
#  Survival rate is assumed the same across subadult 1, subadult 2 and adult
means <- array(NA, c(3,1,2))

adult_female_survival_control <- adult_female_survival %>%
	dplyr::filter(Pop == 1) 

means[2:3,1,1] <- mean(adult_female_survival_control$Mean, na.rm = TRUE)
means[2:3,1,2] <- 100


# Starting population size of population unit (vector of values for each age class)
n1 <- numeric()
n1[1] <- count_dat$Estimate_CALFMF[2]*.5
n1[2] <- count_dat$Estimate_ADULTMF[2]*.15*count_dat$SexRatio[2]
n1[3] <- count_dat$Estimate_ADULTMF[2]*.85*count_dat$SexRatio[2]


#  Survival estimates by year of penned and control (not penned) population units separately
sdat <- adult_female_survival %>%
	mutate(dau = 1, 
		age = 3,
		pop = Pop,
		mu = ifelse(Mean!= 1, Mean, 0.99),
		tau = 1/ (SD * SD)) %>%
		#tau = 500) %>%
	dplyr::filter(!is.na(mu)) %>%
	left_join(yr_df, by = c("Year" = "yrs")) %>%
	dplyr::select(dau, yr = yr_idx, age, pop, mu, tau)
sdat$tau[sdat$tau == "Inf"] <- mean(sdat$tau[!is.infinite(sdat$tau)])
ns <- nrow(sdat)	


#  Recruitment estimates by year of control (not penned) population unit
rdat <- adult_female_recruit  %>%
	mutate(dau = 1, 
		age = NA,
		pop = Pop,
		mu = Mean,
		tau = 1/ (SD * SD)) %>%
		#tau = 0.1) %>%
	dplyr::filter(!is.na(mu)) %>%
	left_join(yr_df, by = c("Year" = "yrs")) %>%
	dplyr::select(dau, yr = yr_idx, age, pop, mu, tau) 
nr <- nrow(rdat)	


#  Minium counts of all individuals (only females in this set up) combined over penned and control (not penned) 
#   population units 
mindat <- count_dat %>%
	mutate(dau = 1, 
		age = NA,
		pop = NA,
		mu = SurveryCount_ADULTMF,
		tau = NA) %>%
	dplyr::filter(mu > 4) %>%
	left_join(yr_df, by = c("Year" = "yrs")) %>%
	dplyr::select(dau, yr = yr_idx, age, pop, mu, tau)

calf_mindat <- count_dat %>%
	mutate(dau = 1, 
		age = NA,
		pop = NA,
		mu = SurveryCount_CALFMF,
		tau = NA) %>%
	dplyr::filter(mu > 4) %>%
	left_join(yr_df, by = c("Year" = "yrs")) %>%
	dplyr::select(dau, yr = yr_idx, age, pop, mu, tau) 
nmin <- nrow(mindat)

abundat <- count_dat %>%
	mutate(dau = 1, 
		age = NA,
		pop = NA,
		mu = Estimate_ADULTMF,
		tau = 1/(SD_Estimate_ADULTMF * SD_Estimate_ADULTMF)) %>%
	dplyr::filter(mu > 4) %>%
	left_join(yr_df, by = c("Year" = "yrs")) %>%
	dplyr::select(dau, yr = yr_idx, age, pop, mu, tau) 
na <- nrow(abundat)

calf_abundat <- count_dat %>%
	mutate(dau = 1, 
		age = NA,
		pop = NA,
		mu = Estimate_CALFMF,
		tau = 1/(SD_Estimate_CALFMF * SD_Estimate_CALFMF)) %>%
	dplyr::filter(mu > 4) %>%
	left_join(yr_df, by = c("Year" = "yrs")) %>%
	dplyr::select(dau, yr = yr_idx, age, pop, mu, tau) 



## ----Run model in JAGS, results='hide', message=FALSE, warning=FALSE-------------------------------------------------------------------------
#  Gather data inputs in a list
ipm_dat <- list(nyr = nyr,
	nmin = nmin,
	ns = ns,
	nr = nr,
	na = na,
	#wolfpen_ind = wolfpen_ind,
	means = means,
	meanr = meanr,
	n1 = n1, 
	adult_sex_ratio = adult_sex_ratio,
	mindat = mindat, 
	calf_mindat = calf_mindat, 
	abundat = abundat,
	calf_abundat = calf_abundat,
	sdat = sdat, 
	rdat = rdat)


#  Initial values for N to avoid parent node erros
Nst <- array(10, c(nyr,3,1))
ipm_inits <- function(){ 
		list(N = Nst)}
	
	
#  Model parameters to monitor
model_parms <- c("lambda","logla",
	"totN", "N", "S", "R", "r1", "r2", 
	#"wolfpen_eff_r", "wolfpen_eff_s",
	"sa_yr", "r_yr", 
	"totCalves", "totAdults",	
	"totAdultsMF", "totCalvesMF", "totNMF",
	"geom_mean_lambda",
	"geom_mean_lambda_pre", "geom_mean_lambda_post",
	"diff_geom_mean_lambda_post_to_pre","diff_geom_mean_lambda_post_to_pre_iwolf",
	"geom_mean_lambda_post_iwolf",
	"mean_surv_pre", "mean_surv_post",
	"mean_r_pre", "mean_r_post",
	"mean_r3_pre", "mean_r3_post")


nt <- 3
nb <- 8000
nc <- 3
nad <- 5000
ni <- 50000

out_rnd_eff <- jagsUI::jags(ipm_dat, 
	inits = ipm_inits,
	model_parms,
	model.file = here::here("jags","QT_IPM_rnd_eff_w_abund_and_mincount_CL.txt"),
	n.chains = nc, 
	n.iter = ni,
	n.burnin = nb,
	n.thin = nt,
	n.adapt = nad)

###save output, using piggyback because too large for github, then delete file so it doesn't upload
# ##need a github token
# ##get yours here:
# browse_github_pat()
# github_token()
# ##Add to your r environ 
#  
# ##bring up R envir using this command below
# edit_r_environ()
# 
# #add this line
# GITHUB_TOKEN="***ADD TOKEN HERE***"
#
#
###OR for 1 time use add it here and run
#Sys.setenv(GITHUB_TOKEN="***ADD TOKEN HERE***")

save(out_rnd_eff, file = here::here("QT", "qt_out_rnd_eff.Rdata"))
pb_upload(here::here("QT", "qt_out_rnd_eff.Rdata"), 
          repo = "ctlamb/KZ_QT_IPM", 
          tag = "v0.0.1",
          overwrite=TRUE)
unlink(here::here("QT", "qt_out_rnd_eff.Rdata"))



## ----MODEL DIAGNOSTICS,eval=FALSE, results='hide', message=FALSE, warning=FALSE--------------------------------------------------------------
## MCMCtrace(out_rnd_eff,
##         params = "totAdult",
##         ISB = FALSE,
##         pdf = FALSE)
## 
## MCMCtrace(out_rnd_eff,
##         params = "r1",
##         ISB = FALSE,
##         pdf = FALSE)
## 
## MCMCtrace(out_rnd_eff,
##         params = "R",
##         ISB = FALSE,
##         pdf = FALSE)
## 
## 
## MCMCtrace(out_rnd_eff,
##         params = "lambda",
##         ISB = FALSE,
##         pdf = FALSE)
## 
## MCMCsummary(out_rnd_eff)%>%filter(n.eff<400)
## 


## ----Plot results-abundance, results='hide', message=FALSE, warning=FALSE--------------------------------------------------------------------
res_df <- data.frame(rbind(yr_df,yr_df),
                     est=c(out_rnd_eff$mean$totN, out_rnd_eff$mean$lambda),
                     q2.5=c(out_rnd_eff$q2.5$totN, out_rnd_eff$q2.5$lambda),
                     q97.5=c(out_rnd_eff$q97.5$totN, out_rnd_eff$q97.5$lambda),
                     param=rep(c("totN", "Lambda"), each=nyr))

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
  geom_vline(xintercept = 2016.5)+
  labs(x="Year",title="Population Growth Trajectory")



res_df%>%rbind(data.frame(rbind(yr_df),
                          est=c(out_rnd_eff$mean$totNMF),
                          q2.5=c(out_rnd_eff$q2.5$totNMF),
                          q97.5=c(out_rnd_eff$q97.5$totNMF),
                          param=rep(c("totNMF"), each=nyr)))%>%
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
  geom_vline(xintercept = 2016.5)+
    labs(x="Year",title="Population Abundance and Trajectory")+
  expand_limits(y=0)

#ggsave(here::here("plots", "abundance_MF.png"), width=5, height=8)


## ----Plot results-vital rates, results='hide', message=FALSE, warning=FALSE------------------------------------------------------------------
#R
pop_df_r <- data.frame(estimate=out_rnd_eff$mean$R,
                       lower=out_rnd_eff$q2.5$R,
                       upper=out_rnd_eff$q97.5$R,
                       pop=rep(c("Free"), each=nyr),
                       param="Recruitment",
                       yrs=rep(yrs, times=1))%>%
  mutate(lower=case_when(lower<0~0,
                         TRUE~lower))

#S
pop_df_s <- data.frame(estimate=c(out_rnd_eff$mean$S),
                       lower=c(out_rnd_eff$q2.5$S),
                       upper=c(out_rnd_eff$q97.5$S),
                       pop=rep(c("Free"), each=nyr),
                       param="Survival",
                       yrs=rep(yrs, times=1))



# pop_df <- data.frame(estimate=pop_df_s$estimate/(1-pop_df_r$estimate),
#                      lower=pop_df_s$lower/(1-pop_df_r$lower),
#                      upper=pop_df_s$upper/(1-pop_df_r$upper),
#                      pop=rep(c("Free", "Pen"), each=nyr),
#                      yrs=rep(yrs, times=2))

pop_df <- data.frame(estimate=pop_df_s$estimate+pop_df_r$estimate,
                     lower=pop_df_s$lower+pop_df_r$lower,
                     upper=pop_df_s$upper+pop_df_r$upper,
                     pop=rep(c("Free"), each=nyr),
                     yrs=rep(yrs, times=1))

pop_df2 <- rbind(pop_df_r, pop_df_s)

ggplot(data=pop_df2%>%
         filter(yrs>2000),
       aes(x = yrs, y = estimate, fill=pop,ymin=lower, ymax=upper)) +
  geom_line(aes(color=pop)) +
  #geom_cloud(aes( ymin=lower, ymax=upper),steps=20, max_alpha = 1,se_mult=1.96)+
  geom_ribbon(alpha=0.4)+
  theme_ipsum()+
  labs(x="Year", y="Rate", title="Annual vital rates")+
  facet_wrap(vars(param), scales="free_y")+
    geom_vline(xintercept = 2016.5)

#ggsave(here::here("plots", "vital_rates.png"), width=8, height=5)


## ----Plot results-vital rates vs raw data, results='hide', message=FALSE, warning=FALSE------------------------------------------------------
calc.vr <- pop_df2%>%
         filter(yrs>2001)

raw.vr <- rbind(adult_female_survival%>%mutate(param="Survival",
                                         pop=case_when(Pop_desc%in%"Control"~"Free",TRUE~as.character(Pop_desc)),
                                         estimate=Mean,
                                         lower=NA,
                                         upper=NA,
                                         yrs=Year)%>%
  dplyr::select(estimate,lower,upper,pop,param,yrs),  ###ADD AFS
  
  adult_female_recruit%>%mutate(param="Recruitment",
                                         pop=case_when(Pop_desc%in%"Control"~"Free",TRUE~as.character(Pop_desc)),
                                         estimate=Mean,
                                         lower=NA,
                                         upper=NA,
                                         yrs=Year)%>%
   dplyr::select(estimate,lower,upper,pop,param,yrs))  ###ADD FREE RECRUIT
  
          
  
calc.vr%>%mutate(class="modelled")%>%
  rbind(raw.vr%>%mutate(class="raw"))%>%
ggplot(aes(x = yrs, y = estimate, fill=class,ymin=lower, ymax=upper)) +
  #geom_cloud(steps=20, max_alpha = 1,se_mult=1.96)+
    geom_ribbon(alpha=0.4)+
  geom_line(aes(color=class)) +
  theme_ipsum()+
  theme_ipsum()+
  labs(x="Year", y="Rate", title="Comparing modelled vs raw vital rate")+
  facet_wrap(vars(param,pop), scales="free_y")+
    geom_vline(xintercept = 2016.5)

#ggsave(here::here("plots", "QT_vital_ratesfit.png"), width=8, height=5)

write_csv(raw.vr, here::here("data", "QT", "vitalrate_validation_QT.csv"))


## ----Plot results-counts vs raw data, fig.height=6, fig.width=6, message=FALSE, warning=FALSE, results='hide'--------------------------------

calc.abund <- data.frame(rbind(yr_df),
                          est=c(out_rnd_eff$mean$totAdults),
                          q2.5=c(out_rnd_eff$q2.5$totAdults),
                          q97.5=c(out_rnd_eff$q97.5$totAdults),
                          param=rep(c("TotN"), each=nyr))


raw.abund <-count_dat%>%mutate(param="TotN",
                                         est=Estimate_ADULTMF*.64,
                                         q2.5=NA,
                                         q97.5=NA,
                                         yrs=Year,
                                         yr_idx=NA)%>%
  select(colnames(calc.abund))%>%
mutate(param="TotN")

out_rnd_eff$mean$totNMF[19]
out_rnd_eff$q2.5$totNMF[19]
out_rnd_eff$q97.5$totNMF[19]


calc.abund%>%mutate(class="modelled")%>%
  rbind(raw.abund%>%as.data.frame()%>%mutate(class="observed"))%>%
ggplot(aes(x = yrs, y = est, fill=class)) +
    geom_cloud(aes( ymin=q2.5, ymax=q97.5),steps=20, max_alpha = 1,se_mult=1.96)+
  geom_line(aes(color=class)) +
  geom_point(aes(color=class)) +
  theme_ipsum()+
  theme_ipsum()+
  labs(x="Year", y="Abundance", title="Comparing modelled vs raw abundance")+
    geom_vline(xintercept = 2016.5)+
  expand_limits(y=0)




## ----Plot results-summarise vital rates, results='asis'--------------------------------------------------------------------------------------
summary.s <- tribble(
  ~pop,~s, ~s.lower, ~s.upper,
"pre-mgmt",out_rnd_eff$mean$mean_surv_pre, out_rnd_eff$q2.5$mean_surv_pre,out_rnd_eff$q97.5$geom_mean_lambda_pre,
"post-mgmt", out_rnd_eff$mean$mean_surv_post, out_rnd_eff$q2.5$mean_surv_post, out_rnd_eff$q97.5$mean_surv_post)


summary.r <- tribble(
  ~pop,~r, ~r.lower, ~r.upper,
"pre-mgmt",out_rnd_eff$mean$mean_r_pre, out_rnd_eff$q2.5$mean_r_pre, out_rnd_eff$q97.5$mean_r_pre,
"post-mgmt", out_rnd_eff$mean$mean_r_post, out_rnd_eff$q2.5$mean_r_post, out_rnd_eff$q97.5$mean_r_post)

summary.r3 <- tribble(
  ~pop,~r.ad, ~r.ad.lower, ~r.ad.upper,
"pre-mgmt", out_rnd_eff$mean$mean_r3_pre, out_rnd_eff$q2.5$mean_r3_pre, out_rnd_eff$q97.5$mean_r3_pre,
"post-mgmt", out_rnd_eff$mean$mean_r3_post, out_rnd_eff$q2.5$mean_r3_post, out_rnd_eff$q97.5$mean_r3_post)

summary.vr <- summary.s%>%
  left_join(summary.r)%>%
  mutate_if(is.numeric,function(x) round(x,2))


summary.vr <- summary.s%>%
  left_join(summary.r)%>%
  left_join(summary.r3)%>%
  mutate_if(is.numeric,function(x) round(x,2))

summary.vr$Years <- c("2002-2015", "2016-2020")

summary.vr <- summary.vr%>%
  mutate(`s 95% CI`=paste(s.lower,s.upper, sep="-"),
         `r 95% CI`=paste(r.lower,r.upper, sep="-"),
         `r.ad 95% CI`=paste(r.ad.lower,r.ad.upper, sep="-"))%>%
  select(pop, Years, s,`s 95% CI`, r,`r 95% CI`,r.ad,`r.ad 95% CI`)

gt(summary.vr)%>%
  tab_header(
    title = md("Vital Rates")
  ) 


## ----Plot results-summarise growth rates, results='asis'-------------------------------------------------------------------------------------


summary.l <- tribble(
  ~pop,~l, ~l.lower, ~l.upper,
"pre-mgmt",out_rnd_eff$mean$geom_mean_lambda_pre, out_rnd_eff$q2.5$geom_mean_lambda_pre,out_rnd_eff$q97.5$geom_mean_lambda_pre,
"post-mgmt", out_rnd_eff$mean$geom_mean_lambda_post, out_rnd_eff$q2.5$geom_mean_lambda_post, out_rnd_eff$q97.5$geom_mean_lambda_post,
"post-mgmt_iwolf", out_rnd_eff$mean$geom_mean_lambda_post_iwolf, out_rnd_eff$q2.5$geom_mean_lambda_post_iwolf, out_rnd_eff$q97.5$geom_mean_lambda_post_iwolf)%>%
  mutate_if(is.numeric,function(x) round(x,2))


summary.l$Years <- c("2002-2015", "2016-2020", "2017-2020")
colnames(summary.l) <- c("Group", "Lambda", "Lamba.Lower", "Lambda.Upper", "Years")

summary.l <- summary.l%>%
  mutate(`95% CI`=paste(Lamba.Lower,Lambda.Upper, sep="-"))%>%
  dplyr::select(Group, Years, Lambda,`95% CI`)

gt(summary.l)%>%
    tab_header(
    title = md("Population Growth Rates")
  ) 


## ----effects, results='asis'-----------------------------------------------------------------------------------------------------------------

summary.effect <- tribble(
  ~pop,~lambda.dif, ~lower, ~upper,
"pre vs post",out_rnd_eff$mean$diff_geom_mean_lambda_post_to_pre, out_rnd_eff$q2.5$diff_geom_mean_lambda_post_to_pre,out_rnd_eff$q97.5$diff_geom_mean_lambda_post_to_pre,
"pre vs post_iwolf",out_rnd_eff$mean$diff_geom_mean_lambda_post_to_pre_iwolf, out_rnd_eff$q2.5$diff_geom_mean_lambda_post_to_pre_iwolf,out_rnd_eff$q97.5$diff_geom_mean_lambda_post_to_pre_iwolf)%>%
    mutate_if(is.numeric,function(x) round(x,2))



gt(summary.effect)%>%
  tab_header(
    title = md("Treatment Effect")
  ) 

S <- ggs(out_rnd_eff$samples)%>%
  filter(Parameter%in%c("diff_geom_mean_lambda_post_to_pre","diff_geom_mean_lambda_post_to_pre_iwolf"))%>%
  mutate(Parameter=case_when(Parameter%in%"diff_geom_mean_lambda_post_to_pre"~"All Wolf Control (2016-2020)",
                             Parameter%in%"diff_geom_mean_lambda_post_to_pre_iwolf"~"Intense Wolf Control (2017-2020)"
                             ))

ggplot(S, aes(x = value,fill=Parameter)) +
  geom_density(alpha=0.5) +
  theme_ipsum()+
  theme_ipsum()+
  labs(x="Change in Population Growth", y="Posterior Sample Density", title="Treatment Effects on Population Growth")+
    geom_vline(xintercept = 0)

