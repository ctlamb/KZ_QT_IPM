library(dplyr)
library(jagsUI)
l


# Read in raw data

#  Survival
adult_female_survival <- read.csv("C:/Users/saraw/OneDrive/Documents/GitHub/KZ-and-QT-IPMs/data/KZ/adult_female_survival_KZ.csv")

# Recruitment
adult_female_recruit <- read.csv("C:/Users/saraw/OneDrive/Documents/GitHub/KZ-and-QT-IPMs/data/KZ/adult_female_recruit_KZ.csv")

# Sex ratio
adult_sex_ratio <- read.csv("C:/Users/saraw/OneDrive/Documents/GitHub/KZ-and-QT-IPMs/data/KZ/adult_sexratio_KZ.csv")

# Mincounts and abundance estimates and SD's (both adult and calf)
count_dat <- read.csv("C:/Users/saraw/OneDrive/Documents/GitHub/KZ-and-QT-IPMs/data/KZ/count_KZ.csv")

#  Counts in pen (both adult and calf)
adult_female_pencount_in <-  read.csv("C:/Users/saraw/OneDrive/Documents/GitHub/KZ-and-QT-IPMs/data/KZ/adult_female_pencount_in_KZ.csv")
subadult1_female_pencount_in <-  read.csv("C:/Users/saraw/OneDrive/Documents/GitHub/KZ-and-QT-IPMs/data/KZ/sa1_pencount_in_KZ.csv")
calf_female_pencount_recruit <-  read.csv("C:/Users/saraw/OneDrive/Documents/GitHub/KZ-and-QT-IPMs/data/KZ/calf_female_pencount_recruit_KZ.csv")


#  Years of study
yrs <-  seq(from = 1996, to = 2020, by = 1)
nyr <- length(yrs)
yr_idx <- seq(from = 1, to = nyr, by = 1)
yr_df <- as.data.frame(cbind(yrs, yr_idx))

#  Wolf and pen indicator per year
wolfpen_ind <- adult_female_pencount_in$Wolf_control_pen

#  Create vectors of sex ratio
adult_sex_ratio <- adult_sex_ratio$Mean

#  Create vectors of known numbers of adult females and calves going into and coming out  of pen
nPen <- adult_female_pencount_in$Mean
nPenSA <- subadult1_female_pencount_in$Mean
nPenC <- calf_female_pencount_recruit$Mean

#  Recruitment 
#adult_female_recruit <- adult_female_recruit %>%	
	#mutate(Mean=Mean/2, SD=SD/2)

meanr <- array(NA, c(1,1,2)) 
meanr[1,1,1] <- mean(adult_female_recruit$Mean, na.rm = TRUE)
meanr[1,1,2] <- 100


# Survival - need estimates for both penned and control (not penned) population units
means <- array(NA, c(1,2,2))

adult_female_survival_control <- adult_female_survival %>%
	dplyr::filter(Pop == 1) 
means[1,1,1] <- mean(adult_female_survival_control$Mean, na.rm = TRUE)
means[1,1,2] <- 100

adult_female_survival_pen <- adult_female_survival %>%
	dplyr::filter(Pop == 2) 
means[1,2,1] <- mean(adult_female_survival_pen$Mean, na.rm = TRUE)
means[1,2,2] <- 100


# Starting population size of population unit (vector of values for each age class)
n1 <- numeric()
n1[1] <- count_dat$SurveryCount_CALFMF[1]*.5
n1[2] <- (count_dat$SurveryCount_ADULTMF[1]*0.64)*.15  
n1[3] <- (count_dat$SurveryCount_ADULTMF[1]*0.64)*.85  


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
nmin <- nrow(mindat)

calf_mindat <- count_dat %>%
	mutate(dau = 1, 
		age = NA,
		pop = NA,
		mu = SurveryCount_CALFMF,
		tau = NA) %>%
	dplyr::filter(mu > 4) %>%
	left_join(yr_df, by = c("Year" = "yrs")) %>%
	dplyr::select(dau, yr = yr_idx, age, pop, mu, tau) 
nminc <- nrow(calf_mindat)

abundat <- count_dat %>%
	mutate(dau = 1, 
		age = NA,
		pop = NA,
		mu = Estimate_ADULTMF2,
		tau = 1/(SD_Estimate_ADULTMF2 * SD_Estimate_ADULTMF2)) %>%
	dplyr::filter(mu > 4) %>%
	left_join(yr_df, by = c("Year" = "yrs")) %>%
	dplyr::select(dau, yr = yr_idx, age, pop, mu, tau) 
na <- nrow(abundat)

calf_abundat <- count_dat %>%
	mutate(dau = 1, 
		age = NA,
		pop = NA,
		mu = Estimate_CALFMF2,
		tau = 1/(SD_Estimate_CALFMF2 * SD_Estimate_CALFMF2)) %>%
	dplyr::filter(mu > 4) %>%
	left_join(yr_df, by = c("Year" = "yrs")) %>%
	dplyr::select(dau, yr = yr_idx, age, pop, mu, tau) 
nac <- nrow(calf_abundat)



#  Gather data inputs in a list
ipm_dat <- list(nyr = nyr,
	nmin = nmin,
	ns = ns,
	nr = nr,
	na = na,
	nac = nac,
	nminc = nminc,
	#wolfpen_ind = wolfpen_ind,
	means = means,
	meanr = meanr,
	n1 = n1, 
	adult_sex_ratio = adult_sex_ratio,
	mindat = mindat, 
	calf_mindat = calf_mindat, 
	abundat = abundat,
	calf_abundat = calf_abundat,
	nPenC = nPenC,
	nPen = nPen,
	nPenSA = nPenSA,
	sdat = sdat, 
	rdat = rdat)


#  Initial values for N to avoid parent node erros
Nst <- array(10, c(nyr,3,2))
Nst[,,2] <- NA
ipm_inits <- function(){ 
		list(N = Nst)}
	
	
#  Model parameters to monitor
model_parms <- c("lambda","logla",
	"totN", "N", "S", "R", "r1", #"r2", 
	#"wolfpen_eff_r", "wolfpen_eff_s",
	"sa_yr", "r_yr", 
	"totCalves", "totAdults",	
	"totAdultsMF", "totCalvesMF", "totNMF") #,
	#"geom_mean_lambda",
	#"geom_mean_lambda_pre", "geom_mean_lambda_post",
	#"diff_geom_mean_lambda_post_to_pre","diff_geom_mean_lambda_post_to_pre_iwolf",
	#"geom_mean_lambda_post_iwolf",
	#"mean_surv_pre", "mean_surv_post",
	#"mean_r_pre", "mean_r_post")


nt <- 3
nb <- 8000
nc <- 3
nad <- 5000
ni <- 50000

out_rnd_eff <- jagsUI::jags(ipm_dat, 
	inits = ipm_inits,
	model_parms,
	model.file = "C:/Users/saraw/OneDrive/Documents/GitHub/KZ-and-QT-IPMs/jags/KZ_IPM_rnd_eff_ABUND.txt",
	n.chains = nc, 
	n.iter = ni,
	n.burnin = nb,
	n.thin = nt,
	n.adapt = nad)


mcmcplots::mcmcplot(out_rnd_eff$samples, parms = c("r1"))


mcmcplots::mcmcplot(out_rnd_eff$samples, parms = c("totAdultsMF"))
mcmcplots::mcmcplot(out_rnd_eff$samples, parms = c("totCalvesMF"))



##ABUNDANCE
res_df <- data.frame(yr_df,
                     est=out_rnd_eff$mean$totAdultsMF,
                     q2.5=out_rnd_eff$q2.5$totAdultsMF, 
                     q97.5=out_rnd_eff$q97.5$totAdultsMF, 
                     param=rep("Model_Est_Adults_MF", 25))
tmp <- count_dat %>% 
	cbind(yr_idx = yr_df$yr_idx) %>% 
	dplyr::select(Estimate_ADULTMF, Estimate_ADULTMF2)
res_df <- res_df %>% 
	cbind(tmp)

ggplot(res_df,aes(x = yrs, y = est, ymin=q2.5, ymax=q97.5)) +
	geom_cloud(steps=20, max_alpha = 1,se_mult=1.96, fill = "#EBCC2A")+
	geom_line(size = 1) +
	geom_point(aes(x = yrs, y = Estimate_ADULTMF), shape = 17, col = "#3B9AB2", size = 4) +
	#geom_point(aes(x = yrs, y = Estimate_ADULTMF2), shape = 19, col = "#F21A00", size = 4) +
	#theme_ipsum()+
	theme_bw()+
	ylab("Abdundance estimate")+
	xlab("Year") 
