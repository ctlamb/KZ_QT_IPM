KZ and QT IPM results
================
Sara Williams, Hans Martin, and Clayton Lamb
06 October, 2020

\#See folders KZ and QT for the IPM’s for each herd

\#\#Load Data

``` r
library(ggmcmc)
library(jagsUI)
library(knitr)
library(gt)
library(ggallin)
library(hrbrthemes)
library(tidyverse)
library(rjags)
library(MCMCvis)
library(ggpubr)
library(readxl)
library(piggyback)
library(here)

pb_download("kz_out_rnd_eff.Rdata", 
            repo = "ctlamb/KZ_QT_IPM",
            dest = here::here())
load(here::here("kz_out_rnd_eff.Rdata"))
kz <- out_rnd_eff
unlink(here::here("kz_out_rnd_eff.Rdata"))


yrs <-  seq(from = 1995, to = 2020, by = 1)
nyr <- length(yrs)
yr_idx <- seq(from = 1, to = nyr, by = 1)
kz_yr_df <- as.data.frame(cbind(yrs, yr_idx))



pb_download("qt_out_rnd_eff.Rdata", 
            repo = "ctlamb/KZ_QT_IPM",
            dest = here::here())
load(here::here("qt_out_rnd_eff.Rdata"))
qt <- out_rnd_eff
unlink(here::here("qt_out_rnd_eff.Rdata"))


yrs <-  seq(from = 2001, to = 2020, by = 1)
nyr <- length(yrs)
yr_idx <- seq(from = 1, to = nyr, by = 1)
qt_yr_df <- as.data.frame(cbind(yrs, yr_idx))



##counts
qt_count_dat <- read_xlsx(here::here("data", "QT", "Count_summary_QT.xlsx"))
kz_count_dat <- read_xlsx(here::here("data", "KZ", "Count_summary_KZ.xlsx"))

##raw vital rates
qt_vr <- read_csv(here::here("data", "QT", "vitalrate_validation_QT.csv"))
kz_vr <- read_csv(here::here("data", "KZ", "vitalrate_validation_KZ.csv"))
```

\#\#ABUNDANCE

``` r
res_df <- data.frame(kz_yr_df,
                     est=c(kz$mean$totNMF),
                     q2.5=c(kz$q2.5$totNMF),
                     q97.5=c(kz$q97.5$totNMF),
                     param="Total M+F",
                     herd="Klinse-Za")%>%
  rbind(
    data.frame(qt_yr_df,
                     est=c(qt$mean$totNMF),
                     q2.5=c(qt$q2.5$totNMF),
                     q97.5=c(qt$q97.5$totNMF),
                     param="Total M+F",
                     herd="Quintette")
    
  )

#res_df <- res_df[-(nyr+1),]

ggplot(res_df,aes(x = yrs, y = est, ymin=q2.5, ymax=q97.5, fill=herd)) +
  #geom_cloud(steps=50, max_alpha = 1,se_mult=1.96)+
  geom_ribbon(alpha=0.4)+
  geom_line() +
  geom_point() +
  theme_ipsum()+
  theme(legend.position = "none")+
  ylab("Abundance")+
  xlab("Year")+
  facet_wrap(vars(herd), scales="free_y")+
  labs(x="Year",title="Population Trajectory")+
    expand_limits(y=0)+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text = element_text(size=10),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15))+
    geom_vline(data=data.frame(herd=unique(res_df$herd), year=c(2012.5,2015.5)),aes(xintercept = year),linetype="dashed")
```

![](README_files/figure-gfm/Plot%20results-abundance-1.png)<!-- -->

``` r
ggsave(here::here("plots", "abundance_MF.png"), width=8, height=5)
```

\#\#ABUNDANCE FIT

``` r
fit_df <- res_df%>%
  mutate(type="modelled")%>%
  rbind(data.frame(kz_yr_df,
                     est=kz_count_dat$Estimate_ADULTMF + kz_count_dat$Estimate_CALFMF,
                     q2.5=NA,
                     q97.5=NA,
                     param="Total M+F",
                     herd="Klinse-Za",
                     type="Observed_Estimate")%>%
  rbind(
     data.frame(qt_yr_df,
                     est=qt_count_dat$Estimate_ADULTMF + qt_count_dat$Estimate_CALFMF,
                     q2.5=NA,
                     q97.5=NA,
                     param="Total M+F",
                     herd="Quintette",
                     type="Observed_Estimate")
    
  )%>%
    rbind(data.frame(kz_yr_df,
                     est=kz_count_dat$`KZ total count_MAX`,
                     q2.5=NA,
                     q97.5=NA,
                     param="Total M+F",
                     herd="Klinse-Za",
                     type="Observed_Mincount"))%>%
  rbind(
     data.frame(qt_yr_df,
                     est=qt_count_dat$MinCount_ADULTMF + qt_count_dat$MinCount_CALFMF, 
                     q2.5=NA,
                     q97.5=NA,
                     param="Total M+F",
                     herd="Quintette",
                     type="Observed_Mincount")
    
  )
)

#res_df <- res_df[-(nyr+1),]

ggplot(data=fit_df%>%filter(type%in%"modelled"), aes(x = yrs, y = est, ymin=q2.5, ymax=q97.5, fill=type)) +
  #geom_cloud(steps=20, max_alpha = 1,se_mult=1.96)+
  geom_ribbon(alpha=0.4)+
  geom_line(aes(color=type)) +
  geom_point(aes(color=type)) +
  geom_jitter(data=fit_df%>%filter(!type%in%"modelled"),aes(color=type), size=1) +
  theme_ipsum()+
  ylab("Abundance")+
  xlab("Year")+
  facet_wrap(vars(herd), scales="free_y")+
  labs(x="Year",title="Population Trajectory")+
    expand_limits(y=0)+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text = element_text(size=10),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15))
```

![](README_files/figure-gfm/Plot%20results-abundance%20fit-1.png)<!-- -->

``` r
ggsave(here::here("plots", "abundancefit_MF.png"), width=9, height=5)
```

\#\#VITAL RATES

``` r
#R
pop_df_r <- rbind(
                  data.frame(estimate=c(
                                  kz$mean$totCalvesP/kz$mean$totAdultsP),
                       lower=c(
                                  kz$q2.5$totCalvesP/kz$q2.5$totAdultsP),
                       upper=c(
                                  kz$q97.5$totCalvesP/kz$q97.5$totAdultsP),
                       pop=rep(c("KZ-Wolf + Pen"), each=length(kz$mean$R)),
                       param="Recruitment",
                       yrs=rep(kz_yr_df$yrs, times=1)),
                  
                    data.frame(estimate=c(
                                  kz$mean$R[,1]),
                       lower=c(
                                  kz$q2.5$R[,1]),
                       upper=c(
                                  kz$q97.5$R[,1]),
                       pop=rep(c("KZ-Wolf"), each=length(kz$mean$R)),
                       param="Recruitment",
                       yrs=rep(kz_yr_df$yrs, times=1)),
                  
                  data.frame(estimate=qt$mean$R,
                       lower=qt$q2.5$R,
                       upper=qt$q97.5$R,
                       pop=rep(c("QT-Wolf"), each=length(qt$mean$R)),
                       param="Recruitment",
                       yrs=rep(qt_yr_df$yrs, times=1)))%>%
  mutate(lower=case_when(lower<0~0,
                         TRUE~lower))

#S
pop_df_s <- data.frame(estimate=c(kz$mean$S[,1],
                                  kz$mean$S[,2]),
                       lower=c(kz$q2.5$S[,1],
                               kz$q2.5$S[,2]),
                       upper=c(kz$q97.5$S[,1],
                               kz$q97.5$S[,2]),
                       pop=rep(c("KZ-Wolf", "KZ-Wolf + Pen"), each=length(kz$mean$S)/2),
                       param="Survival",
                       yrs=rep(kz_yr_df$yrs, times=2))%>%
  rbind(
    data.frame(estimate=c(qt$mean$S),
                       lower=c(qt$q2.5$S),
                       upper=c(qt$q97.5$S),
                       pop="QT-Wolf",
                       param="Survival",
                       yrs=qt_yr_df$yrs)
  )



mod_vr <- rbind(pop_df_r, pop_df_s)%>%
         filter(yrs>1995)%>%
         mutate(estimate=case_when(param=="Recruitment" & yrs<=2015 & pop=="KZ-Wolf + Pen"~NA_real_, 
                                   param=="Survival" & yrs<=2014 & pop=="KZ-Wolf + Pen"~NA_real_, 
                                   TRUE~estimate))

ggplot(mod_vr,
       aes(x = yrs, y = estimate, fill=pop,ymin=lower, ymax=upper)) +
  geom_line(aes(color=pop)) +
  geom_ribbon(alpha=0.4)+
  theme_ipsum()+
  labs(x="Year", y="Rate", title="Annual vital rates")+
    geom_vline(data=data.frame(pop=unique(mod_vr$pop), yrs=c(2013.5,2012.5,2015.5), param="Recruitment")%>%
                 rbind(data.frame(pop=unique(mod_vr$pop), yrs=c(2013.5,2012.5,2015.5), param="Survival")),aes(xintercept = yrs),linetype="dashed")+
    facet_wrap(vars(param, pop), scales="free_y")
```

![](README_files/figure-gfm/Plot%20results-vital%20rates-1.png)<!-- -->

\#\#COMPARE TO RAW VITAL RATE DATA

``` r
kz_vr%>%mutate(pop=case_when(pop%in%"Free"~"KZ-Wolf",
                             pop%in%"Pen"~"KZ-Wolf + Pen"))%>%
  rbind(qt_vr%>%mutate(pop=case_when(pop%in%"Free"~"QT-Wolf")))%>%
  mutate(type="Observed")%>%
  rbind(mod_vr%>%mutate(type="Modelled"))%>%
  ggplot(aes(x = yrs, y = estimate, fill=type, ymin=lower, ymax=upper)) +
  geom_line(aes(color=type)) +
  geom_ribbon(alpha=0.4)+
  theme_ipsum()+
  labs(x="Year", y="Rate", title="Annual vital rates")+
    geom_vline(data=data.frame(pop=unique(mod_vr$pop), yrs=c(2013.5,2012.5,2015.5), param="Recruitment")%>%
                 rbind(data.frame(pop=unique(mod_vr$pop), yrs=c(2013.5,2012.5,2015.5), param="Survival")),aes(xintercept = yrs),linetype="dashed")+
    facet_wrap(vars(param, pop), scales="free_y")
```

![](README_files/figure-gfm/Plot%20results-vital%20rates%20vs%20raw%20data-1.png)<!-- -->

``` r
ggsave(here::here("plots", "vitalratefit_F.png"), width=9, height=5)
```

\#\#WOLF EFFECT

``` r
##Intense Wolf Control== (2017-2020)



S <- ggs(qt$samples)%>%
  filter(Parameter%in%c("diff_geom_mean_lambda_post_to_pre","diff_geom_mean_lambda_post_to_pre_iwolf"))%>%
  mutate(Parameter=case_when(Parameter%in%"diff_geom_mean_lambda_post_to_pre"~"All Wolf Control",
                             Parameter%in%"diff_geom_mean_lambda_post_to_pre_iwolf"~"Intense Wolf Control"
                             ),
         Herd="Quintette")


S2 <- ggs(kz$samples)%>%
  filter(Parameter%in%c("wolf_eff","wolf_eff_iwolf"))%>%
    mutate(Parameter=case_when(Parameter%in%"wolf_eff"~"All Wolf Control",
                             Parameter%in%"wolf_eff_iwolf"~"Intense Wolf Control"
                             ),
         Herd="Klinse-Za")

S3 <- rbind(S,S2)

ggplot(S3, aes(x = value,fill=Parameter)) +
  geom_density(alpha=0.5) +
  theme_ipsum()+
  theme_ipsum()+
  labs(x="Change in Population Growth", y="Posterior Samples", title="Effect of Wolf Control on Population Growth")+
    geom_vline(xintercept = 0, linetype="dashed")+
  facet_wrap(vars(Herd), ncol=2)+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text = element_text(size=10),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15))
```

![](README_files/figure-gfm/wolf-1.png)<!-- -->

``` r
ggsave(here::here("plots", "wolf_effect.png"), width=9, height=5)
```

\#\#KZ EFFECT

``` r
S <- ggs(kz$samples)%>%
  filter(Parameter%in%c("wolf_eff","pen_eff", "wolf_eff_iwolf","pen_eff_iwolf"))

a <- S%>%
  pivot_wider(id_cols=c("Iteration", "Chain"), names_from = "Parameter", values_from = "value")%>%
  mutate(pen_eff=pen_eff+wolf_eff,
         pen_eff_iwolf=pen_eff_iwolf+wolf_eff_iwolf)%>%
  pivot_longer(-c("Iteration", "Chain"))%>%
  filter(name%in%c("pen_eff", "wolf_eff"))%>%
  mutate(name=case_when(name%in%"pen_eff"~"Pen + Wolf",
                        name%in%"wolf_eff"~"Wolf"))%>%
  rename(Treatment=name)%>%
ggplot(aes(x = value,fill=Treatment)) +
  geom_density(alpha=0.5) +
  theme_ipsum()+
  theme_ipsum()+
  labs(x="Change in Population Growth", y="Posterior Samples", title="Klinse-Za Treatment Effects")+
    geom_vline(xintercept = 0, linetype="dashed")+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text = element_text(size=10),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15))


pop.sim <-data.frame(yrs=rep(c(2014:2020), times=3),
                     est=c(kz$mean$simTotC, kz$mean$simTotP, kz$mean$simTotBAU),
                     q2.5=c(kz$q2.5$simTotC, kz$q2.5$simTotP, kz$q2.5$simTotBAU),
                     q97.5=c(kz$q97.5$simTotC, kz$q97.5$simTotP, kz$q97.5$simTotBAU),
                     Treatment=rep(c("Wolf", "Pen + Wolf", "Control"), each=length(kz$q97.5$simTotC)))



b <- ggplot(pop.sim%>%mutate(Treatment=fct_relevel(Treatment,"Pen + Wolf","Wolf","Control")),
       aes(x = yrs, y = est, ymin=q2.5, ymax=q97.5, fill=Treatment, linetype=Treatment)) +
  geom_cloud(steps=20, max_alpha = 0.8,se_mult=1.96)+
  #geom_ribbon(alpha=0.2)+
  geom_line() +
  theme_ipsum()+
  labs(x="Year", y="Abundance", title="Simulated Female Abundance", subtitle="Same starting abundance, treatment applied to all females")+
  expand_limits(y=0)+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text = element_text(size=10),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15))




ggarrange(a,b,labels="AUTO")
```

![](README_files/figure-gfm/kz%20effect-1.png)<!-- -->

``` r
ggsave(here::here("plots", "KZ_effect_and_sim.png"), width=11, height=5)
```

\#\#Wolf vs Pen effect

``` r
S <- ggs(kz$samples)%>%
  filter(Parameter%in%c("wolf_eff_proportion","pen_eff_proportion"))

S%>%
    mutate(Treatment=case_when(Parameter%in%"pen_eff_proportion"~"Pen",
                        Parameter%in%"wolf_eff_proportion"~"Wolf"))%>%
ggplot(aes(x = value,fill=Treatment)) +
  geom_density(alpha=0.5) +
  theme_ipsum()+
  labs(x="Proportion", y="Posterior Samples", title="Proportion Increase Attributable to Treatment")+
  xlim(0,1)+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text = element_text(size=10),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15))
```

![](README_files/figure-gfm/Wolf%20vs%20Pen%20effect-1.png)<!-- -->

``` r
ggsave(here::here("plots", "KZ_effect_prop.png"), width=7, height=5)
```

\#\#Summarize population growth

``` r
summary.l <- tribble(
  ~pop,~period,~l, ~l.lower, ~l.upper,
"Klinse-Za","pre-mgmt",kz$mean$geom_mean_lambda_prepen, kz$q2.5$geom_mean_lambda_prepen,kz$q97.5$geom_mean_lambda_prepen,
"Klinse-Za","post-mgmt", kz$mean$geom_mean_lambda_postpen, kz$q2.5$geom_mean_lambda_postpen, kz$q97.5$geom_mean_lambda_postpen,
"Quintette","pre-mgmt",qt$mean$geom_mean_lambda_pre, qt$q2.5$geom_mean_lambda_pre,qt$q97.5$geom_mean_lambda_pre,
"Quintette","post-mgmt", qt$mean$geom_mean_lambda_post, qt$q2.5$geom_mean_lambda_post, qt$q97.5$geom_mean_lambda_post,
"Quintette","post-mgmt_iwolf", qt$mean$geom_mean_lambda_post_iwolf, qt$q2.5$geom_mean_lambda_post_iwolf, qt$q97.5$geom_mean_lambda_post_iwolf
)%>%
  mutate_if(is.numeric,function(x) round(x,2))


summary.l$Years <- c("1996-2012", "2014-2020","2002-2015", "2016-2020", "2017-2020")
colnames(summary.l) <- c("Herd", "Period", "Lambda", "Lamba.Lower", "Lambda.Upper", "Years")


summary.l <- summary.l%>%
  mutate(`95% CI`=paste(Lamba.Lower,Lambda.Upper, sep="-"))%>%
  select(Herd, Period, Years, Lambda,`95% CI`)

kable(summary.l, caption="Population Growth Rates")
```

| Herd      | Period           | Years     | Lambda | 95% CI    |
| :-------- | :--------------- | :-------- | -----: | :-------- |
| Klinse-Za | pre-mgmt         | 1996-2012 |   0.89 | 0.88-0.89 |
| Klinse-Za | post-mgmt        | 2014-2020 |   1.07 | 1.04-1.09 |
| Quintette | pre-mgmt         | 2002-2015 |   0.94 | 0.9-0.97  |
| Quintette | post-mgmt        | 2016-2020 |   1.01 | 0.92-1.1  |
| Quintette | post-mgmt\_iwolf | 2017-2020 |   1.08 | 0.99-1.16 |

Population Growth Rates

\#\#Summarize effect of treatments

``` r
summary.effect <- tribble(
  ~pop,~period,~lambda.dif, ~lower, ~upper,
"Klinse-Za","wolf + pen",kz$mean$diff_geom_mean_lambda_post_to_pre, kz$q2.5$diff_geom_mean_lambda_post_to_pre,kz$q97.5$diff_geom_mean_lambda_post_to_pre,
"Klinse-Za","wolf",kz$mean$wolf_eff, kz$q2.5$wolf_eff,kz$q97.5$wolf_eff,
"Klinse-Za","pen",kz$mean$pen_eff, kz$q2.5$pen_eff,kz$q97.5$pen_eff,
"Quintette","wolf",qt$mean$diff_geom_mean_lambda_post_to_pre, qt$q2.5$diff_geom_mean_lambda_post_to_pre,qt$q97.5$diff_geom_mean_lambda_post_to_pre,
"Quintette","iwolf",qt$mean$diff_geom_mean_lambda_post_to_pre_iwolf, qt$q2.5$diff_geom_mean_lambda_post_to_pre_iwolf,qt$q97.5$diff_geom_mean_lambda_post_to_pre_iwolf)%>%
    mutate_if(is.numeric,function(x) round(x,3))



gt(summary.effect)%>%
  tab_header(
    title = md("Treatment Effect")
  ) 
```

<!--html_preserve-->

<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#tgeshwfedm .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#tgeshwfedm .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#tgeshwfedm .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#tgeshwfedm .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 4px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#tgeshwfedm .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#tgeshwfedm .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#tgeshwfedm .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#tgeshwfedm .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#tgeshwfedm .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#tgeshwfedm .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#tgeshwfedm .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#tgeshwfedm .gt_group_heading {
  padding: 8px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
}

#tgeshwfedm .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#tgeshwfedm .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#tgeshwfedm .gt_from_md > :first-child {
  margin-top: 0;
}

#tgeshwfedm .gt_from_md > :last-child {
  margin-bottom: 0;
}

#tgeshwfedm .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#tgeshwfedm .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 12px;
}

#tgeshwfedm .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#tgeshwfedm .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#tgeshwfedm .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#tgeshwfedm .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#tgeshwfedm .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#tgeshwfedm .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#tgeshwfedm .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#tgeshwfedm .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#tgeshwfedm .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#tgeshwfedm .gt_left {
  text-align: left;
}

#tgeshwfedm .gt_center {
  text-align: center;
}

#tgeshwfedm .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#tgeshwfedm .gt_font_normal {
  font-weight: normal;
}

#tgeshwfedm .gt_font_bold {
  font-weight: bold;
}

#tgeshwfedm .gt_font_italic {
  font-style: italic;
}

#tgeshwfedm .gt_super {
  font-size: 65%;
}

#tgeshwfedm .gt_footnote_marks {
  font-style: italic;
  font-size: 65%;
}
</style>

<div id="tgeshwfedm" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">

<table class="gt_table">

<thead class="gt_header">

<tr>

<th colspan="5" class="gt_heading gt_title gt_font_normal" style>

Treatment Effect

</th>

</tr>

<tr>

<th colspan="5" class="gt_heading gt_subtitle gt_font_normal gt_bottom_border" style>

</th>

</tr>

</thead>

<thead class="gt_col_headings">

<tr>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

pop

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

period

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">

lambda.dif

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">

lower

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">

upper

</th>

</tr>

</thead>

<tbody class="gt_table_body">

<tr>

<td class="gt_row gt_left">

Klinse-Za

</td>

<td class="gt_row gt_left">

wolf + pen

</td>

<td class="gt_row gt_right">

0.180

</td>

<td class="gt_row gt_right">

0.153

</td>

<td class="gt_row gt_right">

0.207

</td>

</tr>

<tr>

<td class="gt_row gt_left">

Klinse-Za

</td>

<td class="gt_row gt_left">

wolf

</td>

<td class="gt_row gt_right">

0.095

</td>

<td class="gt_row gt_right">

0.027

</td>

<td class="gt_row gt_right">

0.160

</td>

</tr>

<tr>

<td class="gt_row gt_left">

Klinse-Za

</td>

<td class="gt_row gt_left">

pen

</td>

<td class="gt_row gt_right">

0.085

</td>

<td class="gt_row gt_right">

0.026

</td>

<td class="gt_row gt_right">

0.149

</td>

</tr>

<tr>

<td class="gt_row gt_left">

Quintette

</td>

<td class="gt_row gt_left">

wolf

</td>

<td class="gt_row gt_right">

0.071

</td>

<td class="gt_row gt_right">

\-0.038

</td>

<td class="gt_row gt_right">

0.179

</td>

</tr>

<tr>

<td class="gt_row gt_left">

Quintette

</td>

<td class="gt_row gt_left">

iwolf

</td>

<td class="gt_row gt_right">

0.138

</td>

<td class="gt_row gt_right">

0.047

</td>

<td class="gt_row gt_right">

0.226

</td>

</tr>

</tbody>

</table>

</div>

<!--/html_preserve-->

\#\#Summarize vital rates

``` r
##define fn
gm_mean = function(a){prod(a)^(1/length(a))}

summary.s <- tribble(
  ~pop, ~period, ~s, ~s.lower, ~s.upper,
"Quintette","pre-mgmt",qt$mean$mean_surv_pre, qt$q2.5$mean_surv_pre, qt$q97.5$geom_mean_lambda_pre,
"Quintette","post-mgmt", qt$mean$mean_surv_post, qt$q2.5$mean_surv_post, qt$q97.5$mean_surv_post,
"Klinse-Za","pre-mgmt",kz$mean$mean_surv_pre, kz$q2.5$mean_surv_pre, kz$q97.5$geom_mean_lambda_pre,
"Klinse-Za","post-mgmt (wolf+pen)", kz$mean$mean_surv_pen, kz$q2.5$mean_surv_pen, kz$q97.5$mean_surv_pen,
"Klinse-Za","post-mgmt (wolf)", kz$mean$mean_surv_wolf, kz$q2.5$mean_surv_wolf, kz$q97.5$mean_surv_wolf)


summary.r <- tribble(
  ~pop,~period,~r, ~r.lower, ~r.upper,
"Quintette","pre-mgmt",gm_mean(qt$mean$R[2:16]), gm_mean(qt$q2.5$R[2:16]),gm_mean(qt$q97.5$R[2:16]),
"Quintette","post-mgmt", qt$mean$mean_r_post, qt$q2.5$mean_r_post, qt$q97.5$mean_r_post,
"Klinse-Za","pre-mgmt",gm_mean(kz$mean$R[2:19,1]), gm_mean(kz$q2.5$R[2:19,1]),gm_mean(kz$q97.5$R[2:19,1]),
"Klinse-Za","post-mgmt (wolf+pen)", kz$mean$mean_r_pen, kz$q2.5$mean_r_pen, kz$q97.5$mean_r_pen,
"Klinse-Za","post-mgmt (wolf)",kz$mean$mean_r_wolf, kz$q2.5$mean_r_wolf, kz$q97.5$mean_r_wolf)

summary.vr <- summary.s%>%
  left_join(summary.r)%>%
  mutate_if(is.numeric,function(x) round(x,2))
```

    ## Joining, by = c("pop", "period")

``` r
summary.vr$Years <- c("2002-2015", "2016-2020", "1995-2012", "2014-2020", "2013-2020")

summary.vr <- summary.vr%>%
  mutate(`s 95% CI`=paste(s.lower,s.upper, sep="-"),
         `r 95% CI`=paste(r.lower,r.upper, sep="-"))%>%
  dplyr::select(pop, period,Years, s,`s 95% CI`, r,`r 95% CI`)



colnames(summary.vr) <- c("Group", "Period", "Years", "AF Survival","95% CI", "Recruitment","r95% CI")


gt(summary.vr)%>%
  tab_header(
    title = md("Vital Rates")
  ) 
```

<!--html_preserve-->

<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#ohzkmxaroi .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#ohzkmxaroi .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#ohzkmxaroi .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#ohzkmxaroi .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 4px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#ohzkmxaroi .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ohzkmxaroi .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#ohzkmxaroi .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#ohzkmxaroi .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#ohzkmxaroi .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#ohzkmxaroi .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#ohzkmxaroi .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#ohzkmxaroi .gt_group_heading {
  padding: 8px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
}

#ohzkmxaroi .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#ohzkmxaroi .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#ohzkmxaroi .gt_from_md > :first-child {
  margin-top: 0;
}

#ohzkmxaroi .gt_from_md > :last-child {
  margin-bottom: 0;
}

#ohzkmxaroi .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#ohzkmxaroi .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 12px;
}

#ohzkmxaroi .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#ohzkmxaroi .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#ohzkmxaroi .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#ohzkmxaroi .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#ohzkmxaroi .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ohzkmxaroi .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#ohzkmxaroi .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#ohzkmxaroi .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#ohzkmxaroi .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#ohzkmxaroi .gt_left {
  text-align: left;
}

#ohzkmxaroi .gt_center {
  text-align: center;
}

#ohzkmxaroi .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#ohzkmxaroi .gt_font_normal {
  font-weight: normal;
}

#ohzkmxaroi .gt_font_bold {
  font-weight: bold;
}

#ohzkmxaroi .gt_font_italic {
  font-style: italic;
}

#ohzkmxaroi .gt_super {
  font-size: 65%;
}

#ohzkmxaroi .gt_footnote_marks {
  font-style: italic;
  font-size: 65%;
}
</style>

<div id="ohzkmxaroi" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">

<table class="gt_table">

<thead class="gt_header">

<tr>

<th colspan="7" class="gt_heading gt_title gt_font_normal" style>

Vital Rates

</th>

</tr>

<tr>

<th colspan="7" class="gt_heading gt_subtitle gt_font_normal gt_bottom_border" style>

</th>

</tr>

</thead>

<thead class="gt_col_headings">

<tr>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

Group

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

Period

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

Years

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">

AF Survival

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

95% CI

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">

Recruitment

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

r95% CI

</th>

</tr>

</thead>

<tbody class="gt_table_body">

<tr>

<td class="gt_row gt_left">

Quintette

</td>

<td class="gt_row gt_left">

pre-mgmt

</td>

<td class="gt_row gt_left">

2002-2015

</td>

<td class="gt_row gt_right">

0.84

</td>

<td class="gt_row gt_left">

0.81-0.97

</td>

<td class="gt_row gt_right">

0.12

</td>

<td class="gt_row gt_left">

0.07-0.17

</td>

</tr>

<tr>

<td class="gt_row gt_left">

Quintette

</td>

<td class="gt_row gt_left">

post-mgmt

</td>

<td class="gt_row gt_left">

2016-2020

</td>

<td class="gt_row gt_right">

0.90

</td>

<td class="gt_row gt_left">

0.86-0.94

</td>

<td class="gt_row gt_right">

0.18

</td>

<td class="gt_row gt_left">

0.14-0.21

</td>

</tr>

<tr>

<td class="gt_row gt_left">

Klinse-Za

</td>

<td class="gt_row gt_left">

pre-mgmt

</td>

<td class="gt_row gt_left">

1995-2012

</td>

<td class="gt_row gt_right">

0.74

</td>

<td class="gt_row gt_left">

0.71-0.89

</td>

<td class="gt_row gt_right">

0.21

</td>

<td class="gt_row gt_left">

0.12-0.32

</td>

</tr>

<tr>

<td class="gt_row gt_left">

Klinse-Za

</td>

<td class="gt_row gt_left">

post-mgmt (wolf+pen)

</td>

<td class="gt_row gt_left">

2014-2020

</td>

<td class="gt_row gt_right">

0.91

</td>

<td class="gt_row gt_left">

0.91-0.91

</td>

<td class="gt_row gt_right">

0.24

</td>

<td class="gt_row gt_left">

0.24-0.24

</td>

</tr>

<tr>

<td class="gt_row gt_left">

Klinse-Za

</td>

<td class="gt_row gt_left">

post-mgmt (wolf)

</td>

<td class="gt_row gt_left">

2013-2020

</td>

<td class="gt_row gt_right">

0.85

</td>

<td class="gt_row gt_left">

0.79-0.91

</td>

<td class="gt_row gt_right">

0.21

</td>

<td class="gt_row gt_left">

0.17-0.25

</td>

</tr>

</tbody>

</table>

</div>

<!--/html_preserve-->
