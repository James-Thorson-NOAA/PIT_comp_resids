#### hacky code from a different experiment, repurposed to
#### explore residuals in composition data
## Cole | 2021-07

### Run and compare models
library(tidyverse)
library(r4ss)
library(snowfall)
library(ggplot2)
theme_set(theme_bw())
packageVersion('r4ss') #  '1.42.0'
source('functions.R')


### ------------------------------------------------------------
### Run bootstrapped data sets back through the model (matches
### perfectly).
## For each boostrap data set, run a retrospective analysis
Nreps <- 498 # messed up starter and put 500 total
reps <- 1:Nreps
Npeels <- 0
peels <- 0:-Npeels
## Setup to run parallel, saving a single core free.
cpus <- parallel::detectCores()-2
sfStop()
sfInit( parallel=TRUE, cpus=cpus)
sfExportAll()
## Run the models
run_model(499:500, model.name='BSAI_FHS', clean.files=FALSE)


### Do the PIT calculations

## read in each simulated data set and the original. Could read
## in from full ss_new but that'd be slower... this still takes a
## few minutes at least
lc <- array(NA, dim=c(85, 48, 499))
## ac <- array(NA, dim=c(1085, 42,499))
for(boot in 0:498){
  wd <- paste0('runs/BSAI_FHS/boot_', boot)
  dat <- SS_readdat(file.path(wd, 'data.ss'), verbose=FALSE)
  ## ac[,,boot+1] <- dat$agecomp[,-(1:9)] %>% as.matrix
  lc[,,boot+1] <- dat$lencomp[,-(1:6)] %>% as.matrix
}
## So first element is the observed data, rest are
## simulated.
## saveRDS(lc, file='results/lencomp_boots.RDS')
lc <- readRDS('results/lencomp_boots.RDS')


## Need to get Pearson residuals, it's long here so pivot wider
replist <- SS_output('runs/BSAI_FHS/boot_0', covar=FALSE,
                     printstats=FALSE, verbose=FALSE)
##saveRDS(replist, 'results/replist.RDS')
replist <- readRDS('results/replist.RDS')

lc.pearson <- replist$lendbase %>%
  select(Yr, FltSvy=Fleet, Gender=Sex,
         Nsamp=Nsamp_in, Pearson, Bin) %>%
  pivot_wider(names_from=Bin, names_prefix='l',
              values_from=Pearson) %>%
  ## Drop negative years and ghost fleets
  filter(FltSvy>0 & Yr>0 & Gender==1) %>%
  arrange(Yr, FltSvy) %>% select(-Gender)

## Now calculate PIT resids
lc.pit <- array(NA, dim=c(dim(lc)[1:2]))
set.seed(1214124)
jit <- function(x) x+runif(length(x),-.5,.5)
for(y in 1:nrow(lc.pit)){ # loop year
  for(b in 1:ncol(lc.pit)){ # loop length bins
    ## P(obs data>simulated data)
    lc.pit[y,b] <- qnorm(mean(jit(lc[y,b,1])>jit(lc[y,b,-1])))
  }
}
## Get back into data.frame and match to Pearson. Note that sex=1
## only and I'm not sure if I did this right
lc.pit <- dat$lencomp %>%
  select(Yr, FltSvy, Nsamp) %>%
  cbind(lc.pit[, 1:length(dat$lbin_vector)]) %>%
  filter(Yr>0 & FltSvy>0) %>%
  arrange(Yr, FltSvy)
names(lc.pit ) <- names(lc.pearson)

## Make sure the two are lined up. There are supposed to be
## differences in Gender and Nsamp, and Seas gets changed too so
## ignoring that. SOmeone should check this carefully
stopifnot(all.equal(lc.pit[,1:2], lc.pearson[,1:2], check.attributes=FALSE))
## Which Nsamp do I want to use here? Probably the one from
## Pearson, need to match or will break pivoting below
lc.pit$Nsamp <- lc.pearson$Nsamp

## Merge together for comparing the two types
lc.all <- bind_rows(cbind(type='PIT', lc.pit),
                    cbind(type='Pearson', lc.pearson)) %>%
  pivot_longer(-c(type, Yr, FltSvy, Nsamp), names_to='bin',
               values_to='residual') %>%
  mutate(FltSvy=paste0('Fleet',FltSvy))

g <- ggplot(lc.all, aes(residual, fill=type)) +
  geom_histogram(position='identity', alpha=.5) + facet_wrap('FltSvy')
ggsave("plots/resids_by_type.png", g, width=7, height=3)

## Scatter plots, cast wide first
lc.all.wide <- lc.all %>%
  pivot_wider(id_cols=c(Yr, FltSvy, Nsamp, bin),
              names_from='type',
              values_from='residual')
g <- ggplot(lc.all.wide, aes(Pearson, PIT)) + geom_point()  +facet_wrap("FltSvy")
ggsave("plots/resids_scatter.png", g, width=7, height=5)

## Look at some bins closely
g <- lc.all.wide%>%  filter(bin %in% c('l6', 'l20', 'l52')) %>%
  ggplot(aes(Pearson, PIT, size=Nsamp, color=FltSvy)) +
  geom_hline(yintercept=0) + geom_vline(xintercept=0)+
  facet_wrap('bin')+ geom_point(alpha=.7)
ggsave("plots/compare_resids_bins.png", g, width=7, height=3)
## Look at some years closely
g <- lc.all.wide%>%  filter(Yr %in% c(1973, 1983, 1999, 2020)) %>%
  ggplot(aes(Pearson, PIT, size=Nsamp, color=FltSvy)) +
  geom_hline(yintercept=0) + geom_vline(xintercept=0)+
  facet_wrap('Yr')+ geom_point(alpha=.7)
ggsave("plots/compare_resids_years.png", g, width=7, height=5)





### Part two is to look at the behavior of Pearson residuals with
### the correctly specified model

## Loop through and read in the comp residuals for all boottrap
## iterations.
lc.list <- list()
for(boot in 0:498){ # 0 is the original file
  ff <- paste0('runs/BSAI_FHS/boot_',boot)
  out <- SS_output(ff, covar=FALSE, printstats=FALSE, verbose=FALSE)
  lc.list[[boot+1]] <- select(out$lendbase, Yr, Fleet, Sex, Bin,
                              Pearson, Obs, Exp, Nsamp_adj) %>% cbind(boot=boot)
}
lcboot <- bind_rows(lc.list)
saveRDS(lcboot, 'results/lcboot.RDS')


## Quick output prep for plotting and analysis
lcboot <- readRDS('results/lcboot.RDS')
lcboot <- mutate(lcboot, Fleet=paste0('fleet',Fleet)) %>%
  ## Drop half the years for easier plotting
  filter(Yr %% 2 == 0)
lcboot0 <- filter(lcboot, boot==0)
lcboot <- lcboot %>% filter(boot!=0)


## Make some big picture plots
g <- ggplot(lcboot, aes(Obs, Pearson, color=Nsamp_adj)) +
  geom_point(alpha=.25) + scale_x_log10() +
  facet_wrap("Fleet", scales='free_x',  ncol=1)
ggsave("plots/null_vs_obs.png", g, width=7, height=5)
g <- lcboot %>%
  ggplot(aes(x=factor(Yr), y=Pearson)) + geom_boxplot() +
  coord_cartesian(ylim=c(-2.5,5))+
  facet_wrap("Fleet", scales='free_x',  ncol=1) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_point(data=lcboot0, mapping=aes(x=factor(Yr), y=Pearson),
             col=2, pch='-')
ggsave("plots/null_distribution_years.png", g, width=7, height=7)
g <- lcboot %>%
  ggplot(aes(x=factor(Bin), y=Pearson)) + geom_boxplot() +
  coord_cartesian(ylim=c(-2.5,5))+ facet_wrap("Fleet", ncol=1) +
  theme(axis.text.x = element_text(angle = 90))+
    geom_point(data=lcboot0, mapping=aes(x=factor(Bin), y=Pearson),
             col=2, pch='-')
ggsave("plots/null_distribution_bins.png", g, width=7, height=7)


## Zoom in a bit? These don't really help I don't think
lcboot2 <- group_by(lcboot, Yr, Fleet, Sex, Bin) %>%
  summarize(mean=mean(Pearson), sd=sd(Pearson), max=max(Pearson),
            range=range(Pearson),
            pvalue=ks.test(Pearson, 'pnorm', mean=0, sd=1)$p.value, .groups='drop')
## ggplot(lcboot2, aes(pvalue)) + geom_histogram()
##   scale_color_gradient(lim=c(0,1))
ggplot(lcboot2, aes(Yr, Bin, color=range)) + geom_point() +
  scale_color_gradient(low='red',  high='blue') +
  facet_wrap('Fleet')


