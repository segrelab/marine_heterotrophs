### import and organize growth profile data ###
# This file formats and combines the growth data from four experiments whose results were combined in this project (GE2-5).
# load required packages
library(readxl)
library(here)
library(tidyverse)

## organize growth data ##

# read in data
ge2 <- read_excel(here('data', 'ge2_all copy.xlsx')) %>% mutate(exp = 'ge2')
ge3 <- read_excel(here('data', 'ge3_all copy.xlsx')) %>% mutate(exp = 'ge3')
ge4 <- read_excel(here('data', 'ge4_all copy.xlsx')) %>% mutate(exp = 'ge4')
ge5 <- read_excel(here('data', 'ge5_all copy.xlsx')) %>% mutate(exp = 'ge5')
tax <- read_excel(here('data', 'strain_tax copy.xlsx'))

# combine data, make variables consistent
allGE <- bind_rows(ge2, ge3, ge4, ge5) %>% 
  dplyr::select(sample, plate_id, well, tmpt, OD600, medium, strain, exp) %>%
  mutate(strain = replace(strain, strain %in% c('none', 'blank'), 'negative')) %>% #rename neg ctrls, remake grp to include negative
  unite('grp', medium, strain, sep = '.', remove = FALSE) %>%
  unite('grp2', exp, grp, sep = '.', remove = FALSE)
rm(ge2, ge3, ge4, ge5)

# add taxonomy info to allGE
allGE <- allGE %>% inner_join(tax, by = 'strain') 
# make the media factors to order data in plots
allGE$medium <- factor(allGE$medium, 
                       levels = c('difcoMB', 'HMBcmpt', 'HMBpep', 'HMBaa', 'HMBlips', 'HMBoligo', 'HMBorg', 'HMBntrl', 'HMBamisug', 'HMBacdsug', 'HMB--'))

# clean up the data in allGE
cleanGE <- allGE %>%
  # remove PR1white and PR1red, did not finish resequencing to identify
  filter(!strain %in% c('PR1red', 'PR1white')) %>%
  # I made a mistake loading these samples, see lab notebook
  filter(!(strain %in% c('atlcs', 'fors') & exp == 'ge3')) %>%
  # remove aberrant time points, OD spikes/valleys likely to be technical artifacts 
  filter(!(grp %in% c('HMBcmpt.bdim', 'HMB--.bdim') & tmpt == 162)) %>%
  filter(!(grp == 'HMBacdsug.fang' & tmpt == 189)) %>%
  filter(!(grp %in% c('HMBntrl.kuro', 'HMBamisug.kuro', 'HMBacdsug.kuro') & tmpt == 189))

# new list of strains to loop through
clean_str <- unique(cleanGE$strain) %>% sort()
# don't include no-cell control in loop
clean_str <- clean_str[!clean_str == 'negative']

# find maxOD
# find the average maxOD across technical replicates, and the time point at which it's reached
maxOD <- cleanGE %>% 
  filter(tmpt < 265) %>% # only 2 experiments went longer than 264 hours
  group_by(grp2, tmpt) %>% mutate(avgMaxOD = mean(OD600)) %>% dplyr::select(grp2, tmpt, avgMaxOD) %>% 
  ungroup() %>% group_by(grp2) %>% filter(avgMaxOD == max(avgMaxOD)) %>% distinct() %>% rename(t_maxOD = tmpt)

# the maxOD df has the first time point at which the average max OD was reached for every strain
maxOD <- maxOD %>% group_by(grp2) %>% filter(t_maxOD == min(t_maxOD)) 

# add the t_maxOD and avgMaxOD to each row of cleanGE
cleanGEwMaxOD <- cleanGE %>% inner_join(maxOD, by = 'grp2')

# create df with just the t_maxOD
dOD <- cleanGEwMaxOD %>% filter(tmpt == t_maxOD) %>% 
  dplyr::select(plate_id, well, grp2, exp, grp, medium, strain, t_maxOD, OD600, avgMaxOD) %>%
  rename(ODatMaxAvg = OD600)

# create a second df with just the info at t=0
ODt0 <- cleanGE %>% filter(tmpt == 0) %>% dplyr::select(grp2, 'plate_id', well, OD600) %>% rename(ODt0 = OD600)

# combine the info for maxOD and t=0, calculate the change in OD for each between t=0 and the OD at the time at which average OD was highest
dOD <- dOD %>% 
  inner_join(ODt0, by = c('grp2', 'plate_id', 'well')) %>% 
  mutate(dOD = ODatMaxAvg - ODt0) %>% ungroup() 

rm(ODt0, maxOD)

## determine positive growth ##

# Dunn's test for all pairs of strains on each medium
media <- unique(dOD$medium)
media.dunn <- tibble()

for (m in media){
  dat <- dOD %>% filter(medium == m)
  med.dn <- FSA::dunnTest(dOD ~ strain, dat, method = 'none', two.sided = TRUE)
  x <- med.dn[[2]] %>% separate(Comparison, c("strain1", "strain2"), sep = ' - ') %>% mutate(medium = m)
  media.dunn <- bind_rows(media.dunn, x)
}

# select only the comparisons with the negative control
s1neg <- media.dunn %>% filter(strain1 == 'negative') %>% 
  dplyr::select(medium, strain2, Z, P.unadj) %>%
  rename(strain = strain2)
s2neg <- media.dunn %>% filter(strain2 == 'negative') %>% 
  dplyr::select(medium, strain1, Z, P.unadj) %>%
  rename(strain = strain1)
media.sigVsNeg <- bind_rows(s1neg, s2neg)

rm(m, dat, med.dn, x, s1neg, s2neg, media.dunn)

# correct p.values for multiple comparisons
media.adjPval <- p.adjust(media.sigVsNeg$P.unadj, method = "BH")
media.sigVsNeg$P.adj <- media.adjPval

dOD.pval <- dOD %>% 
  ungroup() %>% 
  dplyr::select(grp, medium, strain) %>% 
  distinct() %>% 
  inner_join(media.sigVsNeg, by = c('medium', 'strain')) %>%
  mutate(growth = ifelse(P.adj < .05, 'yes', 'no')) %>%
  dplyr::select(grp, growth)

rm(media.adjPval)

# create a new df in which strains that don't grow have dOD set to zero
dODfinal <- dOD %>% 
  inner_join(dOD.pval, by = 'grp') %>% 
  dplyr::select(exp, grp2, grp, strain, medium, dOD, growth) %>%
  mutate(dODadj = ifelse(growth == 'no', 0, dOD)) %>%
  mutate(dODadj = ifelse(is.na(dODadj), 0, dODadj))

rm(dOD.pval, media.sigVsNeg, negs, allGE, cleanGEwMaxOD, media, clean_str)

# psych and jeju are the only two strains that do not grow on any medium, remove from clustering
# see format_growth_profiles.Rmd
# remove jeju and psych, get mean of adjusted dODs
dODfinal <- dODfinal %>% 
  filter(!strain %in% c("jeju", "psych")) %>%
  group_by(strain, medium) %>% 
  mutate(deltaOD = mean(dODadj)) %>% 
  ungroup()
# remove jeju and psych, cleanGE used for maxGR
cleanGE <- cleanGE %>%
  filter(!strain %in% c("jeju", "psych"))

# spread data 
deltaOD <- dODfinal %>% 
  dplyr::select(strain, medium, deltaOD) %>% distinct() %>%
  pivot_wider(names_from = medium, values_from = deltaOD) %>%
  dplyr::select(strain, difcoMB, HMBcmpt, HMBpep, HMBaa, HMBlips, HMBoligo, 
                HMBorg, HMBntrl, HMBamisug, HMBacdsug, `HMB--`) %>%
  ungroup()






