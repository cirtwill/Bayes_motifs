#### TEST MOTIFS AND PERSISTENCE 
require(ggpmisc)
require(tidyverse)

color_sc <- c("#CC79A7", "#E69F00", "#56B4E9")
labs_bp <- c(
    "0" = "0 %",
    "0.4" = "40 %",
    "0.8" = "80 %")

data<-read.table(file = '3sp_roles_participation.tsv', sep = '\t', header = TRUE)

cc <- c(0.02, 0.1, 0.18)
bp <- c(0,0.4,0.8)

tab2 <- data %>%
    filter(in_Degree > 0, Size %in% c(50,100), Basal_p %in% bp, 
           Connectance %in% cc) %>%
    mutate(preyav_length=cut(PATL, breaks=c(1.9, 3, 4, 5, 6), 
                             labels=c("2","3","4","5"))) %>% 
    select(Network, Species, Persistence, m12, m36, m38, m6, 
           preyav_length, Connectance, Basal_p, Size) %>%
    mutate(omni=m38, chain=m12, dcomp=m36,acomp=m6) %>%
    gather("motif","counts",4:7) %>%
    group_by(Network, Species, Persistence, preyav_length, Basal_p,
             omni, chain, dcomp, acomp) %>%
    summarise(total=sum(counts)) %>%
    ungroup() %>%
    mutate(frac_omni = omni/total, frac_chain=chain/total,
           frac_dcomp=dcomp/total, frac_acomp=acomp/total) %>%
    gather("norm_motif", "frac", 11:14) 

ggplot(tab2, aes(x=frac, y=Persistence, color=as.factor(Basal_p))) +
    #geom_point(size=1, alpha=0.5) +
    facet_grid(~norm_motif) +
    theme(axis.text.x = element_text(angle = 90, size=7)) +
    scale_color_manual(name = "Basal prop.", labels = labs_bp, 
                       values=color_sc) +
    stat_smooth(method="lm", se=FALSE) +
    stat_fit_glance(method="lm", geom="text", label.x=c(0.7,0.7,0.7),
                    label.y = c(0.58, 0.56, 0.54),
                    aes(label=sprintf('italic(p)~"="~%.3f',
                                      stat(..p.value..))),
                    parse=TRUE,size = 3)

### RAW COUNTS
tab2 <- data %>%
    filter(in_Degree > 0, Size %in% c(50,100), Basal_p %in% bp, 
           Connectance %in% cc) %>%
    mutate(preyav_length=cut(PATL, breaks=c(1.9, 3, 4, 5, 6), 
                             labels=c("2","3","4","5"))) %>% 
    select(Size, Connectance, Basal_p, Network, Species, Persistence, 
           m12, m36, m38, m6) %>% 
    rename(omni=m38, chain=m12, dcomp=m36,acomp=m6) %>%
    gather("motif","counts",7:10) 

ggplot(tab2, aes(x=counts, y=Persistence, color=as.factor(Basal_p))) +
    #geom_point(size=1, alpha=0.5) +
    facet_grid(~motif) +
    theme(axis.text.x = element_text(angle = 90, size=7)) +
    scale_color_manual(name = "Basal prop.", labels = labs_bp, 
                       values=color_sc) +
    stat_smooth(method="lm", se=FALSE) +
    stat_fit_glance(method="lm", geom="text", label.x=c(2700,2700,2700),
                    label.y = c(0.76, 0.74, 0.72),
                    aes(label=sprintf('italic(p)~"="~%.3f', 
                                      stat(..p.value..))),
                    parse=TRUE,size = 2.5) 
ggsave("figures/motif_stats/websize_100/raw_motif.pdf", width=8, height=6.2)

