require(ggpmisc)
require(tidyverse)
require(broom)
require(mblm)
require(viridis)
require(lme4)

# to get correct labels in fig. legend
labs_bp <- c(
    "0" = "0.1",
    "0.2" = "0.18",
    "0.4" = "0.26",
    "0.6" = "0.34",
    "0.8" = "0.42",
    "1" = "0.5")

# to get nice motif labels
prop.labs <- c("Apparent comp.", "Three-sp. chain", "Direct comp.", "Omnivory")
names(prop.labs) <- c("frac_acomp", "frac_chain", "frac_dcomp", "frac_omni")

# data file
df<-read.table(file = '3sp_roles_participation.tsv', sep = '\t', header = TRUE)

tab_t <- df %>%
    filter(in_Degree > 0) %>%   # remove basal sp.
    unite(NW, Network, Connectance, Size, remove=FALSE) %>%  # unite info about NW, C, S in new column NW
    select(NW, Species, Persistence, m12, m36, m38, m6, 
           Connectance, Basal_p, Size) %>%  # get rid of some columns not necessary
    mutate(omni=m38, chain=m12, dcomp=m36,acomp=m6) %>%  # rename for my sake
    gather("motif","counts",4:7) %>%   
    group_by(NW, Species, Persistence, Basal_p, Connectance, Size,
             omni, chain, dcomp, acomp) %>%
    summarise(total=sum(counts)) %>% # create total of each motif (grouped by above)
    ungroup() %>%
    mutate(frac_omni = omni/total, frac_chain=chain/total,
           frac_dcomp=dcomp/total, frac_acomp=acomp/total) %>%
    gather("norm_motif", "frac", 12:15) # create two columns of motif and frac

# lmer model function
lmer_per_frac <- function(df) with(df,lmer(Persistence~frac+(1|Size:Connectance)+(1|NW)))

# tidyverse extraction 
extract_stats <- function(input) {
    input %>% summary %>% coef %>%
        as_tibble() %>%
        select(1) %>%
        rename(coefs=1) %>%
        mutate(type=c("intercept", "slope")) %>%
        return()
}

# dark nesting of data to increase efficiency, call of lmer function
tab_prop <- tab_t %>% 
    group_by(norm_motif, Basal_p) %>%
    nest() %>%
    ungroup() %>%
    mutate(stats = map(data, lmer_per_frac)) %>%
    #pull(stats) %>% `[[`(1)
    mutate(params = map(stats, extract_stats)) %>%
    unnest(params) %>%
    select(-stats, -data) %>%
    pivot_wider(names_from = type, values_from = "coefs")

# get lines from slopes and intercepts
get_lines <- tab_prop %>%
    mutate(max_frac=c(rep(0.55,6), rep(1,18)), min_frac=0) %>% # create max and min proportion of role, first 6 in max are omnivory which only reaches 0.55
    mutate(min_persistence = (min_frac*slope + intercept)) %>% #min and below max per. for each motif and basal p respectively
    mutate(max_persistence = (max_frac*slope + intercept)) %>%
    gather(type, Persistence, 7:8) %>%
    mutate(frac=ifelse(type=="min_persistence", 0,1)) %>% # to get min and max frac for each motif and basal p (yes could be done neater since info already exists)
    mutate(frac=ifelse(norm_motif=="frac_omni" & type=="max_persistence", 0.55, frac)) %>%
    select(Basal_p, norm_motif, Persistence, frac)

ggplot(get_lines, aes(x=frac, y=Persistence, color=as.factor(Basal_p))) +
    geom_line(size=1) +
    facet_grid(~norm_motif, labeller=labeller(norm_motif=prop.labs)) +
    scale_color_viridis_d(name = "Basal species\nextinction\nprobability", 
                          labels = labs_bp) +
    labs(x="Proportion of role", y="Probability of persistence")   +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size=7)) 
#ggsave("figures/2021/overleaf/prop_lmer_allCS.pdf", width=8, height=6.2)
