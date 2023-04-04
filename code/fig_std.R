require(ggpmisc)
require(tidyverse)
require(broom)
require(mblm)
require(viridis)
require(lme4)

#df <- read.table(file = 'empirical_3sp_roles_participation_nonlinear.tsv', sep = '\t', header = TRUE)
df<-read.table(file = '../data/3sp_roles_participation.tsv', sep = '\t', header = TRUE)


labs_bp <- c(
    "0" = "0.1",
    "0.2" = "0.18",
    "0.4" = "0.26",
    "0.6" = "0.34",
    "0.8" = "0.42",
    "1" = "0.5")

prop.labs <- c("Apparent comp.", "Three-sp. chain", "Direct comp.", "Omnivory")
names(prop.labs) <- c("frac_acomp", "frac_chain", "frac_dcomp", "frac_omni")

## modification of df
tab_t <- df %>% 
    filter(in_Degree > 0, Size %in% c(50,100),
           Connectance==0.02) %>%
    filter(in_Degree > 0) %>% 
    unite(NW, Network, Connectance, Size, remove=FALSE) %>%
    select(NW, Species, Persistence, m12, m36, m38, m6, 
           Connectance, Basal_p, Size) %>%
    mutate(omni=m38, chain=m12, dcomp=m36, acomp=m6) %>%
    gather("motif","counts", 4:7) %>%
    group_by(NW, Species, Persistence, Basal_p, Connectance, Size,
             omni, chain, dcomp, acomp) %>%
    summarise(total=sum(counts)) %>% 
    ungroup() %>%
    mutate(frac_omni = omni/total, frac_chain=chain/total,
           frac_dcomp=dcomp/total, frac_acomp=acomp/total) %>%
    gather("norm_motif", "frac", 12:15) 


# tab_temp%>%filter(NW=="kongsfjorden_0.024_260", Species=="sp10", Basal_p==1.0)
#lmer_per_frac <- function(df) with(df,lmer(Persistence~frac+(1|Size:Connectance)+(1|NW)))
lmer_per_frac <- function(df) with(df,lmer(Persistence~frac+(1|NW)))

## extract stats (hopefully std...)
extract_stats_std <- function(input) {
    input %>% summary %>% coef %>%
        as_tibble() %>%
        select(2) %>%
        rename(coefs=1) %>%
        mutate(type=c("std_error")) %>%
        return()
}

## extract stats (this I know is correct, extracts intercept and slope)
extract_stats <- function(input) {
    input %>% summary %>% coef %>%
        as_tibble() %>%
        select(1) %>%
        rename(coefs=1) %>%
        mutate(type=c("intercept", "slope")) %>%
        return()
}

## modification of df with std-info
tab_prop_std <- tab_t %>% 
    group_by(norm_motif, Basal_p) %>%
    nest() %>%
    ungroup() %>%
    mutate(stats = map(data, lmer_per_frac)) %>%
    #pull(stats) %>% `[[`(1)
    mutate(params = map(stats, extract_stats_std)) %>%
    unnest(params) %>%
    select(-stats, -data) %>%
    #slice(-seq(2, 48, 2)) %>%
    rename_with(.cols = 3, ~"std_error") %>%
    select(-"type")

## modification of df with other stat info (should be correct)
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

tab_all <- tab_prop %>%
    left_join(tab_prop_std, by=c("Basal_p", "norm_motif"))

# manipulations to get lines for all different combinations...
get_lines <- tab_all %>%
    mutate(intercept_max = intercept + std_error, intercept_min = intercept - std_error) %>%
    select(-"std_error") %>%
    mutate(max_frac=c(rep(0.55,12), rep(1,36)), min_frac=0) %>%
    mutate(val_min_frac = (min_frac*slope + intercept)) %>%
    mutate(val_max_frac = (max_frac*slope + intercept)) %>%
    mutate(val_min_frac_std_max = (min_frac*slope + intercept_max)) %>%
    mutate(val_max_frac_std_max = (max_frac*slope + intercept_max)) %>%
    mutate(val_min_frac_std_min = (min_frac*slope + intercept_min)) %>%
    mutate(val_max_frac_std_min = (max_frac*slope + intercept_min)) %>%
    gather(type, Persistence, 9:14) %>%
    mutate(frac=ifelse(type=="val_min_frac" | type=="val_min_frac_std_max" | 
                           type=="val_min_frac_std_min", 0,1)) %>%
    mutate(frac=ifelse(norm_motif=="frac_omni" & frac==1, 0.55, frac)) %>%
    select(Basal_p, norm_motif, Persistence, frac, type) %>%
    mutate(type=ifelse(type=="val_min_frac" | type=="val_max_frac", "original", type)) %>%
    mutate(type=ifelse(type=="val_min_frac_std_min" | type=="val_max_frac_std_min", "std_min", type)) %>%
    mutate(type=ifelse(type=="val_min_frac_std_max" | type=="val_max_frac_std_max", "std_max", type)) 

ggplot(get_lines) +
    geom_line(aes(x=frac, y=Persistence, color=as.factor(Basal_p), linetype=type), size=1) +
    facet_grid(~norm_motif, labeller=labeller(norm_motif=prop.labs)) +
    scale_color_viridis_d(name = "Basal species\nextinction\nprobability", 
                          labels = labs_bp) +
    labs(x="Proportion of role", y="Probability of persistence")   +
    scale_linetype_manual(values=c("solid", "dashed", "dashed")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size=7)) 

# ggplot(get_lines, aes(x=frac, y=Persistence, color=as.factor(Basal_p))) +
#     stat_smooth(method="loess", se=TRUE, aes(fill=as.factor(Basal_p)), alpha=0.3) +
#     facet_grid(~norm_motif, labeller=labeller(norm_motif=prop.labs)) +
#     scale_color_viridis_d(name = "Basal species\nextinction\nprobability", 
#                           labels = labs_bp) +
#     labs(x="Proportion of role", y="Probability of persistence")   +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 90, size=7)) 

# ggplot(get_lines, aes(x=frac, y=Persistence, color=as.factor(Basal_p))) +
#     stat_summary(geom="ribbon", fun.min="min", fun.max="max", aes(fill=as.factor(Basal_p))) +
#     facet_grid(~norm_motif, labeller=labeller(norm_motif=prop.labs)) +
#     scale_color_viridis_d(name = "Basal species\nextinction\nprobability", 
#                           labels = labs_bp) +
#     labs(x="Proportion of role", y="Probability of persistence")   +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 90, size=7)) 
# 
# ggplot(get_lines) +
#     geom_ribbon(aes(x=frac, y=Persistence), fill=as.factor(Basal_p))
#     stat_summary(geom="ribbon", fun.min="min", fun.max="max", aes(fill=as.factor(Basal_p), shape=type)) +
#     facet_grid(~norm_motif, labeller=labeller(norm_motif=prop.labs)) +
#     scale_color_viridis_d(name = "Basal species\nextinction\nprobability", 
#                           labels = labs_bp) +
#     labs(x="Proportion of role", y="Probability of persistence")   +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 90, size=7)) 
# 
# 
# 
# ggplot(get_lines, aes(x=frac, y=Persistence, color=as.factor(Basal_p))) +
#     geom_line(size=1) +
#     facet_grid(~norm_motif, labeller=labeller(norm_motif=prop.labs)) +
#     #scale_color_manual(name = "Basal species\nextinction\nprobability", labels = labs_bp,
#     #                   values=color_sc) +
#     scale_color_viridis_d(name = "Basal species\nextinction\nprobability", 
#                           labels = labs_bp) +
#     labs(x="Proportion of role", y="Probability of persistence")   +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 90, size=7)) 
# ggsave("figures/2021/overleaf/nl_empirical_prop_lmer_allCS.pdf", width=8, height=6.2)
