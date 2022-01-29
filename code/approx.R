require(XRJulia) ## package allowing for Julia interfacing
require(tidyverse) ## for efficient data manipulation and plotting
require(igraph)
library(dplyr)
source("getTL.R") ## to get prey average TL

juliaUsing("Random, Distributions")

## Monte Carlo evaluation of Bayesian network, written in Julia
bayes_ntw_MC_jul <- juliaEval("
                              function bayesianNetworkMonteCarlo(A, pBaseline, a, b, reps)
                              sA = sum(A, dims = 2)
                              nRows = length(sA)
                              lastBasal = maximum(findall(sA .== minimum(sA)))[1][1]
                              marginals = zeros(nRows)
                              for rep in 1:1:reps
                              extant = ones(nRows)
                              for i in 1:1:lastBasal
                              extant[i] = rand() < pBaseline[i] ? 0 : 1
                              end
                              for i in (lastBasal + 1):1:size(A, 1)
                              frac = 1 - (extant' * A[i,:])[1] / sum(A[i,:])
                              pExt = pBaseline[i] + (1 - pBaseline[i]) * cdf(Beta(a, b), frac)
                              extant[i] = rand() < pExt ? 0 : 1
                              end
                              marginals = marginals + extant
                              end
                              return(vec(marginals) / reps)
                              end
                              ")

## turn the above Julia code to an R function
bayes_ntw_MC <- JuliaFunction(bayes_ntw_MC_jul)

species <- seq(50, 100, by=10)     # Size of simulated networks, from 50 to 100 in steps of 10
c <- seq(0.02, 0.18, by=0.04)      # Connectance level of simulated networks, from 0.02 to 0.18 in steps of 0.04   
webs <- c(1:100)                   # Web number, total of 100 simulated networks
disturbance <- seq(0,1,by=0.2)     # Simulated level of proportion of basal species disturbed, transforms in a later
                                   # stage to a disturbance level of 0.1 to 0.5 in steps of 0.1



for (i in 1:length(species)){
    sp_no <- species[i]
    for (j in 1:length(c)){
        connectance <- c[j]
        for (k in 1:length(webs)){
            web_no <- webs[k]
            A <- read.table(paste("websR/",sp_no, "/", connectance, "/", "initial_matrix_",
                                  web_no,".csv",sep=''), sep=",", header=TRUE) # topo sorted, acyclic graph
            S <- dim(A)[1]    # no of species in food web A
            A <- as.matrix(A) # rows: predators, columns: prey, lower triangular
            colnames(A)<-rownames(A)
            dd <- distances(graph_from_adjacency_matrix(A), mode="out")
            basals <- rownames(A[rowSums(A)==0,])
            short_tl <- 1 + as.numeric(apply(dd[,basals], 1, min))
            FW <- list()
            FW$M <- t(A) # upper triangular
            FW <- getTL(FW) 
            
            # collect df with each species in- and out-degree, STL
            df <- tibble(
                sp_no = rownames(A), 
                in_degree = rowSums(A),
                out_degree = colSums(A),
                shortest_tl = short_tl
            )
            
            for (l in 1:length(disturbance)){
                dist_level <- disturbance[l]
                df2 <- mutate(df, basal_p = dist_level)
                # get the identity of basal sp, which will be disturbed
                pertube <- df %>%
                    filter(shortest_tl==1) %>%
                    pull(sp_no)
                
                
                per_vec <- match(pertube,rownames(A))   
                Pb <- rep(0.1,S) ## first step baseline extinction probabilities: same for all species 
                Pb[per_vec] <- 0.4*dist_level + 0.1 # the previously selected basal sp get higher baseline extinction prob
                
                # persistence calculations
                persistence <- tibble(
                    sp_no=rownames(A), ## species numbers
                    per = bayes_ntw_MC(A, Pb, 1, 1, 100000) %>% # input: number of species, 
                                                                # probability of extinction vector, 
                                                                # alpha and beta for the functional response of consumers, 
                                                                # 100k evals of BN
                        juliaGet                                ## convert vector from Julia to R format
                )  
                # connect the simulated persistences to df2 
                df2 <- bind_cols(select(persistence, per),df2)  
                
                # collect the results from various disturbance levels into a single df,
                # which are being saved (one df per network id, S, C)
                collect <- df2 %>%
                    select(sp_no, in_degree, out_degree, shortest_tl, per, basal_p) 
                if (l==1) {
                    collect<-df2
                } else{
                    collect <- bind_rows(collect, df2)
                }
                save(collect, file=paste0("Results_restricted/Linear/", sp_no, "/", connectance, "/TESTweb_",web_no, ".RData"))
            }
        }
    }
}

