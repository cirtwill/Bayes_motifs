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

webs=as.character(list.files(path="../data/empirical/global_verts/matrix/", full.names=FALSE, pattern="*.tsv"))
disturbance <- seq(0,1,by=0.2)     # Simulated level of proportion of basal species disturbed, transforms in a later
                                   # stage to a disturbance level of 0.1 to 0.5 in steps of 0.1
for(style in c('Nonlinear','Linear')){
    for (k in 1:length(webs)){
        web_no <- strsplit(webs[k],'.tsv')[[1]][1]
        print(web_no)
        A <- read.table(paste("../data/empirical/global_verts/acyclic_matrix/",
                              web_no,'.tsv',sep=''), sep=",", header=TRUE) # topo sorted, acyclic graph
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
            if(style=='Nonlinear'){
                persistence <- tibble(
                    sp_no=rownames(A), ## species numbers
                    per = bayes_ntw_MC(A, Pb, 3, 3, 100000) %>% # input: number of species, 
                                                                # probability of extinction vector, 
                                                                # alpha and beta for the functional response of consumers, 
                                                                # 100k evals of BN
                        juliaGet                                ## convert vector from Julia to R format
                )  
            } else {
                persistence <- tibble(
                    sp_no=rownames(A), ## species numbers
                    per = bayes_ntw_MC(A, Pb, 1, 1, 100000) %>% # input: number of species, 
                                                                # probability of extinction vector, 
                                                                # alpha and beta for the functional response of consumers, 
                                                                # 100k evals of BN
                        juliaGet                                ## convert vector from Julia to R format
                )                  
            }
            # connect the simulated persistences to df2 
            df2 <- bind_cols(select(persistence, per),df2)  
            
            # collect the results from various disturbance levels into a single df,
            # which are being saved (one df per network id, S, C)
            if (l==1) {
                collect<-df2
            } else{
                collect <- bind_rows(collect, df2)
            }
            print(dim(collect))
            if(style=='Nonlinear'){
                save(collect, file=paste0("Results_restricted/Nonlinear/global_verts/", web_no, ".RData"))            
            } else {
                save(collect, file=paste0("Results_restricted/Linear/global_verts/", web_no, ".RData"))            
            }
        }
    }
}
