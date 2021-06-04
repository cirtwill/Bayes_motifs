require(XRJulia) ## package allowing for Julia interfacing
require(igraph)
require(tidyverse) ## for efficient data manipulation and plotting
source("BN-DFS.R") ## to make the adjacency matrix acyclic
source("getTL.R") ## to get prey average TL, filter good webs

## sort topological and remove cycles
topo_sort_acyclic <- function(A) {
    A<-as.matrix(A)
    ## Make sure producers are the first rows; this is needed for the DFS
    ## algorithm to work (which removes all cycles from the digraph)
    A <- A[order(rowSums(A)),order(rowSums(A))]
    ## Now make the digraph acyclic
    A <- DFS(A)
    
    ## Find a sorting of A's rows and columns to make A lower triangular
    o <- graph_from_adjacency_matrix(A) %>% topo_sort(mode="in")
    ## Return sorted, lower triangular adjacency matrix A
    return(A[o,o])
}

# 50sp webs done, 60sp webs done,
sp_nos <- seq(50,100, by=10)
c <- seq(0.02,0.18,by=0.04)
web_nos <- c(0:649)
for (i in 1:length(sp_nos)){
    sp_no <- sp_nos[i]
    print(c(sp_no,'Species'))
    for (j in 1:length(c)){
        connectance <- c[j]
        print(c(connectance,'Connectance'))
        goodfiles=as.character(list.files(path=paste0("../data/networks/acyclic_webs/",as.character(sp_no),"/",as.character(connectance)), full.names=FALSE, pattern="*.csv"))
        if(length(goodfiles)<100){ # Keep addding files up to 100
            print('incomplete')
            n=length(goodfiles)
            for (k in n:max(web_nos)){
                if(!paste("initial_matrix_",k,".csv",sep='')%in%goodfiles && n<100){
                    A <- read.table(paste("../data/networks/pre_disturbance/",sp_no, "/", connectance, "/", "initial_net_",
                                          k,".csv",sep=''), sep=",", header=TRUE)
                    rownames(A) <- colnames(A) <- 1:nrow(A) 
                    image(t(apply(A,2,rev)))
                    
                    A <- topo_sort_acyclic(A)
                    image(t(apply(A,2,rev)))

                    # Calculate the PATL to sort good webs
                    FW <- list()
                    FW$M <- t(A)
                    FW <- getTL(FW)
                    preyav_tl=round(FW$TL,2)
                    # Let's only keep the ones with max PATL<6
                    # Also indexing by "good webs" so I know how many to re-run
                    if(max(preyav_tl)<6){
                        if(n<100){
                            write.table(A, file=paste("../data/networks/acyclic_webs/",sp_no, "/", connectance, "/", "initial_matrix_",
                                                k,".csv",sep=''), sep=",")
                            print(c(n,k))
                            n=n+1
                        }
                    }
                }
            }
        }
    }
}




