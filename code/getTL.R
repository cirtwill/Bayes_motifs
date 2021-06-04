# ADD THIS IN MAIN CODE!
# FW <- list()
# FW$M <- t(A)
# FW <- getTL(FW)

normalizeM <- function(M){
    colsum_M <- colSums(M);
    colsum_M[colsum_M==0] <- 1;
    return(t(t(M)/colsum_M));
}

getTL <- function(FW){
    S <- dim(FW$M)[1] # number of species
    ## normalize M such that a predator has no preference over
    ## multiple prey.
    ##Ref: A process-oriented approach to the multispecies functional response
    ##  Koen-Alonso
    ## From energetics to ecosystems: The dynamics and structure of ecological systems Springer (2007)
    M <- normalizeM(FW$M)
    
    ## Since: 
    ##
    ##    	TL[j] = 1 + sum_over_i ( M[i, j] TL[i] )
    ## 
    ## which written in matrix notation is:
    ##
    ##		TL = 1 + t(M) %*% TL
    ## =>  (I - t(M) ) %*% TL = 1    <- this is a vector
    ## =>  TL = (I - t(M))^(-1) %*% rep(1, S);
    ## if (I - t(M)) is not singular, 
    ## then it is fine to do the matrix inversion
    if(det( diag(S)-t(M) )!=0){
        FW$TL <- solve(diag(S)-t(M), rep(1, S))
    } else{
        ## Otherwise, need to use
        ## I + t(M) + t(M)%*% t(M) + t(M)%*%t(M) %*% t(M) + ... to approximate
        ## (I - t(M))^(-1)
        tmp <- diag(S)
        for(i in 1:9){
            tmp <- tmp  %*% t(M) + diag(S)
        }
        FW$TL <- tmp %*% rep(1, S)
    }
    FW$W <- M ## Save the normalized adjacency matrix as weight matrix
    return(FW)
}




