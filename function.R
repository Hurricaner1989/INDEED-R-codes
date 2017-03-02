# Functions defined to help demo.R
#
# Reference:
# [1] Zuo, Yiming, Yi Cui, Cristina Di Poto, Rency S. Varghese, Guoqiang Yu, Ruijiang Li, and 
#     Habtom W. Ressom. "INDEED: Integrated differential expression and differential network 
#     analysis of omic data for biomarker discovery." Methods 111 (2016): 12-20.
#
# Copyright 2016-2017, Yiming Zuo.

## Create log likelihood error function 
loglik_ave <- function(data, theta){
    loglik <- c()
    loglik <- log(det(theta)) - sum(diag(var(data) %*% theta))
    return(-loglik)
}

## Draw error curve
choose_rho <- function(data, n_fold, rho) {
    # randomly shuffle the data
    Data <- data[sample(nrow(data)), ]
    # create n_fold equally size folds
    folds <- cut(seq(1, nrow(Data)), breaks = n_fold, labels = FALSE)
    # tune parameters
    d <- ncol(Data)
    loglik_cv <- c()
    loglik_rho <- c()
    pb <- txtProgressBar(min = 0, max = length(rho), style = 3) # create progress bar
    for(i in 1 : length(rho)){
        Sys.sleep(0.1)
        # perform n_fold cross validation
        loglik <- c()
        for(j in 1 : n_fold){
            # segement your data by fold using the which() function 
            testIndexes <- which(folds == j, arr.ind = TRUE)
            testData <- Data[testIndexes, ]
            trainData <- Data[-testIndexes, ]
            # use test and train data partitions however you desire...
            cov <- var(trainData) # compute the covariance matrix
            pre<- glasso(cov, rho = rho[i])
            loglik <- c(loglik, loglik_ave(testData, pre$wi))  
        }
        loglik_cv <- c(loglik_cv, sum(loglik) / n_fold)
        loglik_rho <- c(loglik_rho, sd(loglik) / sqrt(n_fold))
        setTxtProgressBar(pb, i) # update progress bar
    }
    close(pb)
    plot(rho, loglik_cv, main = "Error curve using corss validation",
         xlab = expression(lambda), ylab = "Error")
    lines(rho, loglik_cv)
    error <- list("log.cv" = loglik_cv, "log.rho" = loglik_rho)
    return(error)
}

## Compute partial correlation
compute_par <- function(pre_inv) {
    p <- nrow(pre_inv)
    pc <- matrix(0, p, p) # partial correlation
    for(i in 1:(p - 1)){
        for(j in (i + 1) : p){ 
            pc[i, j] = -pre_inv[i, j] / sqrt(pre_inv[i, i] * pre_inv[j, j])
            pc[j, i] = pc[i, j]
        }
    }
    return(pc)
}

## Permutation to build differential network 
permutation <- function(m, p, n_CIRR, n_HCC, data_CIRR, data_HCC, rho_CIRR_opt, rho_HCC_opt) {
    diff_p <- array(0, dim = c(m, p, p))
    pb <- txtProgressBar(min = 0, max = m, style = 3)
    for(t in 1 : m) {
        Sys.sleep(0.1)
        data_CIRR_p <- matrix(0, n_CIRR, p)
        for(i in 1 : p) {
            data_CIRR_p[, i] <- data_CIRR[sample(n_CIRR), i]
        }
        data_HCC_p <- matrix(0, n_HCC, p)
        for(i in 1 : p) {
            data_HCC_p[, i] <- data_HCC[sample(n_HCC), i]
        }
        per_CIRR <- glasso(var(data_CIRR_p), rho = rho_CIRR_opt)
        per_HCC <- glasso(var(data_HCC_p), rho = rho_HCC_opt)
        pc_CIRR_p <- compute_par(per_CIRR$wi)
        pc_HCC_p <- compute_par(per_HCC$wi)
        diff_p[t, , ] <- pc_CIRR_p - pc_HCC_p
        # update progress bar
        setTxtProgressBar(pb, t)
    } 
    close(pb)
    return(diff_p)
}

## Calculate the positive and negative threshold based on the permutation result
permutation_thres <- function(thres_left, thres_right, p, diff_p) {
    pc_thres_p <- matrix(0, p, p)
    pc_thres_n <- matrix(0, p, p)
    for(i in 1 : (p-1)){
        for(j in (i + 1) : p){ 
            pc_thres_n[i, j] <- quantile(diff_p[, i, j], probs = thres_left)
            pc_thres_n[j, i] <- pc_thres_n[i, j]
            pc_thres_p[i, j] <- quantile(diff_p[, i, j], probs = thres_right)
            pc_thres_p[j, i] <- pc_thres_p[i, j]
        }
    }
    pc_thres <- list("positive" = pc_thres_p, "negative" = pc_thres_n)
    return(pc_thres)
}

## Calculate differential network score
compute_dns <- function(pc_binary, z_score){
    # get adjacent matrix
    diff_d <- abs(pc_binary)
    # set diagonal elements to 1
    diag(diff_d) <- 1
    # compute differential network score for each row
    dns <- apply(diff_d, 1, function(x, y = z_score) sum(y[which(x == 1)]))
    return(dns)
}
