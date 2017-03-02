# A demo of dwgLASSO to prioritize metabolite list based on INDEED
#
# Reference:
# [1] Zuo, Yiming, Yi Cui, Cristina Di Poto, Rency S. Varghese, Guoqiang Yu, Ruijiang Li, and 
#     Habtom W. Ressom. "INDEED: Integrated differential expression and differential network 
#     analysis of omic data for biomarker discovery." Methods 111 (2016): 12-20.
#
# Copyright 2016-2017, Yiming Zuo.

## Start clean
rm(list = ls())

## Make sure current working directory is INDEED package
getwd()

## Load data
# metabolomic GU cohort
load("Data/M_GU.RData")
# # metabolomic TU cohort
# load("Data/M_TU.RData")

p <- nrow(Met_HCC_GU)
n_HCC <- ncol(Met_HCC_GU)
n_CIRR <- ncol(Met_CIRR_GU)

## Load R package
# glasso
library("glasso") # need install

## Load functions
source("function.R")

## Z-transform the data for group-specific normalization
data_CIRR <- scale(t(Met_CIRR_GU)) # n1*p CIRR
data_HCC <- scale(t(Met_HCC_GU)) # n2*p HCC
cov_CIRR <- var(data_CIRR)
cov_HCC <- var(data_HCC)

set.seed(100)
## Apply grapihcal LASSO given data_CIRR and data_HCC
# CIRR first
n_fold <- 5 # number of folds
rho = exp(seq(log(0.6), log(0.01), length.out = 20))
# draw error curve 
error <- choose_rho(data_CIRR, n_fold, rho)
# chosse optimal rho 
rho[error$log.cv == min(error$log.cv)] # rho based on minimum rule
abline(v = rho[error$log.cv == min(error$log.cv)], col = "red", lty = 3)
# one standard error rule
abline(h = min(error$log.cv) + error$log.rho[error$log.cv == min(error$log.cv)], col = "blue") 
rho
#################################################################################################
# select rho that are under blue line and on the right of the red line, in this case, rho = 0.253
#################################################################################################
rho_CIRR_opt <- 0.253
# perform gLASSO
pre_CIRR <- glasso(cov_CIRR, rho = rho_CIRR_opt) 
rm(n_fold, rho, error)
# examine the precision matrix
thres <- 1e-3
sum(abs(pre_CIRR$wi) > thres) 
pre_CIRR$wi[1:10, 1:10]
# compute partial correlation
pc_CIRR <- compute_par(pre_CIRR$wi)
# examine the partial correlation matrix
sum(abs(pc_CIRR) > thres) 
pc_CIRR[1:10, 1:10]
rm(pre_CIRR, thres)

## HCC second
n_fold <- 5 # number of folds
rho = exp(seq(log(0.6), log(0.01), length.out = 20))
# draw error curve 
error <- choose_rho(data_HCC, n_fold, rho)
# chosse optimal rho 
rho[error$log.cv == min(error$log.cv)] # rho based on minimum rule
abline(v = rho[error$log.cv == min(error$log.cv)], col = "red", lty = 3)
# one standard error rule
abline(h = min(error$log.cv) + error$log.rho[error$log.cv == min(error$log.cv)], col = "blue") 
rho
#################################################################################################
# select rho that are under blue line and on the right of the red line, in this case, rho = 0.253
#################################################################################################
rho_HCC_opt <- 0.253
# perform gLASSO
pre_HCC <- glasso(cov_HCC, rho = rho_HCC_opt) 
rm(n_fold, rho, error)
# examine the precision matrix
thres <- 1e-3
sum(abs(pre_HCC$wi) > thres) 
pre_HCC$wi[1:10,1:10]
# compute partial correlation
pc_HCC <- compute_par(pre_HCC$wi)
# examine the partial correlation matrix
sum(abs(pc_HCC) > thres) 
pc_HCC[1:10,1:10]
rm(pre_HCC, thres)

## Build differential networks
diff <- pc_CIRR - pc_HCC
thres = 1e-3
sum(abs(diff) > thres)
diff[1:10, 1:10]

## Permutation test
# This step takes some time. If you have done this before and no parameters are needed to change, 
# you can simply load the saved result
m <- 1000 # permutation times
diff_p <- permutation(m, p, n_CIRR, n_HCC, data_CIRR, data_HCC, rho_CIRR_opt, rho_HCC_opt)
save(diff_p, file = "Data/diff_p_M_GU.RData")
rm(m, thres, rho_CIRR_opt, rho_HCC_opt)
# load the saved result
load("Data/diff_p_M_GU.RData")
# calculate differential network connections
thres_left <- 0.025
thres_right <- 0.975
pc_thres <- permutation_thres(thres_left, thres_right, p, diff_p)
rm(thres_left, thres_right)
# get binary matrix
pc_thres_p <- pc_thres$positive
pc_thres_n <- pc_thres$negative
pc_binary <- matrix(0, p, p)
pc_binary[diff < pc_thres_n] <- -1
pc_binary[diff > pc_thres_p] <- 1
sum(diff < pc_thres_n)
sum(diff > pc_thres_p)
pc_binary[1:10, 1:10]
rowSums(abs(pc_binary)) # node degree for differential networks
rm(diff_p)

## Convert adjacent matrix into edge list
edge <- matrix(0, (sum(diff < pc_thres_n) + sum(diff > pc_thres_p)) / 2, 3)
k <- 1
for(i in 1:(nrow(pc_binary) - 1)){
    for (j in (i + 1) : nrow(pc_binary)){
        if(pc_binary[i, j] != 0){
            edge[k, 1] <- i
            edge[k, 2] <- j
            edge[k, 3] <- pc_binary[i, j]
            k <- k + 1
        }
    }
}
edge_dn <- data.frame("Met1" = edge[, 1], "Met2" = edge[, 2], "Weight" = edge[, 3])
write.csv(edge_dn, file = "Met_dn.csv", quote = F, row.names = F)
rm(i, j, k, edge, edge_dn, pc_thres, rho_CIRR_opt, rho_HCC_opt, thres)

## Load in the significant list
sig_l <- read.table("Data/pvalue_M_GU.csv", header = TRUE, sep = ",", stringsAsFactors = F)
# trasfer p-value to z-score
z_score <- abs(qnorm(1 - sig_l$p.value / 2)) 
# calculate differntial network score
dn_score <- compute_dns(pc_binary, z_score)
dns_idx <- sort(dn_score, index.return = T, decreasing = T)$ix # prioritize the list
rowSums(abs(pc_binary))[dns_idx] # the node degree for each metabolite

## save top 10 metabolite after prioritization based on INDEED
write.table(sig_l[dns_idx[1:10],], file = "top10table_INDEED.csv", sep=",", quote = F, 
            row.names = F, col.names = F)