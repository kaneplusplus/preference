# R code to generate the data sets in the preference package.

# ORIGINAL CODE TO GENERATE DATA. NO LONGER USED...
# # Unstratified summary data
# imap <- data.frame(matrix(NA, nrow=4, ncol=5))
# colnames(imap) <- c("mean", "sd", "n", "trt", "arm")
# imap$mean <- c(47.57283, 50.58991, 46.16386, 45.51061)
# imap$sd <- c(10.6162889, 4.8686232, 9.7385777, 9.9744353)
# imap$n <- c(50, 22, 76, 64)
# imap$trt <- as.factor(c("HPV", "Pap", "HPV", "Pap"))
# imap$arm <- as.factor(c("C", "C", "R", "R"))
# save(imap, file="../data/imap.rda")
# 
# # Stratified summary data (stratified by STAI score)
# imap_strat <- data.frame(matrix(NA, nrow=8, ncol=6))
# colnames(imap_strat) <- c("mean", "sd", "n", "trt", "arm", "stratum")
# imap_strat$mean <- c(54.3337592, 41.0085061, 51.5041323, 49.8965119, 
#   52.0043692, 42.179696, 53.8427658, 41.1461439)
# imap_strat$sd <- c(5.5000253, 10.5061654, 6.1202786, 3.8040966, 6.5603227, 
#   9.6234748, 5.2579869, 9.0525148)
# imap_strat$n <- c(24, 25, 10, 11, 30, 44, 22, 42)
# imap_strat$trt <- as.factor(c("HPV", "HPV", "Pap", "Pap", "HPV", "HPV", 
#   "Pap", "Pap"))
# imap_strat$arm <- as.factor(c("C", "C", "C", "C", "R", "R", "R", "R"))
# imap_strat$stratum <- as.factor(c(1, 2, 1, 2, 1, 2, 1, 2))
# save(imap_strat, file="../data/imap_strat.rda")

imap <- read.csv("imap.csv")[,-1]
imap$arm[imap$arm == FALSE] <- "choice"
imap$arm[imap$arm == "TRUE"] <- "random"
save(imap, file="../data/imap.rda")

imap_summary <- read.csv("imap_summary.csv")[,-1]
imap_summary$arm[imap_summary$arm == FALSE] <- "choice"
imap_summary$arm[imap_summary$arm == "TRUE"] <- "random"
save(imap_summary, file="../data/imap_summary.rda")

imap_stratified_summary <- read.csv("imap_stratified_summary.csv")[,-1]
imap_stratified_summary$arm[imap_stratified_summary$arm == FALSE] <- "choice"
imap_stratified_summary$arm[imap_stratified_summary$arm == "TRUE"] <- "random"
save(imap_stratified_summary, file="../data/imap_stratified_summary.rda")
