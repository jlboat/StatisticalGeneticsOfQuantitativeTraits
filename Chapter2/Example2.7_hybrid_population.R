
continuity_corrected_chisq <- function(n_obs, n_exp){
    return((((n_obs-n_exp)^2))/n_exp)
}

chisq_stat <- 0

for (i in c(44,59)){
    chisq_stat = chisq_stat + continuity_corrected_chisq(i,50)
}
print(chisq_stat)
