library(ggplot2)
library(dplyr)
library(tidyverse)


n <- 1000

sim_out <- lapply(1:100, function(x){
  
  sim_out1 <- lapply(1:n, function(j){
    
    dummy_genes <- seq(1,3989)
    
    sample_out <- sort(sample(dummy_genes, 564))
    
    diff_sample_out <- diff(sample_out)
    
    diff_sample_out <- replace(diff_sample_out, diff_sample_out != 1, 0)
    
    diff_sample_out <- c(0, diff_sample_out)
    
    diff.sum <- c(rep(0, length(diff_sample_out)))
    diff.sum
    x <- 0
    
    for(i in 1:length(diff_sample_out)) {
      
      
      if(i < x) { next }
      
      else {
        if(diff_sample_out[i] == 1) {
          d <- -1
          x <- i
          
          
          while(diff_sample_out[x] == 1) { 
            d <- d + 1
            x <- x + 1
            
            if(x > length(diff_sample_out)) { break }
          } 
          diff.sum[(i-1):(i+d)] <- d+2
          
        } else {
          diff.sum[i] <- 0
        }
      }
    }
    
    
    diff.sum <- replace(diff.sum, diff.sum == 0, 1)
    
    out <- 
      (table(diff.sum)/as.numeric(names(table(diff.sum))))/sum(table(diff.sum)/as.numeric(names(table(diff.sum))))
    
    #out <- bind_rows(table(1:50), out)
    #out[-1,]
    
    
  })
  
  
  sim_out1 <- bind_rows(sim_out1)
  
})

sim_out <- bind_rows(sim_out)



zeros_list <- lapply(1:ncol(sim_out), function(x){0})
names(zeros_list) <- names(sim_out)
sim_out <- sim_out %>% replace_na(zeros_list)

# Results of Simulation
sapply(sim_out, mean)

#      1            2            3            4            5            6            7            8            9           10
# 8.586554e-01 1.213703e-01 1.716458e-02 2.418907e-03 3.363765e-04 4.638356e-05 6.945458e-06 1.022033e-06 1.036373e-07 2.109705e-08 


### Bonferroni-Holm Adjustment for alpha-Testniveau
# total of 294 cluster of consecutive essential genes:
#  1: 195
#  2: 52
#  3: 17
#  4: 14
#  5: 4
#  6: 1
#  7: 3
#  8: 2
#  9: 2
# 10: 1
# 11: 1
# 14: 1
# 42: 1
# ---------
#     294

# Global Testniuveau:
Alpha <- 0.05

# B.-H. Adjustment:
Alphas_BH <- Alpha/(294:1)


#     42(1)     |      14(1)    |      11(1)    |      10(1)   |      9(2)    |      8(2)    |      7(3)    |       6(1)   |       5(4)   |      4(14)   |     3(17) 
# <2.109705e-08 | <2.109705e-08 | <2.109705e-08 | 2.109705e-08 | 1.036373e-07 | 1.022033e-06 | 6.945458e-06 | 4.638356e-05 | 3.363765e-04 | 2.418907e-03 | 1.716458e-02 
#  <=0.05/294   |  <=0.05/293   |  <=0.05/292   |  <=0.05/291  |  <=0.05/290  |  <=0.05/288  |  <=0.05/286  |  <=0.05/283  |  >0.05/282   |  >0.05/278   |  >0.05/264





