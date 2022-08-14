library(ggplot2)
library(dplyr)
library(tidyverse)

# Determine all IS of KEGG genes
# Determine Insertion Orientation (plus, minus) 
# Determine IS frame (+1, +2, +3)

# Set working directory
setwd("set/working/directory")


# load azo genes of KEGG gene database
kegg.genes <- readRDS("./0_Data/KEGG_azo_all_genes_data_frame.rds")

kegg.genes <- as_tibble(kegg.genes)


#### SM+N microaerob with index 1, 2, 3


# Laod Input Library (SM+N aerob)
singl.index <- c("Ferm_I_01", "Ferm_I_02", "Ferm_I_03")

for(i in 1:length(singl.index)) {
  
  TN5 <- singl.index[i]
  
  # Load library specific IS 
  Ins_Data <- read_csv(paste("./0_Data/Tn5_Ferm_I_01-I_03/Tn5_", TN5, 
                             "_list_ins_combined.csv", sep=""))
  
  
  
  translation_data <- lapply(1:nrow(kegg.genes), function(x){
    
    tmp_1 <- filter(Ins_Data, 
                    Insertion >=  kegg.genes$gene_start[x] & 
                      Insertion <=  kegg.genes$gene_end[x]
    ) 
    
    if(nrow(tmp_1)==0){
      
      tmp_1 <- NULL
      
    }else{
      
      tmp_2 <- as_tibble(cbind(kegg.genes[x,],
                               tmp_1
      ))
      
      tmp_1 <- tmp_2 %>% 
        mutate(
          Modulo = case_when( Orientation == "plus" & complement == "plus" ~ (gene_end - Insertion + 1) %% 3,
                              Orientation == "minus" & complement == "minus" ~ (Insertion - gene_start + 1) %% 3
                              
          )) %>%
        mutate(
          Frame = case_when( Modulo == 0 ~ "+1 frame", 
                             Modulo == 1 ~ "+3 frame",
                             Modulo == 2 ~ "+2 frame"),
          Translation = case_when(Modulo == 0 ~ "yes",
                                  Modulo == 1 ~ "no",
                                  Modulo == 2 ~ "no")
        ) %>% select(
          -Modulo, -nt_seq_length
        )
      
    }
    
    return(tmp_1)
    
    
  })
  
  print(TN5)
  
  translation_data <- do.call(rbind, translation_data)
  
  translation_data <- rename(translation_data, IS_orientation = Orientation)
  
  # Save translation data
  write_csv(translation_data, 
            file = paste("./99_Results/2_Find_Translationframe/Translation_IS_",
                         TN5, ".csv", sep=""))
  
  
}




