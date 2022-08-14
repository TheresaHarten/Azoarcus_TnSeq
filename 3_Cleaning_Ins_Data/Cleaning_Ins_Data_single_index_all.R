library(ggplot2)
library(dplyr)
library(tidyverse)




# Working directory
setwd("set/working/directory/")


# Choose input data

# (VM-Malat)
singl.index <- c("Ferm_I_10", "Ferm_I_11", "Ferm_I_12")

# (SM+N aerob)
singl.index <- c("Ferm_I_01", "Ferm_I_02", "Ferm_I_03")

# (SM+N microaerob)
singl.index <- c("Ferm_I_04", "Ferm_I_05", "Ferm_I_06")

# (SM-N microaerob, NF)
singl.index <- c("Ferm_I_07", "Ferm_I_08", "Ferm_I_09")





for(i in 1:length(singl.index)) {
  
  Tn5 <- singl.index[i]
  
  # Load translation Data
  azo.transl.data <- 
    read.csv(paste("./99_Results/2_Find_Translationframe/Translation_IS_", Tn5, 
                   ".csv", sep=""))
  
  
  # Add "gene" information 
  azo.transl.data <- 
    azo.transl.data %>% 
    mutate(gene = case_when(grepl("azo_", entry) == TRUE ~ "RNA", 
                            T ~ "protein")) 
  
  
  # Laod KEGG Data (List of 4054 azo genes)
  kegg.genes <- read_csv("./0_Data/KEGG/KEGG_gene_data.csv")
  
  
  # protein coding genes
  kegg.protein.genes <- kegg.genes %>% filter(gene == "protein")
  # protein coding gens
  azo.protein.data <- azo.transl.data %>% filter(gene == "protein")
  
  
  # RNA coding genes
  kegg.RNA.genes <- kegg.genes %>% filter(gene == "RNA")
  # RNA coding genes
  azo.RNA.data <- azo.transl.data %>% filter(gene == "RNA")
  
  
  # Find non-hit genes
  non_hit_genes <- kegg.genes %>% 
    mutate(dummy = case_when(
      entry %in% unique(azo.transl.data$entry) ~ "in",
      T ~ "out"
    )
    ) %>% 
    filter(dummy == "out") %>%
    select(-dummy) %>% 
    mutate(gene = case_when(grepl("azo_", entry) == TRUE ~ "RNA", 
                            T ~ "protein"))
  
  
  # protein coding non-hit genes
  non_hit_proteins <- non_hit_genes %>% filter(gene == "protein")
  # RNA coding non-hit genes
  non_hit_RNA <- non_hit_genes %>% filter(gene == "RNA")

  
  
  #Exclude insertions with "+1 frame" (truncated protein synthesis possible)
  
  tmp_1 <- filter(azo.protein.data,
                  is.na(Frame) | 
                    Frame != "+1 frame")
  
  length(unique(tmp_1$entry))
  
  
  # Exclude insertions located in the STOP codon
  
  tmp_2 <- filter(tmp_1,
                  (complement == "plus" &
                     (Insertion != gene_end &
                        Insertion != gene_end - 1 &
                        Insertion != gene_end - 2)
                  ) |
                    (complement == "minus" & 
                       (Insertion != gene_start &
                          Insertion != gene_start + 1 &
                          Insertion != gene_start + 2 
                       )
                    )
  )
  
  length(unique(tmp_2$entry))
  
  
  # Combine DNA gene Data with RNA gene Data
  
  tmp_3 <- bind_rows(tmp_2, azo.RNA.data)
  
  length(unique(tmp_3$entry))
  length(unique(azo.RNA.data$entry))
  

  
  # Filter "sens.Ins.last.9"
  
  tmp_3a <- filter(tmp_3, 
                   (complement == "plus" & 
                      IS_orientation == "minus") |
                     (complement == "plus" &
                        IS_orientation == "plus" &
                        Insertion < gene_end -8))
  tmp_3b <- filter(tmp_3,
                   (complement == "minus" &
                      IS_orientation == "plus") |
                     (complement == "minus" &
                        IS_orientation == "minus" &
                        Insertion > gene_start +8))
  
  length(unique(tmp_3a$entry)) + length(unique(tmp_3b$entry))
  
  # Exclude anti-sense insertions in the first 9 bp of a gene
  # Filter "anti.sens.Ins.first.9"
  
  tmp_4a <- filter(tmp_3a,
                   IS_orientation == "plus" | 
                     (IS_orientation == "minus" &
                        Insertion > gene_start +8))
  tmp_4b <- filter(tmp_3b,
                   IS_orientation == "minus" |
                     (IS_orientation == "plus" &
                        Insertion < gene_end -8))                   
  
  length(unique(tmp_4a$entry)) + length(unique(tmp_4b$entry))
  
  
  # bind rows and sort entries by insertion (increasing)
  tmp_5 <- bind_rows(tmp_4a, tmp_4b)
  tmp_5 <- as_tibble(tmp_5)
  tmp_5 <- arrange(tmp_5, Insertion)
  
  length(unique(tmp_5$entry))
  

  ins.data <- tmp_5
  
  write_csv(ins.data, path = (
    paste("./99_Results/3_Cleaning_Ins_Data/Cleaned_IS_", 
          Tn5, ".csv", sep="")
  ))
  
  
  def.zero.ins.genes <- kegg.genes %>%
    mutate(dummy = case_when(
      entry %in% unique(ins.data$entry) ~ "out",
      T ~ "in"
    )
    ) %>%
    filter(dummy == "in") %>%
    select(-dummy) %>%
    mutate(ess_criterium = case_when(
      entry %in% unique(non_hit_genes$entry) ~ "non.hit", 
      T ~ "bio.irrelevant"
    )
    )
  
  length(unique(def.zero.ins.genes$entry))
  
  # Defined zero.ins. Gene speichern
  write_csv(def.zero.ins.genes, path = (
    paste("./99_Results/3_Cleaning_Ins_Data/def_zero_ins_genes_", 
          Tn5, ".csv", sep="")
  ))
  
}

  
  

