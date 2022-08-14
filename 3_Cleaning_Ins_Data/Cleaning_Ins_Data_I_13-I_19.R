library(ggplot2)
library(dplyr)
library(tidyverse)



# Working directory
setwd("set/working/directory/")

# Definiere Data Input MAster Library
Tn5 <- "Library_I_13-I_19"

# Load Translation Data
azo.transl.data <- 
  read.csv(paste("./99_Results/2_Find_Translationframe/Translation_IS_", Tn5, 
                 ".csv", sep=""))


# Add gene info
azo.transl.data <- 
  azo.transl.data %>% 
  mutate(gene = case_when(grepl("azo_", entry) == TRUE ~ "RNA", 
                          T ~ "protein")) 



# Load KEGG Data
kegg.genes <- 
  readRDS("./0_Data/KEGG_azo_all_genes_data_frame.rds")

# Add "gene" Info 
kegg.genes <- kegg.genes %>% 
  mutate(gene = case_when(grepl("azo_", entry) == TRUE ~ "RNA", 
                          T ~ "protein"))

# Save kegg.genes
write_csv(kegg.genes, path = ("./0_Data/KEGG/KEGG_gene_data.csv"))


# Filter protein coding genes
kegg.protein.genes <- kegg.genes %>% filter(gene == "protein")
# Filter protein coding genes
azo.protein.data <- azo.transl.data %>% filter(gene == "protein")


# Filter RNA coding genes
kegg.RNA.genes <- kegg.genes %>% filter(gene == "RNA")
# Filter RNA coding genes
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


# Filter non_hit_genes for protein coding genes
non_hit_proteins <- non_hit_genes %>% filter(gene == "protein")


# Filter non_hit_genes for RNA coding genes
non_hit_RNA <- non_hit_genes %>% filter(gene == "RNA")



#Exclude insertions with "+1 frame" (truncated protein synthesis possible)

table(azo.protein.data$Frame)

tmp_1 <- filter(azo.protein.data,
                is.na(Frame) | 
                  Frame != "+1 frame")


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


# Combine DNA gene Data with RNA gene Data

tmp_3 <- bind_rows(tmp_2, azo.RNA.data)


# Exclude sense insertions in the last 9 bp of a gene

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


 
# Exclude anti-sense insertions in the first 9 bp of a gene 

tmp_4a <- filter(tmp_3a,
                 IS_orientation == "plus" | 
                   (IS_orientation == "minus" &
                      Insertion > gene_start +8))
tmp_4b <- filter(tmp_3b,
                 IS_orientation == "minus" |
                   (IS_orientation == "plus" &
                      Insertion < gene_end -8))                   



# bind rows and sort entries by insertion (increasing)


tmp_5 <- bind_rows(tmp_4a, tmp_4b)
tmp_5 <- as_tibble(tmp_5)
tmp_5 <- arrange(tmp_5, Insertion)


# Insertion Data cleaned 
Ins_Data_clean <- tmp_5

write_csv(Ins_Data_clean, 
  paste("./99_Results/3_Cleaning_Ins_Data/Cleaned_IS_", 
        Tn5, ".csv", sep="")
  )



def.zero.ins.genes <- kegg.genes %>%
  mutate(dummy = case_when(
    entry %in% unique(Ins_Data_clean$entry) ~ "out",
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


# Defined zero.ins. Gene speichern
write_csv(def.zero.ins.genes, path = (
  paste("./Results/3_Cleaning_Ins_Data/def_zero_ins_genes_", 
        Tn5, ".csv", sep="")
  ))


