library(ggplot2)
library(dplyr)
library(tidyverse)




# Working directory
setwd("set/working/directory")

# Data Input
Tn5 <- "Library_I_13-I_19"

# Load Cleaned_IS data
ins.data.clean <- read.csv(paste("./99_Results/3_Cleaning_Ins_Data/Cleaned_IS_", 
                                 Tn5, ".csv", sep=""))



# How to define essential genes?
# Set thresholds
drop_first_ins <- FALSE
distance_1s_2nd <- 100
non_disrup_threshold <- 0.9
non_disrup_threshold_gap <- 0.8
density_threshold <- 1/500



ins.data.f1 <- lapply(unique(ins.data.clean$entry), function(x){
  
  data_tmp <- ins.data.clean %>% filter(entry == x)
  
  # only 1 insertion site
  if(nrow(data_tmp) == 1){
    
    ORF_length <- 
      as.numeric(data_tmp[1, "gene_end"] - data_tmp[1, "gene_start"] + 1)
    
    target_length <- ORF_length - 12
    
    data_out <- as_tibble(
      data.frame(
        entry=x, 
        name = data_tmp[1, "name"],
        gene = data_tmp$gene[1],
        complement = data_tmp$complement[1],
        density = nrow(data_tmp)/target_length,
        target.length = target_length,
        ins_count = nrow(data_tmp),
        non.disrupted = 1,
        drop.first.ins.criterion = drop_first_ins
      )
    )
    
    data_out <- data_out %>% 
      mutate(
        essentiality =  "essential")
    
    data_out <- data_out %>% 
      mutate(
        ess_criterium =  "single Insertion")
    
  }
  
  
  
  if(nrow(data_tmp) > 1){
    
    # looking for non disruptable 5' from start till first/second Insertion
    if(drop_first_ins == TRUE){
      
      
      if(data_tmp$complement[1] == "plus"){
        ins.1st <- sort(data_tmp$Insertion, decreasing = FALSE)[1]
        ins.2nd <- sort(data_tmp$Insertion, decreasing = FALSE)[2]
        
        if(ins.2nd - ins.1st <= distance_1s_2nd){
          dropped_first <- "no"
          non.disrupt <- ins.1st + 1 - data_tmp$gene_start[1]
        }else{
          dropped_first <- "yes"
          non.disrupt <- ins.2nd + 1 - data_tmp$gene_start[1]
        }
        
      }else{
        ins.1st <- sort(data_tmp$Insertion, decreasing = TRUE)[1]
        ins.2nd <- sort(data_tmp$Insertion, decreasing = TRUE)[2]
        
        if(ins.1st - ins.2nd <= distance_1s_2nd){
          dropped_first <- "no"
          non.disrupt <- data_tmp$gene_end[1] + 1 - ins.1st
        }else{
          dropped_first <- "yes"
          non.disrupt <- data_tmp$gene_end[1] + 1 - ins.2nd
        }
        
      }
      
    }else{  ## no dropping of first ins 
      
      dropped_first <- "no"
      if(data_tmp$complement[1] == "plus"){
        non.disrupt <- min(data_tmp$Insertion) +1 - data_tmp$gene_start[1]
      }else{
        non.disrupt <- data_tmp$gene_end[1] + 1 - max(data_tmp$Insertion)
      }
    }
    
    ORF_length <- 
      as.numeric(data_tmp[1, "gene_end"] - data_tmp[1, "gene_start"] + 1)
    
    target_length <- ORF_length - 12
  
    data_out <- as_tibble(
      data.frame(
        entry=x, 
        name = data_tmp[1, "name"],
        gene = data_tmp$gene[1],
        complement = data_tmp$complement[1],
        density = nrow(data_tmp)/target_length,
        target.length = target_length,
        ins_count = nrow(data_tmp),
        non.disrupted = non.disrupt/target_length,
        drop.first.ins.criterion = drop_first_ins
      )
    )
    
    data_out <- data_out %>% 
      mutate(
        essentiality = case_when(
          non.disrupted >= non_disrup_threshold ~ "essential",
          T ~ "non.essential"
        )
      )
    
    if(data_out$essentiality == "essential"){
      data_out <- data_out %>% 
        mutate(
          ess_criterium =  "5' region")
    }else{
      data_out <- data_out %>% 
        mutate(
          ess_criterium =  NA)
    }
    
    
    if(data_out$essentiality == "non.essential"){
      Ins_diffs <- abs(diff(sort(data_tmp$Insertion, decreasing = FALSE)))
      max_Ins_diffs <- max(Ins_diffs)
      
      if(max_Ins_diffs/target_length >= non_disrup_threshold_gap & 
         data_out$density < density_threshold){
        data_out$essentiality <- "essential"
        data_out$ess_criterium <- "gap"
      }
      
    }
    
  }
  
  
  data_out
  
})

ins.data.f1 <- do.call(bind_rows, ins.data.f1)

# Save ins.data.f1
write_csv(ins.data.f1, path = (
  paste("./99_Results/4_Find_Essential_Genes/Ins_Data_f1_", 
        Tn5, ".csv", sep="")
))


ins.data.f1.ess <- ins.data.f1 %>% 
  filter(essentiality == "essential")


def.zero.ins.genes <- 
  read.csv(paste("./99_Results/3_Cleaning_Ins_Data/def_zero_ins_genes_", Tn5, 
                 ".csv", sep=""))



# combine ins.data.f1.ess and def.zero.ins.genes 
ess.genes.tmp1 <- bind_rows(ins.data.f1.ess, def.zero.ins.genes)



# Apply non.disrupt. 5' region criterion for ess.genes.tmp1, too

singl.ins.genes <- lapply(unique(ins.data.clean$entry), function(x){
  
  data_out <- NULL
  
  data_1 <- ins.data.clean %>% filter(entry == x)
  
  if(nrow(data_1) == 1){
    
    if(data_1$complement[1] == "plus"){
      
      non.disrupt.5 <- 
        as.numeric(data_1[1, "Insertion"] - data_1[1, "gene_start"] + 1)
      
      ORF_length <- 
        as.numeric(data_1[1, "gene_end"] - data_1[1, "gene_start"] + 1)
      
      target_length <- ORF_length - 12
      
    }else{
      non.disrupt.5 <-
        as.numeric(data_1[1, "gene_end"] - data_1[1, "Insertion"] + 1)
      
      ORF_length <- 
        as.numeric(data_1[1, "gene_end"] - data_1[1, "gene_start"] + 1)
      
      target_length <- ORF_length - 12
      
    }
    
    data_out <- as_tibble(
      data.frame(
        entry=x, 
        name = data_1[1, "name"],
        gene = data_1$gene[1],
        complement = data_1$complement[1],
        density = 1/target_length,
        target.length = target_length,
        ins_count = 1,
        non.disrupt.5 = non.disrupt.5/target_length
      )
    )
  }
  data_out
})


singl.ins.genes <- do.call(bind_rows, singl.ins.genes)



singl.ins.genes_out <- singl.ins.genes %>% filter(non.disrupt.5 <= non_disrup_threshold)


ess.genes.tmp2 <- ess.genes.tmp1 %>% 
  mutate(dummy = case_when(
    entry %in% singl.ins.genes_out$entry ~ "out",
    T ~ "in"
  )
  ) %>% 
  filter(dummy == "in") %>%
  select(-dummy)


# Laden der kegg.genes (Gen Liste mit allen 4054 Genen)
kegg.genes <- 
  read_csv("./0_Data/KEGG/KEGG_gene_data.csv")


kegg.ess.genes.tmp1 <- kegg.genes %>%
  mutate(dummy = case_when(
    entry %in% ess.genes.tmp2$entry ~ "in",
    T ~ "out"
  )) %>%
  filter(dummy == "in") %>%
  select(-dummy)


kegg.ess.genes.tmp2 <- kegg.ess.genes.tmp1 %>% 
  mutate(ins_count = case_when(
    entry %in% ins.data.f1$entry ~ "TRUE",
    T ~ "FALSE"
  ))


kegg.ess.genes.tmp2$ins_count <- as.logical(kegg.ess.genes.tmp2$ins_count)


kegg.ess.genes.tmp3 <- lapply(kegg.ess.genes.tmp2$entry, function(x){

  data_out <- NULL
  
  data_tmp <- kegg.ess.genes.tmp2 %>% filter(entry == x)
  
  if(data_tmp$ins_count == FALSE){
    data_tmp$ins_count <- as.numeric(0)
  } else {
    data_tmp2 <- ins.data.f1.ess %>% filter(entry == x)
    data_tmp$ins_count <- as.numeric(data_tmp2$ins_count)
  }
  data_out <- as_tibble(
    data.frame(
      entry= x, 
      name = kegg.ess.genes.tmp2$name[kegg.ess.genes.tmp2$entry==x],
      definition = kegg.ess.genes.tmp2$definition[kegg.ess.genes.tmp2$entry==x],
      DNA = kegg.ess.genes.tmp2$complement[kegg.ess.genes.tmp2$entry==x],
      start = kegg.ess.genes.tmp2$gene_start[kegg.ess.genes.tmp2$entry==x],
      end = kegg.ess.genes.tmp2$gene_end[kegg.ess.genes.tmp2$entry==x],
      length = kegg.ess.genes.tmp2$nt_seq_length[kegg.ess.genes.tmp2$entry==x],
      gene = kegg.ess.genes.tmp2$gene[kegg.ess.genes.tmp2$entry==x],
      ins_count = data_tmp$ins_count[1],
      target_length = ((data_tmp$gene_end[1] - data_tmp$gene_start[1])+1)-12,
      density = 1/(((data_tmp$gene_end[1] - data_tmp$gene_start[1])+1)-12),
      ess_criterium = ess.genes.tmp2$ess_criterium[ess.genes.tmp2$entry==x]
    ))
  
  data_out
})

ess.genes <- do.call(rbind, kegg.ess.genes.tmp3)



write_csv(ess.genes, path = (
  paste("./99_Results/4_Find_Essential_Genes/Ess_Genes_", 
        Tn5, ".csv", sep="")
))


