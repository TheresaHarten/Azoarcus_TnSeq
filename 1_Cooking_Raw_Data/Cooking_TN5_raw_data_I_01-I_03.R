library(tibble)
library(dplyr)
library(readr)


# Working directory
setwd("set/working/directory/")


# Convert sequencing start site into insertion site
TN5_data <- list()
i <-0

for(TN5 in c("01", "02", "03")){
  
  i <- i+1
  
  TN5_forward <- read.csv(paste("./0_Data/Tn5_Ferm_I_01-I_03/Tn5_I_", 
                                TN5, "_all_plus.txt", sep=""), header = T, sep="\t")
  TN5_forward <- as_tibble(TN5_forward)
  TN5_forward[,1] <- TN5_forward[,1]+8
  TN5_forward[,"Orientation"] <- "minus"
  
  TN5_reverse <- read.csv(paste("./0_Data/Tn5_Ferm_I_01-I_03/Tn5_I_", 
                                TN5, "_all_minus.txt",sep=""), header = T, sep="\t")
  TN5_reverse <- as_tibble(TN5_reverse)
  TN5_reverse[,1] <- TN5_reverse[,1]-8
  TN5_reverse[,"Orientation"] <- "plus"
  
  TN5_new <- bind_rows(TN5_forward, TN5_reverse) 
  TN5_new <- arrange(TN5_new, Insertion)

  write_csv(TN5_new, 
            path = paste("./0_Data/Tn5_Ferm_I_01-I_03/Tn5_Ferm_I_", 
                         TN5, "_list_ins_combined.csv", sep=""))
  
  TN5_data[[i]] <- TN5_new
  
}

# Assign index to list
names(TN5_data) <- paste("I_", c(01, 02, 03), sep = "")

# Create data frame
Ins_Data1 <- do.call(rbind, TN5_data) 

tmp <- paste(Ins_Data1[[1]], Ins_Data1[[3]], sep="")
table(table(tmp))
rm(tmp)

Ins_Data2 <- Ins_Data1 %>% group_by(Insertion, Orientation) %>%
  summarise(Frequency = sum(Frequency), Orientation = Orientation)

# Unique insertion sites
Ins_Data3 <- distinct(Ins_Data2, Insertion, Orientation, .keep_all = T)

# Sort insertion sites
Ins_Data <- arrange(Ins_Data3, Insertion)

 
# Save insertion sites
write_csv(Ins_Data, 
          path ="./0_Data/Tn5_Ferm_I_01-I_03/Tn5_Ferm_I_01-I_03_list_ins_combined.csv")



