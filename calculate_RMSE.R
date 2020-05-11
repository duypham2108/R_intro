#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# RMSE function
RMSE <- function(ground_truth, student_result){
  d <- c(as.matrix(ground_truth)) - c(as.matrix(student_result))
  mask <- !is.na(d)
  n_valid <- sum(mask)
  rmse <- sqrt(sum(d[mask]^2)/n_valid)
  return(rmse)
}

# Load ground truth
ground_truth <- read.table("../../Downloads//evaluating_gold_proteome_10005.txt",header=TRUE)
# Sorting gene name
ground_truth <- ground_truth[order(ground_truth$Gene_ID),]
# Remove Gene_ID
ground_truth <- ground_truth[,2:length(ground_truth)]
# Sorting sample name
column_names <- c(sort(names(ground_truth)))
ground_truth <- ground_truth[,column_names]



result = RMSE(ground_truth,ground_truth)

list_students <- list.dirs("student_submission/")
list_students <- list_students[2:length(list_students)]

RMSE_score_df <- data.frame(Student = character(),
                           Single_RMSE = numeric(),
                           Multiple_RMSE = numeric())

for (i in 1: length(list_students)){
  Student <- as.character(basename(list_students[i]))
  
  single <- read.table(paste0(list_students[i],"/Single.txt"))
  Single_RMSE <- RMSE(ground_truth,single)
  
  multiple <- read.table(paste0(list_students[i],"/Multiple.txt"))
  Multiple_RMSE <- RMSE(ground_truth,multiple)
  
  RMSE_score_df <- rbind(RMSE_score_df,data.frame(Student,Single_RMSE,Multiple_RMSE))
}

write.csv(RMSE_score_df,"RMSE_results.csv",row.names = F)
print("DONE!")
