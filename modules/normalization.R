#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## Script Preparation ##
# Only edit script where instructed to do so. Indicated by EDIT
##

# Read in raw counts from CSV file.
# EDIT: set path to raw counts file below
counts <- data.frame(read.table(args[1],
                                header = TRUE,
                                sep = ",",
                                row.names = 1))

# Capture Inputs (the string "input" must appear somewhere in the sample name to be identified)
inputs <- counts[, grep(".*input*", colnames(counts), ignore.case = TRUE)]

# Rename Column Names to letters A,B,C, etc.
colnames(inputs) <- LETTERS[1:length(inputs)]

# Remove inputs from counts variable
counts[, grep(".*input*", colnames(counts), ignore.case = TRUE)] <- NULL

# Include required Libraries
library(ggplot2)

pdf(file = paste(args[2], "/Rplots.pdf", sep=""))

#----------------------------------------------------------------------------------

# Normalize all inputs
norm_inputs <- data.frame(matrix(ncol = ncol(inputs), nrow = nrow(inputs)),
                          row.names = rownames(inputs))
colnames(norm_inputs) <- colnames(inputs)

for(n in colnames(inputs))
{
  norm_inputs[n] <- inputs[n] / sum(inputs[n])
}

# Input Replicate Correlation Tests
panel.cor <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(x, y, method = "pearson")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "")
  text(0.5, 0.6, txt)
  # p-value calculation
  p <- cor.test(x, y)$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  txt2 <- paste("p= ", txt2, sep = "")
  if(p<0.01)
  {
    txt2 <- paste("p= ", "<0.01", sep = "")
    text(0.5, 0.4, txt2)
  }
  else
  {
    text(0.5, 0.4, txt2, col = "red")
  }
}
pairs(norm_inputs, upper.panel = panel.cor)

#----------------------------------------------------------------------------------

# EDIT: Remove any replicates not passing correlation test.
#norm_inputs$A <- NULL

# Average input replicates
norm_input_avg <- data.frame(rowSums(norm_inputs) / ncol(norm_inputs),
                             row.names = rownames(inputs))
colnames(norm_input_avg) <- "avg"

ggplot(data = norm_input_avg,aes(x = reorder(rownames(norm_input_avg), -avg), y = avg)) +
  geom_bar(stat = "identity", color = "blue4", fill = "white") + 
  theme(axis.text.x = element_text(angle = 90, size = 6)) +
  labs(title = "Normalized Input Barcode Distribution" ,x = "Sample Names", y = "Normalized Counts")

#----------------------------------------------------------------------------------

#Only average values that are non-zero as to not skew mean
norm_mean_input <- mean(norm_input_avg[norm_input_avg != 0,1])

# Generate bar plot showing barcode distribution
ggplot(data = norm_input_avg,aes(x = reorder(rownames(norm_input_avg), -avg), y = avg)) +
  geom_bar(stat = "identity", color = "blue4", fill = "white") + 
  theme(axis.text.x = element_text(angle = 90, size = 6), legend.position = c(0.85,0.9)) +
  labs(title = "Normalized Input Barcode Distribution" ,x = "Sample Names", y = "Normalized Counts") +
  geom_hline(aes(yintercept = norm_mean_input/5, linetype = "5x"), color = "red", show.legend = TRUE) + 
  geom_hline(aes(yintercept = norm_mean_input/10, linetype = "10x"), color = "green", show.legend = TRUE) +
  scale_linetype_manual(name = "Threshold", values = c(1,1), guide = guide_legend(override.aes = list(color = c("green", "red"))))

#----------------------------------------------------------------------------------

# Create list of BCs under threshold
names_under_thresh <- rownames(subset(norm_input_avg, norm_input_avg < (norm_mean_input/10)))
cat("List of Removed Barcodes\n", names_under_thresh)

# Remove these barcodes from our counts list 
counts <- counts[! rownames(counts) %in% names_under_thresh, ]
inputs <- inputs[! rownames(inputs) %in% names_under_thresh, ]

#----------------------------------------------------------------------------------

#Re-normalize all inputs
norm_inputs <- data.frame(matrix(ncol = ncol(inputs), nrow = nrow(inputs)),
                          row.names = rownames(inputs))
colnames(norm_inputs) <- colnames(inputs)

for(n in colnames(inputs))
{
  norm_inputs[n] <- inputs[n] / sum(inputs[n])
}

# Average input replicates
norm_input_avg <- data.frame(rowSums(norm_inputs) / ncol(norm_inputs),
                             row.names = rownames(inputs))
colnames(norm_input_avg) <- "avg"

#----------------------------------------------------------------------------------

#First, remove samples whose columns sum to zero.
cnt0 <- counts[,-1]
rownames(cnt0) <- counts[,1]
counts <- cnt0
counts_nz <- counts[colSums(counts) != 0]
#List samples removed
if(length(names(counts[colSums(counts) == 0])) != 0) {
  cat("Samples removed for zero counts\n", names(counts[colSums(counts) == 0]))
} else {
  cat("All Samples Retained")
}

#----------------------------------------------------------------------------------

#Normalize each sample to unity
selfnorm_counts <- sweep(counts_nz, 2, colSums(counts_nz), '/')

#Normalize counts to normalized input average
inputnorm_counts <- selfnorm_counts / as.matrix(norm_input_avg)

#Scale counts so each column totals 100
scale_factor <- 100 / colSums(inputnorm_counts)
final_norm_counts <- as.matrix(inputnorm_counts) %*% diag(scale_factor)
dimnames(final_norm_counts) <- list(rownames(inputnorm_counts), colnames(inputnorm_counts))

#Write final normalized to counts to output CSV file
write.csv(file = paste(args[2], "/normcounts.csv", sep=""), final_norm_counts) 

