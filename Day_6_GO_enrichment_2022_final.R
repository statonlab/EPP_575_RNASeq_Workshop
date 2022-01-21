####################################
#####agriGO enrichment analysis#####
####################################

# Required packages
library(tibble)
library(dplyr)
library(stringr)

# Athaliana_447_Araport11.annotation_info.txt should be in the working directory (UPLOAD TO GITHUB)

## The annotation file contains multiple rows for many genes on `locusName`. We need to get rid of the extra rows before merging the annotation with our differential expression results. We subset the annotation table `anno` for those rows that do not contain duplicated locus names.

# load annotation data
anno <- read.csv("Athaliana_447_Araport11.annotation_info.txt", header=TRUE, sep="\t", stringsAsFactors = F)
head(anno)
# keep first row for each unique gene in locusName
anno <- anno[!duplicated(anno[,2]),]
head(anno)

## Format MaxRes_sig dataframe prior to merging tables

# makes first column into a column with "VALUE" as header
library(tibble)
test_MaxRes_sig <- tibble::rownames_to_column(MaxRes_sig, "VALUE")
# remove string ".Araport11.447" from entries in VALUE column for GO formatting
library(dplyr)
library(stringr)
test_MaxRes_sig <- test_MaxRes_sig %>% mutate_at("VALUE", str_replace, ".Araport11.447", "")

## We are ready to merge the tables. We need to convert the DESeq object with our results to a `data.frame` object and select the column on each file that contains the name to use for merging. Also, we need to specify to print all the rows from our results whether they are present in the second file with `all.x = TRUE`, and not to print non-matched lines in the second file with `all.y = FALSE`.

# generating dataframe in format for agriGO enrichment analysis
myRes.anno <- merge(as.data.frame(test_MaxRes_sig), anno, by.x = "VALUE", by.y = "locusName", all.x = TRUE, all.y = FALSE)

## We can  now write results to text file(s) to use for agriGO:

# write annotated results to table
write.table(myRes.anno, file="myRes_sig.txt", sep = "\t", row.names = F)

# print one gene and GO per line and write out the reference annotation file
s <- strsplit(anno$GO, split = ",")
anno.go <- data.frame(locus = rep(anno$locusName, sapply(s, length)), GO = unlist(s))
write.table(anno.go, file="go_ref_anno.txt", sep = " ", quote = F , row.names = F, col.names = F)

# print one gene and GO each line and write out my GO results to file
s <- strsplit(myRes.anno$GO, split = ",")
myRes.go <- data.frame(locus = rep(myRes.anno$VALUE, sapply(s, length)), GO = unlist(s))
write.table(myRes.go, file="go_myRes.txt", sep = " ", quote = F, row.names = F, col.names = F)

## agriGO analysis steps ##

# Navigate to agriGO to perform analysis: http://bioinfo.cau.edu.cn/agriGO/analysis.php
# Under, "1. Select analysis tool:", ensure "Singular Enrichment Analysis (SEA)" is selected.
# Select "Customized annotation" for "2. Select the species:" and paste the content of`go_myRes.txt` into the box provided.
# Select "Customized annotated reference" for "3. Select reference:" and load the reference annotation file (go_ref_anno.txt).
# Finally, click submit and explore the results.
