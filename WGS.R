# Install the required package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages("readr")
BiocManager::install("VariantAnnotation")

install.packages("tidyr")
install.packages("dplyr")


# Load the installed Package
library(tidyr)
library(dplyr)

library(readr)
library(bannerCommenter)

library(VariantAnnotation)
library(stringr)
library(data.table)

library(ggplot2)
##Increase max.print to 2000
options(max.print=2000)
# Reading the contents of TSV file using read_tsv() method
wgs_df<-readr::read_tsv("../Original File/covid-wgs.filter-pass.dbsnp.select_genes.sorted.uniq-coord.tsv")

#Remove Duplicate Column
wgs_df <- subset(wgs_df, select = - AF...50)  

##Rename the left Column
colnames(wgs_df)[colnames(wgs_df) == "AF...8"] <- "AF"

head(wgs_df)


#Rename Sample Names to match Metadata file
names(wgs_df) <- gsub(pattern = "BCG-",replacement = "",x = names(wgs_df))
names(wgs_df) <- gsub(pattern = "JP-",replacement = "",x = names(wgs_df))

colnames(wgs_df)
##Read Metadata
metaData <- read.csv("../WGS_metadata_2023-08-12.csv")

##Filter Metadata with the 100 Samples in VC
SampleName <- colnames(wgs_df)
SampleName <- tail(SampleName,100)
#SampleName <- gsub('[JPBCG-]','',SampleName)

metaData.Filtered <- metaData[metaData$study_number%in%SampleName,]


##Printout MetaData.Filtered
write.csv(metaData.Filtered,paste0(getwd(),"\\metaDataFiltered.csv"),row.names = FALSE)

##Load Metadata after filtering according to the Peak-illness (Manually)

metaData.Filtered <- read.csv("../WGS/metaDataFiltered.csv")





#####################################################################################
##                                  Experiment 1.1                                 ##
##                               Allele Frequency (AF)                             ##
##                                      VS                                         ##
##  AFR_AF (Allele frequency for African populations in the 1000 Genomes dataset)  ##
#####################################################################################

# Read the VCF file
vcf1.1 <- wgs_df

#Filter: Exclude rows without AF Score or AFR_AF Score.
summary(vcf1.1$AF)
#153809
vcf1.1 <- vcf1.1[!(vcf1.1$AF=="."),] 
summary(vcf1.1$AF)
#109440

vcf1.1 <- vcf1.1[!(vcf1.1$AFR_AF=="."),] 
summary(vcf1.1$AFR_AF)
#109401




#Remove Nulls
vcf1.1 <- vcf1.1[!is.na(as.numeric(vcf1.1$AF)), ]
vcf1.1 <- vcf1.1[!is.na(as.numeric(vcf1.1$AFR_AF)), ]

#Convert to Numeric
vcf1.1$AF = as.numeric(as.character(vcf1.1$AF)) 
vcf1.1$AFR_AF = as.numeric(as.character(vcf1.1$AFR_AF)) 

# Create a data frame to store results
results1.1 <- data.frame(VariantID = character(),
                      Chrom = character(),
                      Position = numeric(),
                      Gene = character(),
                      AF = numeric(),
                      Reference_AF = numeric(),
                      Difference = numeric(),
                      CI_Lower = numeric(),
                      CI_Upper = numeric(),
                      stringsAsFactors = FALSE)

# Confidence level for the interval
confidence_level <- 0.95





# Iterate through rows (variants) of the VCF
for (i in 1:nrow(vcf1.1)) {
  af_vcf <- vcf1.1$AF[i]
  af_reference <- vcf1.1$AFR_AF[i]
  
  if (!is.na(af_vcf) && !is.na(af_reference)) {
    difference <- af_vcf - af_reference
    
    # Calculate standard error
    se <- sqrt((af_vcf * (1 - af_vcf)) / nrow(vcf1.1))
    
    # Calculate z-score for the desired confidence level
    z <- qnorm((1 + confidence_level) / 2)
    
    # Calculate confidence interval
    ci_lower <- difference - z * se
    ci_upper <- difference + z * se
    
    # Add results to the data frame
    results1.1 <- rbind(results1.1, data.frame(
      VariantID = vcf1.1$ID[i],
      Chrom = vcf1.1$Chrom[i],
      Position = vcf1.1$Position[i],
      Gene = vcf1.1$Gene[i],
      AF = af_vcf,
      Reference_AF = af_reference,
      Difference = difference,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper
    ))
  }
}

# Write results to a CSV file
write.csv(results1.1, "AF.AFR_AF_results1.1.csv", row.names = FALSE)




# Function to assign unique colors to chromosomes
assign_colors <- function(chromosomes) {
  num_chromosomes <- length(unique(chromosomes))
  rainbow_palette <- rainbow(num_chromosomes)
  return(rainbow_palette[as.factor(chromosomes)])
}



# Assign unique colors to chromosomes
results1.1$Chromosome_Color <- assign_colors(results1.1$Chrom)



# Create a scatter plot of differences with Chromosome on the x-axis and custom colors

p <- ggplot(results1.1, aes(x = Chrom, y = Difference, color = Chrom)) +
  geom_point() +
  labs(title = "Scatter Plot of Differences between Variants VS 1000Genome",
       x = "Chromosome",
       y = "Difference") +
  scale_color_manual(values = unique(results1.1$Chromosome_Color)) +  # Set custom colors
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

# Save the scatter plot
ggsave("scatter_plot_chromosome_differences_1.1.png", plot = p, width = 12, height = 8)







########################################################################################
##                                  Experiment 1.2                                    ##
##                               Allele Frequency (AF)                                ##
##                                      VS                                            ##
## gnomADe_AFR_AF (Allele frequency for African populations in the gnomeADe dataset)  ##
########################################################################################

# Read the VCF file
vcf1.2 <- wgs_df



#Filter: Exclude rows without AF Score or gnomADe_AFR_AF Score.
summary(vcf1.2$AF)
#153809
vcf1.2 <- vcf1.2[!(vcf1.2$AF=="."),] 
summary(vcf1.2$AF)
#109440

vcf1.2 <- vcf1.2[!(vcf1.2$gnomADe_AFR_AF=="."),] 
summary(vcf1.2$gnomADe_AFR_AF)
#6360


#Remove Nulls
vcf1.2 <- vcf1.2[!is.na(as.numeric(vcf1.2$AF)), ]
vcf1.2 <- vcf1.2[!is.na(as.numeric(vcf1.2$gnomADe_AFR_AF)), ]

#Convert to Numeric
vcf1.2$AF = as.numeric(as.character(vcf1.2$AF)) 
vcf1.2$gnomADe_AFR_AF = as.numeric(as.character(vcf1.2$gnomADe_AFR_AF))



# Create a data frame to store results
results1.2 <- data.frame(VariantID = character(),
                         Chrom = character(),
                         Position = numeric(),
                         Gene = character(),
                         AF = numeric(),
                         Reference_gADe_AF = numeric(),
                         Difference = numeric(),
                         CI_Lower = numeric(),
                         CI_Upper = numeric(),
                         stringsAsFactors = FALSE)




# Iterate through rows (variants) of the VCF
for (i in 1:nrow(vcf1.2)) {
  af_vcf <- vcf1.2$AF[i]
  af_reference <- vcf1.2$gnomADe_AFR_AF[i]
  
  if (!is.na(af_vcf) && !is.na(af_reference)) {
    difference <- af_vcf - af_reference
    
    # Calculate standard error
    se <- sqrt((af_vcf * (1 - af_vcf)) / nrow(vcf1.2))
    
    # Calculate z-score for the desired confidence level
    z <- qnorm((1 + confidence_level) / 2)
    
    # Calculate confidence interval
    ci_lower <- difference - z * se
    ci_upper <- difference + z * se
    
    # Add results to the data frame
    results1.2 <- rbind(results1.2, data.frame(
      VariantID = vcf1.2$ID[i],
      Chrom = vcf1.2$Chrom[i],
      Position = vcf1.2$Position[i],
      Gene = vcf1.2$Gene[i],
      AF = af_vcf,
      Reference_gADe_AF = af_reference,
      Difference = difference,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper
    ))
  }
}

# Write results to a CSV file
write.csv(results1.2, "AF.gADe_AF_results_1.2.csv", row.names = FALSE)







# Assign unique colors to chromosomes
results1.2$Chromosome_Color <- assign_colors(results1.2$Chrom)



# Create a scatter plot of differences with Chromosome on the x-axis and custom colors
p <- ggplot(results1.2, aes(x = Chrom, y = Difference, color = Chrom)) +
  geom_point() +
  labs(title = "Scatter Plot of Differences between Variants VS gnomeADe",
       x = "Chromosome",
       y = "Difference") +
  scale_color_manual(values = unique(results1.2$Chromosome_Color)) +  # Set custom colors
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

# Save the scatter plot
ggsave("scatter_plot_chromosome_differences_1.2.png", plot = p, width = 12, height = 8)








########################################################################################
##                                     Experiment 2.1                                 ##
##                       SIFT (Sorting Intolerant from Tolerant)                      ##
##                                       Deleterious                                  ##
########################################################################################

# Read the VCF file
vcf2.1 <- wgs_df

#Filter1: Excluded rows without SIFT Score.
summary(vcf2.1$SIFT)    
#153809
vcf2.1 <- vcf2.1[!(vcf2.1$SIFT=="."),] 
summary(vcf2.1$SIFT)
#1828

# Extract the SIFT Score into a new Column

results2.1 <- vcf2.1
results2.1$SIFT_SCORE <- sapply(str_extract_all(results2.1$SIFT, "(?<=\\()[^)(]+(?=\\))"), paste0, collapse =",")



#Exclude all but Deleterious

results2.1 <- results2.1[results2.1$SIFT %like% c("deleterious"), ]

##Sort according the Sift Score
results2.1 <- results2.1[order(results2.1$SIFT,results2.1$SIFT_SCORE),]






    # Create a Vector with the column names
    Vector2.1 <- c("ID","Chrom","Position","Gene","AF","SIFT","SIFT_SCORE")
    results2.1 <- results2.1[,Vector2.1]  

# Write results to a CSV file
write.csv(results2.1, "SIFT_2.1.csv", row.names = FALSE)



################################TEST PLOTS##############################################


# Count the number of variants per chromosome
chromosome_counts <- table(results2.1$Chrom)

# Convert the counts to a data frame
chromosome_counts_df <- data.frame(Chromosome = names(chromosome_counts), Count = as.integer(chromosome_counts))

# Create a connected curve plot
p <- ggplot(chromosome_counts_df, aes(x = Chromosome, y = Count, group = 1)) +
  geom_line() +
  labs(title = "Chromosome SIFT Count Curve",
       x = "Chromosome",
       y = "Variant Count") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

# Save the plot
ggsave("chromosomes_SIFT1.2.png", plot = p, width = 10, height = 6)





################################END OF TEST#############################################


########################################################################################
##                                     Experiment 2.2                                 ##
##                          PolyPhen (Polymorphism Phenotyping)                       ##
##                                     Probably_Damaging                              ##
########################################################################################


# Read the VCF file
vcf2.2 <- wgs_df


#Filter: Excluded rows without PolyPhen Score.

vcf2.2 <- vcf2.2[!(vcf2.2$PolyPhen=="."),] 
summary(vcf2.2$PolyPhen)
#1775

#Extract the PolyPhen Score into a new Column

results2.2 <- vcf2.2
results2.2$PolyPhen_SCORE <- sapply(str_extract_all(results2.2$PolyPhen, "(?<=\\()[^)(]+(?=\\))"), paste0, collapse =",")



#Exclude all but Possibly_Damaging

results2.2 <- results2.2[results2.2$PolyPhen %like% c("possibly_damaging"), ]

##Sort according the PolyPhen Score
results2.2 <- results2.2[order(results2.2$PolyPhen,results2.2$PolyPhen_SCORE),]






# Create a Vector with the column names
Vector2.2 <- c("ID","Chrom","Position","Gene","AF","PolyPhen","PolyPhen_SCORE")
results2.2 <- results2.2[,Vector2.2]  

# Write results to a CSV file
write.csv(results2.2, "PolyPhen_2.2.csv", row.names = FALSE)




# Count the number of variants per chromosome
chromosome_counts <- table(results2.2$Chrom)

# Convert the counts to a data frame
chromosome_counts_df <- data.frame(Chromosome = names(chromosome_counts), Count = as.integer(chromosome_counts))

# Create a connected curve plot
p <- ggplot(chromosome_counts_df, aes(x = Chromosome, y = Count, group = 1)) +
  geom_line() +
  labs(title = "Chromosome PolyPhen Count Curve",
       x = "Chromosome",
       y = "Variant Count") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

# Save the plot
ggsave("chromosomes_PolyPhen2.2.png", plot = p, width = 10, height = 6)




########################################################################################
##                                     Experiment 2.3                                 ##
##                                    SIFT  +  PolyPhen                               ##
##                             Deleterious + Probably_Damaging                        ##
########################################################################################


# Read the VCF file
vcf2.3 <- wgs_df


#Filter: Excluded rows without SIFT or PolyPhen Score.
vcf2.3 <- vcf2.3[!(vcf2.3$SIFT=="."),] 
summary(vcf2.3$SIFT)
#1828

vcf2.3 <- vcf2.3[!(vcf2.3$PolyPhen=="."),] 
summary(vcf2.3$PolyPhen)
#1767

#Extract the PolyPhen Score into a new Column

results2.3 <- vcf2.3
results2.3$SIFT_SCORE <- sapply(str_extract_all(results2.3$SIFT, "(?<=\\()[^)(]+(?=\\))"), paste0, collapse =",")
results2.3$PolyPhen_SCORE <- sapply(str_extract_all(results2.3$PolyPhen, "(?<=\\()[^)(]+(?=\\))"), paste0, collapse =",")




#Exclude all but Deleterious

results2.3 <- results2.3[results2.3$SIFT %like% c("deleterious"), ]


#Exclude all but Possibly_Damaging

results2.3 <- results2.3[results2.3$PolyPhen %like% c("possibly_damaging"), ]

##Sort according the Sift Score
results2.3 <- results2.3[order(results2.3$SIFT,results2.3$SIFT_SCORE),]

# Create a Vector with the column names
Vector2.3 <- c("ID","Chrom","Position","Gene","AF","SIFT","SIFT_SCORE","PolyPhen","PolyPhen_SCORE")
results2.3 <- results2.3[,Vector2.3]  

# Write results to a CSV file
write.csv(results2.3, "SIFT_PolyPhen_2.3.csv", row.names = FALSE)



      ################################################################################




#####################################################################################
##                        3.1) For polyphen/sift score file                        ##
##                         - See where SNP are on genome                           ##
##   - Scan other sources for deleterious effect / status (especially if Novel)    ##
#####################################################################################

# Read the VCF file
vcf3.1 <- wgs_df


#Filter: Excluded rows without SIFT or PolyPhen Score.
vcf3.1 <- vcf3.1[!(vcf3.1$SIFT=="."),] 
summary(vcf3.1$SIFT)
#1828

vcf3.1 <- vcf3.1[!(vcf3.1$PolyPhen=="."),] 
summary(vcf3.1$PolyPhen)
#1767

#Extract the PolyPhen Score into a new Column

results3.1 <- vcf3.1
results3.1$SIFT_SCORE <- sapply(str_extract_all(results3.1$SIFT, "(?<=\\()[^)(]+(?=\\))"), paste0, collapse =",")
results3.1$PolyPhen_SCORE <- sapply(str_extract_all(results3.1$PolyPhen, "(?<=\\()[^)(]+(?=\\))"), paste0, collapse =",")




#Exclude all but Deleterious

results3.1 <- results3.1[results3.1$SIFT %like% c("deleterious"), ]


#Exclude all but Possibly_Damaging

results3.1 <- results3.1[results3.1$PolyPhen %like% c("possibly_damaging"), ]

##Sort according the Sift Score
results3.1 <- results3.1[order(results3.1$SIFT,results3.1$SIFT_SCORE),]

# Create a Vector with the column names
#Vector3.1 <- c("ID","Chrom","Position","Ref","Alt","Gene","AF","SIFT","SIFT_SCORE","PolyPhen","PolyPhen_SCORE")
#results3.1 <- results3.1[,Vector3.1]  

# Write results to a CSV file
write.csv(results3.1, "SIFT_PolyPhen_Position3.1.csv", row.names = FALSE)



# Extracting SNP positions
snps <- results3.1[nchar(results3.1$Ref) == 1 & nchar(results3.1$Alt) == 1, ]
head(snps)
snp_locations <- snps[, c("Chrom", "Position")]




install.packages("httr")
install.packages("jsonlite")
library(httr)
library(jsonlite)

#Combine the 'Chrom' and 'Position' columns into a single string that can be used in a query. 
#This typically takes the form of "chr[Chrom]:[Position]"

#Create The Query Parameters
snps$QueryP1 <- gsub("chr", "", snps$Chrom)
snps$QueryP2 <- paste0(snps$Position)


#Create a Function to Query ClinVar [Step1: Getting list of IDs]

query_clinvar <- function(chromosome, position) {
  # Construct the query URL
  base_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
  query <- paste(chromosome, "[Chromosome]+AND+", position, "[Base+Position]", sep="")
  full_url <- paste(base_url, "?db=clinvar&term=", query, "&retmode=json", sep="")
  
  # Make the API request
  response <- GET(full_url)
  
  # Check if the request was successful
  if (status_code(response) == 200) {
    # Parse the response content
    content <- fromJSON(rawToChar(response$content))
    
    # Extract the idlist part
    # idlist <- content$esearchresult$idlist
   
    # Extract the required parts
    count <- content$esearchresult$count
    idlist <- content$esearchresult$idlist
    querytranslation <- content$esearchresult$querytranslation
    
    # Return a list containing count, idlist, and querytranslation
    return(list(count = count, idlist = idlist, querytranslation = querytranslation))
    
  } else {
    warning("API request failed with status code: ", status_code(response))
    return(NULL)
  }
}

#Test API
GET("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=15[Chromosome]+AND+89804583[Base+Position]&retmode=json")

#Query Each SNP and Store the Results

#Qresults <- query_clinvar(snps$QueryP1[133],snps$QueryP2[133])



Qresults <- apply(snps, 1, function(x) query_clinvar(x['QueryP1'], x['QueryP2']))


#convert the Qresults list into dataframe

Qresults_df <- do.call(rbind, Qresults)
Qresults_df <- data.frame(Qresults_df)
# Convert row names to a column if needed
results_dataframe <- data.frame(ID = rownames(results_dataframe), results_dataframe, row.names = NULL)

# Extract the number before [Chromosome] for ex: "15[Chromosome] AND 000089804583[Base Position]"
Qresults_df$Chrom <- gsub("^([XY]|\\d+).*", "\\1", Qresults_df[,3])

# Extract the number between 'AND' and '[Base Position]' and remove leading zeros
Qresults_df$Position <- gsub(".*AND 0*(\\d+).*", "\\1", Qresults_df[,3])

str(Qresults_df)

library(tidyr)

# Assuming df is your dataframe and 'list_column' is the name of the column with lists
Qresults_df <- Qresults_df %>% unnest(count)
Qresults_df <- Qresults_df %>% unnest(querytranslation)


# Determine the maximum length of lists in the idlist column
max_length <- max(sapply(Qresults_df$idlist, length))

# Create new columns for each element up to the max length
for(i in 1:max_length) {
  Qresults_df[paste("id", i, sep = "")] <- sapply(Qresults_df$idlist, function(x) {
    if(length(x) >= i) return(x[i])
    else return(NA)
  })
}

#remove idlist
Qresults_df$idlist <- NULL




#1- Flatten the first DataFrame which contains multiple ID columns (Qresult_df)

API1_df <- Qresults_df %>%
  pivot_longer(
    cols = starts_with("id"),  # Select columns that start with 'id'
    names_to = "ID_name",       # This temporary column can be ignored or removed later
    values_to = "id"            # The consolidated ID values
  ) %>%
  filter(!is.na(id)) %>%       # Remove rows where the ID is NA
  select(-ID_name)             # Remove the temporary column







#Create a Function to Query ClinVar [Step2: Getting summary of the retirieved IDs]


query_summary <- function(id) {
  # Construct the query URL
  base_url2 <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id="
  full_url2 <- paste(base_url2, id,"&retmode=json", sep="")
  
  
  # Make the API request
  response <- GET(full_url2)
  
  # Check if the request was successful
  if (status_code(response) == 200) {
    # Parse the response content
    content <- fromJSON(rawToChar(response$content))
    
    
    # Extract the required parts
    
    # Check if the ID exists in the content$result
    if (!is.null(content$result[[id]])) {
      # Fetch the required info if the ID exists
      id<- content$result$uid
      accession <- content$result[[id]]$accession
      clin_desc <- content$result[[id]]$clinical_significance$description    } else {
      # Handle the case where the ID is NULL or doesn't exist
      id <- NA
      accession <- NA
      clin_desc <- NA  # or any other placeholder or action you deem appropriate
      
    }
    


    # Return a list containing count, idlist, and querytranslation
    return(list(id = id, accession = accession, clin_desc = clin_desc))  
    
    

    
  } else {
    warning("API request failed with status code: ", status_code(response))
    return(NULL)
  }
}


#Query Each SNP and Store the Results

Q2results_id3 <- query_summary(1279989)



Q2results <- apply(API1_df, 1, function(x) query_summary(x['id']))

#convert the Q2results list into dataframe
Q2results_df <- do.call(rbind, Q2results)
Q2results_df <- data.frame(Q2results_df)
# Assuming df is your dataframe and 'list_column' is the name of the column with lists
Q2results_df <- Q2results_df %>% unnest(accession)
Q2results_df <- Q2results_df %>% unnest(id)
Q2results_df <- Q2results_df %>% unnest(clin_desc)

API2_df <- data.frame(Q2results_df)


#Handling both API1 and API2 dataframes in one

  #Merge the dataframes using a full outer join
library(dplyr)

API_ALL <- API1_df %>%
  full_join(API2_df, by = "id")

# The full_join function keeps all rows from both dataframes.
# Columns from df1 and df2 are included with NA values if there is no match.

write.csv(API_ALL,"API_ALL.CSV",row.names=FALSE)

#Write the First result which contains nulls in ID in the first API (Qresult_df)
# Filter rows where all id1, id2, and id3 are NA or NULL
API_NULL <- Qresults_df %>%
  filter(is.na(id1) & is.na(id2) & is.na(id3))

# Write the filtered dataframe
write.csv(API_NULL,"API_NULL.CSV",row.names=FALSE)


###############################################################################################



##################################################################################################
##                                  4.1) ACE1 and ACE2 anlysis                                  ##
##   - Take ACE1 and ACE2 activity and see if it can be predictave of SNP in these two genes    ##
##################################################################################################

#Create a vcf var for the unfiltered data

vcf4.1 <- wgs_df

##Read Metadata
metaData_2810 <- read.csv("../WGS/WGS_metadata_2023-10-28.csv")


metaData_2810.Filtered <- metaData_2810[metaData_2810$study_number%in%SampleName,]


##Printout MetaData.Filtered
write.csv(metaData_2810.Filtered,paste0(getwd(),"//metaDataFiltered_2810.csv"),row.names = FALSE)



##Load Metadata after filtering according to the Peak-illness (Manually)

metaData_2810.Filtered <- read.csv("metaDataFiltered_2810.csv")


#Merge Specific Columns from the main WGS file and the metadata provided in 28-10-2023

#Filter the VCF to include ACE1 and ACE2 only
#grep ACE jonny_and_angioedema-gene-panel_biomart.tsv
#ACE2	X	15494566	15607236
#ACE	17	63477061	63498380


vcf4.1_ACE.Filtered <- vcf4.1[vcf4.1$Chrom == "chr17" & (vcf4.1$Position >= 63477061 & vcf4.1$Position <= 63498380) | vcf4.1$Chrom == "chrX" & (vcf4.1$Position >= 15494566 & vcf4.1$Position <= 15607236), ]

#Filter specific coulmn in the vcf4.1_ACE

vcf4.1_ACE.Filtered <- select(vcf4.1_ACE.Filtered, Chrom, Position, Ref, Alt, Allele, AF, all_of(SampleName))

#Merge Phentotypes from Metadata into Gentoypes in vcf1.4

#1- Merge Phenotypes in Metadata

metaData_2810.Filtered$phenotype <- paste0(metaData_2810.Filtered$ACE1_activity, ":", metaData_2810.Filtered$ACE2_activity, ":", metaData_2810.Filtered$severity_4groups_tot)





# Loop through each column (sample) in df_genotype
for (sample_id in SampleName) {
  # Find the corresponding row in df_metadata
  phenotype_data <- metaData_2810.Filtered[metaData_2810.Filtered$study_number == sample_id, "phenotype"]
  
  # Check if phenotype_data is not empty or NA
  if (!is.na(phenotype_data) && length(phenotype_data) > 0) {
    # Concatenate the genotype with the phenotype for each cell in the column
    vcf4.1_ACE.Filtered[, sample_id] <- sapply(vcf4.1_ACE.Filtered[, sample_id], function(genotype) {
      paste(genotype, phenotype_data, sep = ":")
    })
  } else {
    warning(paste("Phenotype data not found for sample:", sample_id))
  }
}


#write Genotype_Phenotype CSV file

write.csv(vcf4.1_ACE.Filtered,"Genotype_Phenotype_ACE.csv",row.names = FALSE)




##################################################################################################
##                                  4.2) ACE1 and ACE2 anlysis                                  ##
##   - Take ACE1 and ACE2 activity and see if it can be predictave of SNP in these two genes    ##
##################################################################################################

#Creating a list of Genotype_Phenotype which is filtered using 2.3 Filter (SIFT+PolyPhen)
#Result2.3 Contains 137 SNPs which confirmed as Deleterious by both Sift and PolyPhen
#vcf4.1_ACE.Filtered Contains 836 SNPs which associated to ACE and ACE2 Genes Location


# Merging df1 with df2 based on 'Chrom' and 'Position'
results4.2 <- merge(vcf4.1_ACE.Filtered, results2.3, by = c("Chrom", "Position"))




################################################################################################



            #################################################################
            ##             5) 1000 genomes and gnomAD analysis             ##
            ##                   - Filter > 0.1 difference                 ##
            #################################################################



#####################################################################################
##                                  Experiment 5.1                                 ##
##                               Allele Frequency (AF)                             ##
##                                      VS                                         ##
##  AFR_AF (Allele frequency for African populations in the 1000 Genomes dataset). ##
##.                                  Confidence Level > 0.1                        ##
#####################################################################################

#Here we want to filter the results of experiment 1.1 to include only difference > 0.1
results5.1 <- results1.1

results5.1 <- results5.1 %>%
  filter(results5.1$Difference > 0.1 | results5.1$Difference < -0.1)


# Write results to a CSV file
write.csv(results5.1, "AF.AFR_AF.Difgt0.1.results5.1.csv", row.names = FALSE)



#Draw Plot
library(ggplot2)
# Function to assign unique colors to chromosomes
assign_colors <- function(chromosomes) {
  num_chromosomes <- length(unique(chromosomes))
  rainbow_palette <- rainbow(num_chromosomes)
  return(rainbow_palette[as.factor(chromosomes)])
}



# Assign unique colors to chromosomes
results5.1$Chromosome_Color <- assign_colors(results5.1$Chrom)



# Create a scatter plot of differences with Chromosome on the x-axis and custom colors

p <- ggplot(results5.1, aes(x = Chrom, y = Difference, color = Chrom)) +
  geom_point() +
  labs(title = "Scatter Plot of Differences between Variants VS 1000Genome with Difference > 0.1",
       x = "Chromosome",
       y = "Difference") +
  scale_color_manual(values = unique(results5.1$Chromosome_Color)) +  # Set custom colors
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

# Save the scatter plot
ggsave("scatter_plot_chromosome_differencesgt0.1_5.1.png", plot = p, width = 12, height = 8)







########################################################################################
##                                  Experiment 5.2                                    ##
##                               Allele Frequency (AF)                                ##
##                                      VS                                            ##
## gnomADe_AFR_AF (Allele frequency for African populations in the gnomeADe dataset)  ##
##                                   Confidence Level > 0.1                           ##
########################################################################################


#Here we want to filter the results of experiment 1.1 to include only difference > 0.1
results5.2 <- results1.2

results5.2 <- results5.2 %>%
  filter(results5.2$Difference > 0.1 | results5.2$Difference < -0.1)


# Write results to a CSV file
write.csv(results5.2, "AF.AFR_AF.Difgt0.1.results5.2.csv", row.names = FALSE)



#Draw Plot
library(ggplot2)
# Function to assign unique colors to chromosomes
assign_colors <- function(chromosomes) {
  num_chromosomes <- length(unique(chromosomes))
  rainbow_palette <- rainbow(num_chromosomes)
  return(rainbow_palette[as.factor(chromosomes)])
}



# Assign unique colors to chromosomes
results5.2$Chromosome_Color <- assign_colors(results5.2$Chrom)



# Create a scatter plot of differences with Chromosome on the x-axis and custom colors

p <- ggplot(results5.2, aes(x = Chrom, y = Difference, color = Chrom)) +
  geom_point() +
  labs(title = "Scatter Plot of Differences between Variants VS gnomADe with Difference > 0.1",
       x = "Chromosome",
       y = "Difference") +
  scale_color_manual(values = unique(results5.2$Chromosome_Color)) +  # Set custom colors
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

# Save the scatter plot
ggsave("scatter_plot_chromosome_differencesgt0.1_5.2.png", plot = p, width = 12, height = 8)


