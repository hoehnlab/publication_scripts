#patient-to-multifasta.R
### Given a study and patient id for cbioportal
### downloads the patient data and creates a multifasta alignment
### of just the SNV variants in those samples for use in phylogenetic work

#manage our packages
if(!"remotes"%in%installed.packages()){
  install.packages("remotes")
}
#using dev version of this package, but there is a stable version if you want
if(!"cbioportalR"%in%installed.packages()){
  remotes::install_github("karissawhiting/cbioportalR")
}
if(!"dplyr"%in%installed.packages()){
  install.packages("dplyr")
}
if(!"stringr"%in%installed.packages()){
  install.packages("stringr")
}


library(cbioportalR)
library(dplyr)
library(stringr)


#set names for download and for output aln

#cbioportal study ID for glass data
study_id <- "difg_glass"

#this patient ID in cbioportal 
patient_id <- "GLSS-HF-4F0A"

#what you want the output to be called--just named it patient id and aln.fasta
fasta_file <- paste0(patient_id,".aln.fasta")


#data download stuff for cbioportal
set_cbioportal_db("public") #set the db
#get the sample ids available for this patient
samples <- available_samples(study_id = study_id)
samples<-samples[samples$patientId==patient_id,]
sample_ids<-samples$sampleId
#print something helpful for the user
print(paste("Found", length(sample_ids), "samples"))
print(sample_ids)


#get mutation data for these samples
print("Downloading mutation data...")
mutations <- get_mutations_by_sample(
  sample_id = sample_ids,
  study_id = study_id
)

#Now we do the filtering

#Get rid of all non-snv mutations
snvs <- mutations %>%
  filter(nchar(referenceAllele) == 1 & nchar(variantAllele) == 1) %>%
  filter(referenceAllele != "-" & variantAllele != "-")


#next we make the alignment
# we need to get unique mutation sites
#  such that we consider any mutated site in any sample
unique_sites <- snvs %>%
  select(chr, startPosition, referenceAllele) %>%
  distinct() %>%
  arrange(chr, startPosition)

#get unique samples (3 for this patient, but written generically here)
unique_samples <- unique(snvs$sampleId)

# Create a matrix for the alignment
# Rows = samples + reference, Columns = mutation sites
alignment_matrix <- matrix(nrow = length(unique_samples) + 1, 
                           ncol = nrow(unique_sites))
rownames(alignment_matrix) <- c(paste0(patient_id,"-PN"), unique_samples) #using PN for paired normal

#populate with the reference alleles first
for(i in 1:nrow(unique_sites)) {
  alignment_matrix[, i] <- unique_sites$referenceAllele[i]
}

#now add mutations for each sample, if they have them 
for(i in 1:nrow(snvs)) {
  mutation <- snvs[i, ]
  sample_name <- mutation$sampleId
  #find mutation in table
  site_index <- which(unique_sites$chr == mutation$chr & 
                        unique_sites$startPosition == mutation$startPosition)
  #update matrix
  alignment_matrix[sample_name, site_index] <- mutation$variantAllele
}

#flatten table into sequence (fasta multialign)
sequences <- apply(alignment_matrix, 1, paste, collapse = "")
fasta_lines <- character()
for(name in names(sequences)) { #one day I'll do this without a loop but not today
  fasta_lines <- c(fasta_lines, paste0(">", name), sequences[name])
}
#save the output
writeLines(fasta_lines, fasta_file)

