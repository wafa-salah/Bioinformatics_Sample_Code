# Authors: Wafa Salah  and Briana Wilson (MSTP)
# Date: November 2020


#PIPELINE
  # 1) Extract file accession and target_gene from metadata
  # 2) Store each target gene with all the files associated to it (replicates)
  #    in hash_map (Dict) --> allows storing of replicates in one place
  # 3) Read in all replicate files for a given Target_gene and:					   
  #     a) Extract “ensemble_gene_ID” and “TPM” from each file				  
  #     b) Use bioMart tool to generate the hgnc_symbol for each gene_ID 		   
  #     c) Store it all in a dataframe 		  
  #     d) Average the TPMs from the different files and save it 
  #         in a separate column in that dataframe
  #     e) Create and Save it as an R-object file that can be read in later
  # 4) Read in the control data as well as the file from step (3) and compute log2fold change --> Removed for confidentiality reasons

library(biomaRt)
library(dplyr)
library(Dict)

setwd("G:/.../ShinyApp/New_ENCODE/CRISPRi/Experiments")

meta_dt <- read.csv("meta.csv") # metadata for the dataset 
#head(meta_dt)

# (1) --------------------- Editting the target_gene name formatting in metadata --------------
target <- c()
for (item in meta_dt$Targets) {
  item1 <- sub("/targets/","",item)
  target <- append(target,item1)
}

targets <- c()
for (item in target) {
  item1 <- sub("-human/","",item)
  targets <- append(targets,item1)
}

meta_dt <- subset(meta_dt, select = -c(Targets)) #remove the 'Targets' column
meta_dt$Targets <- targets # Add a Targets column with the newly formatted target_gene names
# head(meta_dt)



# (2) ------------------- Storing all replicate files for a target_gene in a hash_map --------------------------

# initialize hash_map with any <string , string> pair
#           hash_map < target_name , file accession >
repl_dict <- dict("key" = "value") #initialized map


# Populate hash_map
rownum = 1
for (gene in meta_dt$Targets) {
  if (gene %in% repl_dict$keys)  {
    value <- repl_dict[meta_dt[rownum, "Targets"]]
    new_value <- paste(value, meta_dt[rownum, "File"], sep = " ")
    repl_dict[meta_dt[rownum, "Targets"]] <- new_value
  }
  
  else {
    repl_dict[meta_dt[rownum, "Targets"]] <- meta_dt[rownum,"File"]
  }
  
  rownum = rownum + 1
}

repl_dict$remove("key")

# check if all is stored properly
for(key in repl_dict$keys){
  print(paste(key, repl_dict[key], sep = ": "))
}


# (3) --------------- Read in files, use bioMart, average the TPM for each target_gene ----------------------------

mart <- useEnsembl(biomart = "ensembl",dataset="hsapiens_gene_ensembl", mirror = "useast")

      # Boolean values that'll keep track of how many replicate files exist for a given gene
      # All target genes have at least two replicate files
hasThree <- FALSE
hasFour <- FALSE

for (key in repl_dict$keys){
  print(key)
  file1 <- substr(repl_dict[key], 1, 11)
  ctrl1 <- read.table(paste0(file1, ".tsv", sep=""))
  ctrl11 <- ctrl1[-c(1),c(1,6)]
  
  # Citation: This particular chunk of code, that submits the querry to biomaRT, is written by Briana Wilson
  BM1 <- getBM(attributes= c("hgnc_symbol","ensembl_gene_id_version"), # retrieves hgnc_symbol given Ensembl_gene_id
               filters = "ensembl_gene_id_version",
               values = ctrl11$V1 ,
               mart = mart)
  
  ctrl11 <- subset(ctrl11, V1%in%(BM1$ensembl_gene_id_version))
  colnames(BM1) <- c("symbol", "gene")
  colnames(ctrl11) <- c("gene","TPM")
  #head(ctrl11)
  cont <- inner_join(as.data.frame(ctrl11), BM1, by=c("gene"))
  cont <- subset(cont, symbol != "")
  
  
  
  file2 <- substr(repl_dict[key], 13, 23)
  ctrl2 <- read.table(paste0(file2, ".tsv", sep=""))
  ctrl22 <- ctrl2[-c(1),c(1,6)]
  
  # Citation: This particular chunk of code, that submits the querry to biomaRT, is written by Briana Wilson
  BM2 <- getBM(attributes= c("hgnc_symbol","ensembl_gene_id_version"), #retrieves hgnc_symbol given Ensembl_gene_id
               filters = "ensembl_gene_id_version",
               values = ctrl22$V1 ,
               mart = mart)
  
  ctrl22 <- subset(ctrl22, V1%in%(BM2$ensembl_gene_id_version))
  colnames(BM2) <- c("symbol", "gene")
  colnames(ctrl22) <- c("gene","TPM")
  #head(ctrl11)
  cont2 <- inner_join(as.data.frame(ctrl22), BM2, by=c("gene"))
  cont2 <- subset(cont2, symbol != "")
  
  # "try" allows code to continue executing even if there is an error
  # safeguard against the inevitable errors that will arise for 
  #   the target_genes with only two replicate files
  try({
    file3 <- substr(repl_dict[key], 25, 35)
    ctrl3 <- read.table(paste0(file3, ".tsv", sep=""))
    ctrl33 <- ctrl3[-c(1),c(1,6)]
    
    # Citation: This particular chunk of code, that submits the querry to biomaRT, is written by Briana Wilson
    BM3 <- getBM(attributes= c("hgnc_symbol","ensembl_gene_id_version"), # retrieves hgnc_symbol given Ensembl_gene_id
                 filters = "ensembl_gene_id_version",
                 values = ctrl33$V1 ,
                 mart = mart)
    
    ctrl33 <- subset(ctrl33, V1%in%(BM3$ensembl_gene_id_version))
    colnames(BM3) <- c("symbol", "gene")
    colnames(ctrl33) <- c("gene","TPM")
    #head(ctrl33)
    cont3 <- inner_join(as.data.frame(ctrl33), BM3, by=c("gene"))
    cont3 <- subset(cont3, symbol != "")
    hasThree = TRUE
  })
  
  
  try({ 
  file4 <- substr(repl_dict[key], 37, 47)
  ctrl4 <- read.table(paste0(file4, ".tsv", sep=""))
  ctrl44 <- ctrl4[-c(1),c(1,6)]
  
  # Citation: This particular chunk of code, that submits the querry to biomaRT, is written by Briana Wilson
  BM4 <- getBM(attributes= c("hgnc_symbol","ensembl_gene_id_version"), # retrieves hgnc_symbol given Ensembl_gene_id
               filters = "ensembl_gene_id_version",
               values = ctrl44$V1 ,
               mart = mart)
  ctrl44 <- subset(ctrl44, V1%in%(BM4$ensembl_gene_id_version))
  colnames(BM4) <- c("symbol", "gene")
  colnames(ctrl44) <- c("gene","TPM")
  #head(ctrl44)
  cont4 <- inner_join(as.data.frame(ctrl44), BM4, by=c("gene"))
  cont4 <- subset(cont4, symbol != "")
  hasFour = TRUE
  })
  
  
  
  if (hasFour) { # there are four replicate files, average them accordingly
      knockdowndf <- data.frame(cont$gene, as.numeric(cont$TPM), cont$symbol,cont2$gene, as.numeric(cont2$TPM), cont2$symbol , cont3$gene, as.numeric(cont3$TPM), cont3$symbol,cont4$gene, as.numeric(cont4$TPM), cont4$symbol)
      colnames(knockdowndf) <- c("gene1", "TPM1", "symbol1", "gene2", "TPM2", "symbol2" , "gene3", "TPM3", "symbol3", "gene4", "TPM4", "symbol4")
      knockdowndf$AverageTPM <- (knockdowndf$TPM1 + knockdowndf$TPM2 + knockdowndf$TPM3 + knockdowndf$TPM4)/4
      #head(knockdowndf)
      saveRDS(knockdowndf, file = paste0("Cleaned_Data2/", key, ".rds", sep = ""))}
  
  else if (hasThree) { # there are four replicate files, average them accordingly
      knockdowndf <- data.frame(cont$gene, as.numeric(cont$TPM), cont$symbol,cont2$gene, as.numeric(cont2$TPM), cont2$symbol , cont3$gene, as.numeric(cont3$TPM), cont3$symbol)
      colnames(knockdowndf) <- c("gene1", "TPM1", "symbol1", "gene2", "TPM2", "symbol2" , "gene3", "TPM3", "symbol3")
      knockdowndf$AverageTPM <- (knockdowndf$TPM1 + knockdowndf$TPM2 + knockdowndf$TPM3)/3
      #head(knockdowndf)
      saveRDS(knockdowndf, file = paste0("Cleaned_Data2/", key, ".rds", sep = ""))}
  
  else { #only two replicate files
      knockdowndf <- data.frame(cont$gene, as.numeric(cont$TPM), cont$symbol,cont2$gene, as.numeric(cont2$TPM), cont2$symbol) 
      colnames(knockdowndf) <- c("gene1", "TPM1", "symbol1", "gene2", "TPM2", "symbol2")
      knockdowndf$AverageTPM <- (knockdowndf$TPM1 + knockdowndf$TPM2)/2
      #head(knockdowndf)
      saveRDS(knockdowndf, file = paste0("Cleaned_Data2/", key, ".rds", sep = ""))
    }
  
  repl_dict$remove(key)
  hasThree = FALSE
  hasFour = FALSE
}



#(4) --------------- Compute log2fold change------------------

## This section of code has been removed for confidentiality reasons!





## bioMart Citation

#Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package
#biomaRt. Steffen Durinck, Paul T. Spellman, Ewan Birney and Wolfgang Huber, Nature Protocols
#4, 1184-1191 (2009).

#BioMart and Bioconductor: a powerful link between biological databases and microarray data
#analysis. Steffen Durinck, Yves Moreau, Arek Kasprzyk, Sean Davis, Bart De Moor, Alvis
#Brazma and Wolfgang Huber, Bioinformatics 21, 3439-3440 (2005).
