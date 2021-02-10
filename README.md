# Bioinformatics_Sample_Code
 
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
