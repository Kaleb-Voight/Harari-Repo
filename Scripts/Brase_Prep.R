library(Seurat)
library(ggplot2)
library(patchwork)

setwd("/Users/KalebVoight/Desktop/Neurons")

# Read Dataset
obj <- readRDS("Gene_Neurons.rds")

# Update Dataset
obj <- UpdateSeuratObject(obj)

## Add NeuronTypes Column for inhib and excite neurons

# Get the current clusters from the Seurat object
Idents(obj) <- obj$seurat_clusters
current_idents <- Idents(obj)

# Map cluster identities to neuron types
neuron_types <- ifelse(current_idents %in% c(0, 1, 2, 6), "Excite", "Inhib")

# Add neuron types as a new column in the metadata
obj <- AddMetaData(obj, metadata = neuron_types, col.name = "NeuronType")

### Workflow

# Extract unique samples
unique_sample <- unique(obj@meta.data$Sample_ID)

# Finding percent expressed using row position in meta data
for (sample in unique_sample) {
  
  # Find row numbers for each sample ID & obj@meta.data$NeuronType == 'Inhib'
  sample_rows <- which(obj@meta.data$Sample_ID == sample)
  
  # Extract all values from gene row in matrix
  all_values <- obj@assays$RNA$counts["MAPT", sample_rows]
  
  # Extract positive values
  positive_values <- all_values[all_values > 0]
  
  # Expression: Yes or No
  positive <- ifelse(all_values > 0, "Yes", "No")
  
  # Calculate percentages
  percentages <- ((length(positive_values) / length(sample_rows)) * 100)
  
  # Assign Percent of each variable
  assign(paste0("percent_", sample), percentages)
  
  # Assign Expression
  assign(paste0("pos_", sample), positive)
}

# Adding a percent expressed / zero expression row to meta.data
Idents(obj) <- obj$Sample_ID
current_idents <- Idents(obj)

# Created vector for each sample percent
pos_percentages <- c(pos_sample1, pos_sample2, pos_sample3, pos_sample4, pos_sample5, pos_sample6, 
  pos_sample7, pos_sample8, pos_sample9, pos_sample11, pos_sample13, pos_sample14, pos_sample15, 
  pos_sample16, pos_sample17, pos_sample19, pos_sample20, pos_sample21, pos_sample22, pos_sample23, 
  pos_sample24, pos_sample25, pos_sample26, pos_sample27, pos_sample29, pos_sample30, pos_sample31, 
  pos_sample32, pos_sample34, pos_sample35, pos_sample36, pos_sample37, pos_sample38, pos_sample39, 
  pos_sample41, pos_sample42, pos_sample43, pos_sample44, pos_sample46, pos_sample47, pos_sample49, 
  pos_sample50, pos_sample51, pos_sample52, pos_sample53, pos_sample54, pos_sample55, pos_sample56, 
  pos_sample57, pos_sample58, pos_sample59, pos_sample63, pos_sample65, pos_sample66, pos_sample67, 
  pos_sample69, pos_sample70, pos_sample72, pos_sample73, pos_sample74, pos_sample76
)

# Initialize the RORB_Exp column to 0 (or you can set it to any default value you want)
obj@meta.data$MAPT_Exp_All <- 0

# Matching row names to sample names to ensure data is ordered correctly
matching_samples <- intersect(rownames(obj@meta.data), names(pos_percentages))

# Adding properly ordered rows to meta data
obj@meta.data$MAPT_Exp_All[rownames(obj@meta.data) %in% matching_samples] <- pos_percentages[matching_samples]

# Created vector for each sample percent
percentages <- c(sample1 = percent_sample1, sample2 = percent_sample2, sample3 = percent_sample3, sample4 = percent_sample4, 
  sample5 = percent_sample5, sample6 = percent_sample6, sample7 = percent_sample7, sample8 = percent_sample8, sample9 = percent_sample9, 
  sample11 = percent_sample11, sample13 = percent_sample13, sample14 = percent_sample14, sample15 = percent_sample15, 
  sample16 = percent_sample16, sample17 = percent_sample17, sample19 = percent_sample19, sample20 = percent_sample20, 
  sample21 = percent_sample21, sample22 = percent_sample22, sample23 = percent_sample23, sample24 = percent_sample24, 
  sample25 = percent_sample25, sample26 = percent_sample26, sample27 = percent_sample27, sample29 = percent_sample29, 
  sample30 = percent_sample30, sample31 = percent_sample31, sample32 = percent_sample32, sample34 = percent_sample34, 
  sample35 = percent_sample35, sample36 = percent_sample36, sample37 = percent_sample37, sample38 = percent_sample38, 
  sample39 = percent_sample39, sample41 = percent_sample41, sample42 = percent_sample42, sample43 = percent_sample43, 
  sample44 = percent_sample44, sample46 = percent_sample46, sample47 = percent_sample47, sample49 = percent_sample49, 
  sample50 = percent_sample50, sample51 = percent_sample51, sample52 = percent_sample52, sample53 = percent_sample53, 
  sample54 = percent_sample54, sample55 = percent_sample55, sample56 = percent_sample56, sample57 = percent_sample57, 
  sample58 = percent_sample58, sample59 = percent_sample59, sample63 = percent_sample63, sample65 = percent_sample65, 
  sample66 = percent_sample66, sample67 = percent_sample67, sample69 = percent_sample69, sample70 = percent_sample70, 
  sample72 = percent_sample72, sample73 = percent_sample73, sample74 = percent_sample74, sample76 = percent_sample76
)

# Create Percent Column
obj@meta.data$MAPT_Per_All <- NA

# Assign percentages based on Sample_ID
obj@meta.data$MAPT_Per_All <- percentages[obj@meta.data$Sample_ID]

# Save File
saveRDS(obj, "All_Gene_Neurons_New.rds")

## Column Remover
obj@meta.data <- obj@meta.data[, !colnames(obj@meta.data) %in% c("RORB_Exp")]

## Column Order Editor

# View the current column order
colnames(obj@meta.data)

# Create a vector of column names
new_order <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mito", "nCount_SCT", "nFeature_SCT", "Sample_ID",
               "Gender", "PS1", "PMI", "MS4", "nTREM2", "Status", "LOAD_ADAD_EOAD", "Status_classic", "ADAD", "LOAD", 
               "EOAD", "nAPOE", "AOD", "FBW", "SCT_snn_res.0.2", "SCT_snn_res.0.4", "SCT_snn_res.0.6", "SCT_snn_res.0.8", 
               "SCT_snn_res.1", "seurat_clusters", "SCT_snn_res.0.1", "SCT_snn_res.0.05", "SCT_snn_res.0.15", "NeuronType", 
               "TAFA2_Per_All", "TAFA2_Per_Excite", "TAFA2_Per_Inhib", 
               "NEGR1_Per_All", "NEGR1_Per_Excite", "NEGR1_Per_Inhib", 
               "RORB_Per_All", "RORB_Per_Excite", "RORB_Per_Inhib", 
               "TAFA2_Exp_All", "TAFA2_Exp_Excite", "TAFA2_Exp_Inhib", 
               "NEGR1_Exp_All", "NEGR1_Exp_Excite", "NEGR1_Exp_Inhib", 
               "RORB_Exp_All", "RORB_Exp_Excite", "RORB_Exp_Inhib")

# Reorder the columns in the metadata
obj@meta.data <- obj@meta.data[, new_order]
