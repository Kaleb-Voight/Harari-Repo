library(Seurat)
library(ggplot2)
library(patchwork)
library(sctransform)
library(dplyr)
library(ggpubr)
library(broom)
library(igraph)

setwd("/Users/KalebVoight/Desktop/Datasets_Scripts")

# Read Dataset
obj <- readRDS("Brase_Neurons.rds")

# Update Dataset
obj <- UpdateSeuratObject(obj)

# Violin Plot
Idents(obj) <- obj$Status

# Grabbing Unique Values from each column
unique_data <- distinct(select(obj@meta.data, Sample_ID, Status, 
                               RORB_Per_All, NEGR1_Per_All, TAFA2_Per_All, 
                               RORB_Per_Inhib, NEGR1_Per_Inhib, TAFA2_Per_Inhib, 
                               RORB_Per_Excite, NEGR1_Per_Excite, TAFA2_Per_Excite))


# Defining groups for loop
groups <- c("RORB_Per_All", "NEGR1_Per_All", "TAFA2_Per_All", 
            "RORB_Per_Inhib", "NEGR1_Per_Inhib", "TAFA2_Per_Inhib", 
            "RORB_Per_Excite", "NEGR1_Per_Excite", "TAFA2_Per_Excite")

## Violin Plot Loop
for (group in groups) {
  p <- ggplot(data = unique_data, aes(x = Status, y = get(group), fill = Status)) + 
    geom_violin(trim = TRUE) + 
    geom_jitter(width = 0.2, alpha = 0.6, size = 1, color = "black") +
    labs(title = paste(group, "Plot"),  
         x = "Status",
         y = group) + 
    scale_y_continuous(limits = c(0, 100)) +
    theme_light()
  
  print(p)
}

# Defining disease groups
disease_groups <- c("Neuro_CO", "Neuro_AD", "Neuro_ADAD", "Neuro_OT", "Neuro_Presympt")

# Defining gene names
genes <- c("RORB", "NEGR1", "TAFA2")

# Defining cell states
states <- c("All", "Inhib", "Excite")

## T-tests
run_t_test <- function(unique_data, genes, states, disease_groups) {
  # Create empty data frame to store results
  t_test_results <- data.frame(
    Gene = character(),
    Cell_State = character(),
    Disease_Group_Comparison = character(),
    t_statistic = numeric(),
    df = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop through each gene and each cell state
  for (gene in genes) {
    for (state in states) {
      
      # Define the control and disease-specific columns dynamically
      control_col <- paste0(gene, "_Per_", state)
      control_data <- unique_data[[control_col]][unique_data$Status == "Neuro_CO"]
      
      # Pull out data for each disease group
      disease_data <- lapply(disease_groups, function(status) {
        unique_data[[paste0(gene, "_Per_", state)]][unique_data$Status == status]
      })
      
      # Perform T-test for each disease group against the control
      for (i in seq_along(disease_groups)) {
        group_data <- disease_data[[i]]
        group_name <- disease_groups[i]
        
        # Perform the t-test
        t_test <- t.test(control_data, group_data)
        
        # Store the results in the results data frame
        t_test_results <- rbind(t_test_results, data.frame(
          Gene = gene,
          Cell_State = state,
          Disease_Group_Comparison = paste("Neuro_CO vs", group_name),
          t_statistic = t_test$statistic,
          df = t_test$parameter,
          p_value = t_test$p.value
        ))
      }
    }
  }
  
  # Return the results data frame
  return(t_test_results)
}

# Run function
t_test_results <- run_t_test(unique_data, genes, states, disease_groups)

## Chi-squared Analysis

# Grabbing Values from each column
chi_data <- select(obj@meta.data, Sample_ID, AOD, Gender, Status, NeuronType, seurat_clusters, 
                   RORB_Exp_All, NEGR1_Exp_All, TAFA2_Exp_All, 
                   RORB_Exp_Inhib, NEGR1_Exp_Inhib, TAFA2_Exp_Inhib, 
                   RORB_Exp_Excite, NEGR1_Exp_Excite, TAFA2_Exp_Excite)

# Expression Data
exp_columns <- c("RORB_Exp_All", "NEGR1_Exp_All", "TAFA2_Exp_All", 
                 "RORB_Exp_Inhib", "NEGR1_Exp_Inhib", "TAFA2_Exp_Inhib", 
                 "RORB_Exp_Excite", "NEGR1_Exp_Excite", "TAFA2_Exp_Excite")

# Disease groups
disease_groups <- c("Neuro_ADAD", "Neuro_AD", "Neuro_OT", "Neuro_Presympt")

# Define the function to run the Chi-squared test for Neuro_CO vs other disease groups
run_chi_squared_test <- function(chi_data, exp_columns, disease_groups) {
  
  # Initialize an empty data frame to store results
  chi_results <- data.frame(
    Expression = character(),  
    Disease_Group_Comparison = character(),  # Store the comparison groups
    Chi_squared = numeric(),   
    df = numeric(),            
    p_value = numeric(),       
    stringsAsFactors = FALSE
  )
  
  # Loop through each expression column and run the Chi-squared test
  for (col in exp_columns) {
    
    # Loop through each disease group and compare with "Neuro_CO"
    for (group in disease_groups) {
      
      # Subset the data for "Neuro_CO" and the current disease group
      subset_data <- chi_data[chi_data$Status %in% c("Neuro_CO", group), ]
      
      # Create chi table
      chi_table <- table(subset_data[[col]], subset_data$Status)
      
      # Format Table 
      chi_table <- chi_table[c("Yes", "No"), c("Neuro_CO", group)]
      
      # Print Table 
      print(chi_table)
      
      # Run the Chi-squared test
      chi_res <- chisq.test(chi_table)
      
      # Store the results in the results data frame
      chi_results <- rbind(chi_results, data.frame(
        Expression = col,
        Disease_Group_Comparison = paste("Neuro_CO vs", group),
        Chi_squared = chi_res$statistic,
        df = chi_res$parameter,
        p_value = chi_res$p.value
      ))
    }
  }
  
  # Return the results data frame
  return(chi_results)
}

# Run function
chi_results <- run_chi_squared_test(chi_data, exp_columns, disease_groups)

## Fisher Test

# Grabbing Values from each column
fish_data <- select(obj@meta.data, Sample_ID, Status, NeuronType, 
                   RORB_Exp_All, NEGR1_Exp_All, TAFA2_Exp_All, 
                   RORB_Exp_Inhib, NEGR1_Exp_Inhib, TAFA2_Exp_Inhib, 
                   RORB_Exp_Excite, NEGR1_Exp_Excite, TAFA2_Exp_Excite)

# Define the function to run the Chi-squared test for Neuro_CO vs other disease groups
run_fish_test <- function(fish_data, exp_columns, disease_groups) {
  
  # Initialize an empty data frame to store results
  fish_results <- data.frame(
    Expression = character(),  
    Disease_Group_Comparison = character(),  # Store the comparison groups
    P_Value = numeric(),   
    Odds_Ratio = numeric(),            
    Null_Odds_Ratio = numeric(),       
    stringsAsFactors = FALSE
  )
  
  # Loop through each expression column and run the Chi-squared test
  for (col in exp_columns) {
    
    # Loop through each disease group and compare with "Neuro_CO"
    for (group in disease_groups) {
      
      # Subset the data for "Neuro_CO" and the current disease group
      subset_data <- chi_data[fish_data$Status %in% c("Neuro_CO", group), ]
      
      # Create chi table
      fish_table <- table(subset_data[[col]], subset_data$Status)
      
      # Format Table 
      fish_table <- fish_table[c("Yes", "No"), c("Neuro_CO", group)]
      
      # Print Table 
      print(fish_table)
      
      # Run the Chi-squared test
      fish_res <- fisher.test(fish_table)
      
      # Store the results in the results data frame
      fish_results <- rbind(fish_results, data.frame(
        Expression = col,
        Disease_Group_Comparison = paste("Neuro_CO vs", group),
        P_Value = fish_res$p.value,
        Odds_Ratio = fish_res$estimate,
        Null_Odds_Ratio = fish_res$null.value
      ))
    }
  }
  
  # Return the results data frame
  return(fish_results)
}

# Run function
fish_results <- run_fish_test(fish_data, exp_columns, disease_groups)

## Chi-squared: Cell State RORB

# Grabbing Values from each column
chi_data_state <- obj@meta.data %>%
  filter(seurat_clusters %in% c('0', '1')) %>%
  select(Sample_ID, Status, seurat_clusters, RORB_Exp_All)

# Cell States
states <- c("0", "1")

# Define the function to run the Chi-squared test for Neuro_CO vs other disease groups
run_chi_squared_test <- function(chi_data_state, states, disease_groups) {
  
  # Initialize an empty data frame to store results
  chi_results_state <- data.frame(
    Cell_State = character(),  
    Disease_Group_Comparison = character(),  # Store the comparison groups
    Chi_squared = numeric(),   
    df = numeric(),            
    p_value = numeric(),       
    stringsAsFactors = FALSE
  )
  
  # Loop through each state (seurat_cluster 0 and 1)
  for (state in states) {
    
    # Loop through each disease group and compare with "Neuro_CO"
    for (group in disease_groups) {
      
      # Subset the data for "Neuro_CO" and the current disease group for the current state
      subset_data <- chi_data_state %>%
        filter(Status %in% c("Neuro_CO", group) & seurat_clusters == state)
      
      # Create chi table
      chi_table <- table(subset_data$RORB_Exp_All, subset_data$Status)
      
      # Format Table 
      chi_table <- chi_table[c("Yes", "No"), c("Neuro_CO", group)]
      
      # Print Table 
      print(chi_table)
      
      # Run the Chi-squared test
      chi_res <- chisq.test(chi_table)
      
      # Store the results in the results data frame
      chi_results_state <- rbind(chi_results_state, data.frame(
        Cell_State = state,
        Disease_Group_Comparison = paste("Neuro_CO vs", group),
        Chi_squared = chi_res$statistic,
        df = chi_res$parameter,
        p_value = chi_res$p.value
      ))
    }
  }
  
  # Return the results data frame
  return(chi_results_state)
}

# Run function
chi_results_state <- run_chi_squared_test(chi_data_state, states, disease_groups)

## Fisher: Cell State RORB

# Grabbing Values from each column
fish_data_state <- obj@meta.data %>%
  filter(seurat_clusters %in% c('0', '1')) %>%
  select(Sample_ID, Status, seurat_clusters, RORB_Exp_All)

# Define the function to run the Chi-squared test for Neuro_CO vs other disease groups
run_fish_test <- function(fish_data_state, states, disease_groups) {
  
  # Initialize an empty data frame to store results
  fish_results_state <- data.frame(
    Cell_State = character(),  
    Disease_Group_Comparison = character(),  # Store the comparison groups
    P_Values = numeric(),   
    Odds_Ratio = numeric(),            
    Null_Odds_Ratio = numeric(),       
    stringsAsFactors = FALSE
  )
  
  # Loop through each state (seurat_cluster 0 and 1)
  for (state in states) {
    
    # Loop through each disease group and compare with "Neuro_CO"
    for (group in disease_groups) {
      
      # Subset the data for "Neuro_CO" and the current disease group for the current state
      subset_data <- fish_data_state %>%
        filter(Status %in% c("Neuro_CO", group) & seurat_clusters == state)
      
      # Create chi table
      fish_table <- table(subset_data$RORB_Exp_All, subset_data$Status)
      
      # Format Table 
      fish_table <- fish_table[c("Yes", "No"), c("Neuro_CO", group)]
      
      # Print Table 
      print(fish_table)
      
      # Run the Chi-squared test
      fish_res <- fisher.test(fish_table)
      
      # Store the results in the results data frame
      fish_results_state <- rbind(fish_results_state, data.frame(
        Cell_State = state,
        Disease_Group_Comparison = paste("Neuro_CO vs", group),
        P_Value = fish_res$p.value,
        Odds_Ratio = fish_res$estimate,
        Null_Odds_Ratio = fish_res$null.value
      ))
    }
  }
  
  # Return the results data frame
  return(fish_results_state)
}

# Run function
fish_results_state <- run_fish_test(fish_data_state, states, disease_groups)

## Linear Model

# Create an empty list to store data frames
model_results <- list()

# Loop through disease groups ("Neuro_AD", "Neuro_CO", etc.)
for (group in disease_groups) {
  # Loop throught cell states (0, 1)
  for (state in states) {
    # Create the variable name for data storage
    var_name <- paste(group, state, sep = "_")
    
    # Filter data and prepare the dataset for modeling
    df <- chi_data %>% 
      filter(seurat_clusters == state) %>%  # Filter by cell state (0 or 1)
      filter(Status %in% c(group, "Neuro_CO")) %>%  # Keep only rows for current disease group and control
      select(RORB_Exp_All, Status, AOD, Gender) %>%  # Select relevant columns
      mutate(
        Gender = factor(Gender), # Treat gender as a factor
        RORB_Exp_All = ifelse(RORB_Exp_All == "Yes", 1, 0),  # Convert RORB expression to binary (1 for "Yes", 0 for "No")
        AOD = scale(as.numeric(as.character(AOD)))[, 1],  # Convert AOD to numeric and scale it
        Status = factor(Status, levels = c("Neuro_CO", group))  # Set factor levels for Status
      )
    
    # Logistic regression model to compare the likelihood of RORB expression across groups
    form <- as.formula("RORB_Exp_All ~ Status + AOD + Gender")  # RORB expression as a function of disease group and age
    
    # Fit logistic regression model (binomial logistic regression)
    model_summary <- glm(form, data = df, family = binomial(link = "logit"))
    
    # Get the results in a tidy format
    res <- broom::tidy(model_summary)
    
    # Store the model results in a list with the variable name
    model_results[[var_name]] <- res
  }
}

# Loop through all results
for (var_name in names(model_results)) {
  # Get the model results for the current disease group and state
  res <- model_results[[var_name]]
  
  # Apply the Benjamini-Hochberg correction to the p-values
  res$p.adjusted <- p.adjust(res$p.value, method = "BH")
  
  # Optionally, you can store the adjusted results back into the list
  model_results[[var_name]] <- res
}

FeaturePlot(obj, features = c("RORB", "FOXP2", "POU6F2", "PCDH20"), cols = c("lightgrey", "blue"))

Idents(obj) <- obj$seurat_clusters
VlnPlot(obj, features = c("RORB", "FAM19A1"))

## Network Analysis

# Compute the nearest neighbors (you can adjust k, which is the number of neighbors)
obj <- FindNeighbors(obj, dims = 1:20)  # Use the first 20 principal components as an example

# Extract the nearest neighbor matrix
nn_graph <- obj@graphs$SCT_nn

# Fetch the expression data for RORB
expression_data <- FetchData(obj, vars = "RORB")

# Binarize expression data: 1 if RORB is expressed, 0 otherwise
expressing_cells <- expression_data$RORB > 0

# Manually assign names to expressing_cells (i.e., cell names from expression_data)
names(expressing_cells) <- rownames(expression_data)

# Extract cell names that express RORB
expressing_cell_ids <- names(expressing_cells)[expressing_cells]

# Check the result
print(expressing_cell_ids)

# Subset the nearest neighbor graph to include only expressing cells and their neighbors
nn_subgraph <- nn_graph[expressing_cell_ids, expressing_cell_ids]

# Convert the neighbor matrix to an igraph object
graph <- graph_from_adjacency_matrix(nn_subgraph, mode = "undirected", diag = FALSE)

# Plot the network
plot(graph, vertex.size = 5, vertex.label.cex = 0.7, vertex.label.color = "black", 
     main = "Network of Cells Expressing RORB")
