library(dplyr)
library(tidyr)
library(tibble)
library(limma)
library(GEOquery)
library(writexl)
library(PCAtools)
library(knitr)
library(ggplot2)

#data retrieval####
#CLL data (getting files) 
getGEOSuppFiles("GSE75122")
untar("GSE75122/GSE75122_RAW.tar", exdir = "dataCLL/")

#RS data
getGEOSuppFiles("GSE171481")
untar("GSE171481/GSE171481_RAW.tar", exdir = "dataRS/")

#reading the files
targetsfile <- read.delim("targets.txt")
require(limma)

targetsinfo <- targetsfile

project <- read.maimages(
  targetsinfo,
  source = 'agilent.median',
  green.only = TRUE,
  other.columns = 'gIsWellAboveBG') #retain info about the background
colnames(project) <- sub('raw\\/', '', colnames(project))  

#background normalisation
project.bgcorrect <- backgroundCorrect(project, method = 'normexp')

#data normalisation (quantile method)
project.bgcorrect.norm <- normalizeBetweenArrays(project.bgcorrect,
                                                 method = 'quantile')
 # QC filtering ####
# filter out:
#   - control probes
#   - probes with no name
#   - probes whose signal falls below background in 3 or more samples
Control <- project.bgcorrect.norm$genes$ControlType==1L
NoSymbol <- is.na(project.bgcorrect.norm$genes$GeneName)
IsExpr <- rowSums(project.bgcorrect.norm$other$gIsWellAboveBG > 0) >= 3
project.bgcorrect.norm.filt <- project.bgcorrect.norm[!Control & !NoSymbol & IsExpr, ]
dim(project.bgcorrect.norm)
dim(project.bgcorrect.norm.filt) 

# for replicate probes, replace values with the mean
# ID is used to identify the replicates
project.bgcorrect.norm.filt.mean <- avereps(project.bgcorrect.norm.filt,
                                            ID = project.bgcorrect.norm.filt$genes$ProbeName)
dim(project.bgcorrect.norm.filt.mean)

#OUTPUT ####
require(PCAtools)
p <- pca(
  project.bgcorrect.norm.filt.mean$E,
  metadata = project.bgcorrect.norm.filt.mean$targets)
biplot(p, colby = 'group', legendPosition = 'bottom')    

# DGE analysis ####
targets <- project.bgcorrect.norm.filt.mean$targets
targets$group <- factor(targets$group, levels = c('RS','CLL'))

# Create the study design
design <- model.matrix(~ 0 + targets$group)
colnames(design) <- c('RS', 'CLL')

# Fit the linear model on the study's data
project.fitmodel <- lmFit(project.bgcorrect.norm.filt.mean, design)

# Applying the empirical Bayes method to the fitted values
# Acts as an extra normalisation step and aims to bring the different probe-wise variances to common values
project.fitmodel.eBayes <- eBayes(project.fitmodel)
names(project.fitmodel.eBayes)

# Make individual contrasts: RS vs CLL
res <- makeContrasts(res = 'RS-CLL', levels = design)
res.fitmodel <- contrasts.fit(project.fitmodel.eBayes, res)
res.fitmodel.eBayes <- eBayes(res.fitmodel)

toptable <- topTable(
  res.fitmodel.eBayes,
  adjust = 'BH',
  coef = 'res',
  number = Inf,
  p.value = 1)

write.csv(toptable, "DEG U-RT1 vs CLL.txt")

# Comparison of apoptotic genes ####
apo_genes <- c("MCL1",
               "BCL2",
               "BCL2L1",
               "CDK9")

# Selecting apo-related genes
apo_mt <- project.bgcorrect.norm.filt.mean[project.bgcorrect.norm.filt.mean$genes$GeneName %in% apo_genes, ]
dim(apo_mt)
head(apo_mt)

# Function to select the transcript isoform having the highest average expression
select_highest_isoform <- function(data) {
  #calculate the average of the expression columns for each row, ignoring NAs
  data$average <- rowMeans(data[,-1], na.rm = TRUE) #assuming the first column is the symbol
  
  #group by symbol and select the row with the highest average for each symbol
  selected <- data %>%
    group_by(symbol) %>%
    filter(average == max(average, na.rm = TRUE)) %>% #ensure to handle NAs in max
    slice(1) %>%  #keep only the first occurrence in case of ties
    ungroup() %>%
    select(-average)  #optionally remove the 'average' column
  
  return(selected)
}

# Selecting highest isoform
expr_data <- data.frame(symbol = apo_mt$genes$GeneName,
                        apo_mt$E)
expr_data <- select_highest_isoform(expr_data)
expr_data <- as.data.frame(expr_data)
rownames(expr_data) <- expr_data$symbol
colnames(expr_data)[-1] <- make.names(apo_mt$targets$group, unique = T)
expr_data <- as.data.frame(expr_data)

#write_xlsx(expr_data, "exp_data_apo.xlsx")

# Performing PCA
meta_apo <- data.frame(samples = colnames(expr_data)[-1],
                       condition = c(rep("U-RT1", 3), rep("CLL", 10)),
                       row.names = colnames(expr_data)[-1])
p_apo <- pca(expr_data[,-1], meta_apo, scale = T)
biplot(p_apo,
       showLoadings = T,
       colby = "condition",
       encircle = T,
       ellipse = F,
       lab = NULL,
       legendPosition = "right")

# Barplot visualisation
# reshaping expr_data to long format
expr_long <- expr_data %>%
  pivot_longer(cols = -symbol, names_to = "samples", values_to = "expression")

# merging with metadata
combined_data <- expr_long %>%
  left_join(meta_apo, by = "samples")

# calculating mean and standard error for each group
summary_data <- combined_data %>%
  group_by(symbol, condition) %>%
  summarise(
    mean_expression = mean(expression, na.rm = TRUE),
    se_expression = sd(expression, na.rm = TRUE) / sqrt(n()),  # Standard Error
    .groups = 'drop'
  )
print(summary_data) #checking SD

# performing t-tests (two groups)
t_test_results <- combined_data %>%
  group_by(symbol) %>%
  summarise(t_test_p_value = t.test(expression ~ condition, var.equal = T)$p.value, .groups = 'drop') #because of similar SD between groups, var.equal = T

# adding significance levels based on p-values
t_test_results <- t_test_results %>%
  mutate(significance = case_when(
    t_test_p_value < 0.001 ~ "***",
    t_test_p_value < 0.01 ~ "**",
    t_test_p_value < 0.05 ~ "*",
    TRUE ~ ""
  ))

print(t_test_results)

# merging t-test results with summary data
summary_data <- summary_data %>%
  left_join(t_test_results, by = "symbol")

# creating bar plot with error bars and jitter
ggplot(summary_data, aes(x = symbol, y = mean_expression, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_expression - se_expression, 
                    ymax = mean_expression + se_expression),
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_jitter(data = combined_data, aes(y = expression, color = condition), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 1, alpha = 0.6) +
  geom_text(aes(label = significance), 
            position = position_dodge(width = 0.8), vjust = -0.5) +  # Adjust vjust for positioning
  labs(title = "Microarray expression data", x = "Gene Symbol", y = "Expression Level") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# t-test results
print(t_test_results)

# making good-looking table
require(kableExtra)
rownames(expr_data) <- NULL
colnames(expr_data)[-1] <- make.names(c(rep("U-RT1", 3), rep("CLL", 10)), unique = T)
kable(expr_data, format = "html", caption = "Microarray expression data") %>%
  kable_styling("striped", full_width = F)
