# Statistical packages on R

Link:
DESeq2 (https://lashlock.github.io/compbio/R_presentation.html)
#
Let us imagine two conditions, S (disease group, separated into Focal Cortical Dysplasia II type "a" and Focal Cortical Dysplasia II type "b") and Condition C (control group), which will be compared to identify the differentially expressed genes
#
**S1 to S19** represent samples from the experimental group labeled *'disease'* group, while **C26 to C29** belong to the *'control'* group. The count data reflects the number of reads mapped to each gene based on the selected reference genome during the alignment process. See the count data example bellow:

![image](https://github.com/user-attachments/assets/ce74fc63-599d-4e10-9311-4374cc142737)

In red --> Focal Cortical Dysplasia II type "a" samples

In blue --> Focal Cortical Dysplasia II type "b" samples

In green --> Control samples
#

The metadata for the samples (coldata) provides details about the labeling of each experimental group. It is typically divided into two columns: 'sample' and 'condition.' The 'sample' corresponds to the identifiers used in the count data, while the 'condition' specifies the experimental group to which each sample belongs. See the coldata example bellow:

![image](https://github.com/user-attachments/assets/2fe40c44-8e29-47cd-8961-c0ee869d2a1e)

#
# **Script:**

library(DESeq2)

file <- read.csv("count_table_name_of_your_file.csv",row.names="gene_id")

coldata <- read.csv("coldata_name_of_your_file.csv",row.names="amostras")

matriz <- as.matrix(file)

dds <- DESeqDataSetFromMatrix(countData = matriz, colData = coldata, design = ~ condition)

data_DGE <- DESeq(dds) 

data_DGE <- DESeq(dds)

result_count <- results(data_DGE,contrast=c("condition","condition1","condition2"))

result_count <- result_count[order(result_count$padj),]

write.csv(result_count,file="name_your_file_containing_DGE.csv")

# Script explanation:

**library(DESeq2)**: This loads the DESeq2 package

**file <- read.csv("count_table_name_of_your_file.csv", row.names = "gene_id")**: This reads a CSV file containing the count data (gene expression data) into a data frame. The **row.names = "gene_id"** ensures that the rows are indexed by gene IDs

**coldata <- read.csv("coldata_name_of_your_file.csv", row.names = "samples")**: This reads a CSV file containing the metadata for the samples (e.g., conditions or experimental design as described before) into a data frame. The **row.names = "samples"** ensures that rows are indexed by the sample names

**matriz <- as.matrix(file)**: This converts the count data into a matrix format required by DESeq2

**dds <- DESeqDataSetFromMatrix(countData = matriz, colData = coldata, design = ~ condition)**: This creates a **DESeqDataSet object** that combines the **count data (countData)**, the **metadata (colData)**, and the **experimental design (design)**. The **~ condition** indicates that the differential expression analysis will compare samples based on the "condition" variable in the metadata (coldata) described at the file

**data_DGE <- DESeq(dds)**: This runs the DESeq2 pipeline, normalizing the data and calculating the statistics for differential expression analysis

**result_count <- results(data_DGE, contrast = c("condition", "fcdiia", "control"))**: This extracts the results for comparing "fcdiia" versus "control" (as specified in the contrast argument). It generates a table with statistics for each gene, including fold change and adjusted p-values. So, the differentially expressed genes will be up-regulated or down-regulated for fcdiia (disease condition) compared to the control (control condition). But, this can be changed by switching **c("condition", "control", "fcdiia"))**, in this case the differentially expressed genes will be up-regulated or down-regulated for control (control condition) compared to the fcdiia (disease condition)

**result_count <- result_count[order(result_count$padj),]**: This sorts the results table in ascending order of the **adjusted p-values (padj)**, highlighting the most statistically significant differentially expressed genes

**write.csv(result_count, file = "name_your_file_containing_DGE.csv")**: This writes the sorted results to a new CSV file

# Preparing the PCA using ggplot2 package:

rld_data_DGE <- rlog(data_DGE)

print(plotPCA(rld_data_DGE, intgroup=c("condition")))

library(ggplot2)

name_of_PCA_file <- plotPCA(rld_dados_DE, intgroup=c("condition"))

name_of_PCA_file <- plot_bia_nomes + geom_label(aes(label = name))

print(name_of_PCA_file)


# Script explanation:

**rld_data_DGE <- rlog(data_DGE)**: This applies a regularized log transformation to the normalized count data stored in **data_DGE**. The transformation stabilizes the variance across samples, making it easier to compare gene expression values between conditions and it is used in PCA analysis

**print(plotPCA(rld_data_DGE, intgroup = c("condition")))**: This generates a PCA plot based on the transformed data **(rld_data_DGE)**, grouping samples by the condition variable specified in the metadata

**library(ggplot2)**: The ggplot2 library is loaded to enhance the visualization of the PCA plot

**name_of_PCA_file <- plotPCA(rld_dados_DE, intgroup = c("condition"))**: A new PCA plot is created using a dataset named **rld_dados_DE**, which is assumed to contain regularized log-transformed data. The plot groups samples based on the condition variable

**name_of_PCA_file <- plot_bia_nomes + geom_label(aes(label = name))**: This line adds labels to the plot using the **geom_label()** function from ggplot2. The **aes(label = name)** specifies that sample names are added as labels to the PCA plot for easier identification

**print(name_of_PCA_file)** <- The final PCA plot, now customized with labels, is printed for visualization











