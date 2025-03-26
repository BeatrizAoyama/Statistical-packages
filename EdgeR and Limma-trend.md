# Statistical package on R

Links:

Limma-Trend (https://bioconductor.org/packages/release/bioc/html/limma.html)

EdgeR (https://bioconductor.org/packages/release/bioc/html/edgeR.html)

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

# **Script for Limma**:

OBS: The edgeR package is necessary to run out Limma-trend package
#
library(edgeR)

library(limma)

countData <- as.matrix(read.csv("count_data.csv",sep=",",row.names="gene_id"))

group <- factor(c("s", "s", "s", "s", "s", "s", "s", "c", "c", "c", "c"))

dgList <- DGEList(counts=countData, group=group)

keep <- filterByExpr(dgList)

dgList <- dgList[keep,,keep.lib.sizes=FALSE]

dgList <- calcNormFactors(dgList)

logCPM <- cpm(dgList, log=TRUE, prior.count=3)

plotMDS(logCPM)

design <- model.matrix(~0+group)

colnames(design) <- c("s", "c")

fit <- lmFit(logCPM, design)

contr <- makeContrasts(s_vs_c = s - c, levels=design)

fit.contr <- eBayes(contrasts.fit(fit,contr))

topTable(fit.contr, coef=1, adjust="BH")

result <- topTable(fit.contr, coef=1, adjust="BH", number=10000)

write.table(result, "name_of_the_file.txt")

# Script explained: 

**library(edgeR) and library(limma)** --> These are R packages used for differential gene expression (DGE) analysis. edgeR is for working with count data, while limma is used for linear modeling and statistical analysis

**countData <- as.matrix(read.csv("count_data.csv", sep = ",", row.names = "gene_id"))** --> The script reads a CSV file containing gene expression count data (rows represent genes, columns represent samples) and converts it into a matrix. The **row.names = "gene_id"** ensures that gene IDs serve as row names

**group <- factor(c("s", "s", "s", "s", "s", "s", "s", "c", "c", "c", "c"))** --> This defines a group factor, where "s" represents the *disease* samples in one experimental condition, and "c" represents the *control* samples

**dgList <- DGEList(counts = countData, group = group)** --> This creates a *DGEList object*, which is a data structure used by edgeR to store the count data and the group assignments for DGE analysis

**keep <- filterByExpr(dgList)** and **dgList <- dgList[keep, , keep.lib.sizes = FALSE]** --> Genes with low or insufficient expression across samples are removed using the **filterByExpr()** function. This ensures only genes with meaningful data are analyzed.

**dgList <- calcNormFactors(dgList)** --> This calculates normalization factors to correct for differences in sequencing depth between samples, ensuring fair comparison

**logCPM <- cpm(dgList, log = TRUE, prior.count = 3)** --> The counts are converted into *log-transformed counts per million (logCPM)*. This stabilizes the variance and makes the data more suitable for linear modeling and it is important for the **Multidimensional Scaling (MDS) plot**

**plotMDS(logCPM)** --> The MDS plot is generated to visualize the relationship or clustering of samples based on gene expression as PCA plot

**design <- model.matrix(~0 + group) and colnames(design) <- c("s", "c")** --> A design matrix is created for linear modeling. The ~0 + group specifies that there will be no intercept, and the levels of the group factor ("s" and "c") will be used directly. 

**fit <- lmFit(logCPM, design)** --> A linear model is fit to the log-transformed data for each gene using the design matrix

**contr <- makeContrasts(s_vs_c = s - c, levels = design)** --> **Since the comparision is the group *disease* "s" against *control* "c", the differentially expressed genes will be up-regulated or down-regulated for fcdiia or fcdiib (disease condition) compared to the control (control condition). 
But, this can be changed by switching **makeContrasts(c_vs_s = c - s, levels = design)**, in this case the differentially expressed genes will be up-regulated or down-regulated for control (control condition) compared to the fcdiia or fcdiib (disease condition)**

**fit.contr <- eBayes(contrasts.fit(fit, contr))** --> The contrast is applied to the *linear model* using **contrasts.fit()**, and *empirical Bayes moderation* is performed with **eBayes()** to stabilize variance estimates

**topTable(fit.contr, coef = 1, adjust = "BH")** --> The **topTable()** function extracts the top differentially expressed genes based on the **specified contrast (coef = 1)**. The p-values are adjusted using the *Benjamini-Hochberg (BH) method* to control the false discovery rate (FDR)

**result <- topTable(fit.contr, coef = 1, adjust = "BH", number = 10000) and write.table(result, "name_of_the_file.txt")** --> The results **(top DEGs)** are saved into a text file. The parameter number = 10000 ensures that up to 10,000 genes are included in the output




