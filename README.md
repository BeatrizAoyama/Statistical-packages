# Statistical packages on R

Link:

Sleuth (https://github.com/pachterlab/sleuth)

OBS: This statistical package only runs out with **Kallisto** aligment software

# **Script:**

library("sleuth")

data_count <- read.table("count_data.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

sleuth_prep <- sleuth_prep(data_count, extra_bootstrap_summary = TRUE)

sleuth_fit <- sleuth_fit(sleuth_prep, ~condition, 'full')

sleuth_fit <- sleuth_fit(sleuth_fit, ~1, 'reduced')

sleuth_lrt <- sleuth_lrt(sleuth_lrt, 'reduced', 'full')

sleuth_lrt <- sleuth_lrt(sleuth_fit, 'reduced', 'full')

sleuth_table2 <- sleuth_results(sleuth_fit, 'reduced:full', 'lrt', show_all = FALSE)

sleuth_significant2 <- dplyr::filter(sleuth_table2, qval <= 0.05)

plot_pca(sleuth_fit, color_by = 'condition')

# Script explanation:

**library("sleuth")** --> This loads the Sleuth package

**data_count <- read.table("count_data.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)** --> The script reads a CSV file containing count data, with the first row as **column headers (header = TRUE)**

**sleuth_prep <- sleuth_prep(data_count, extra_bootstrap_summary = TRUE)** --> This initializes the Sleuth object and prepares the data for downstream analysis. The **extra_bootstrap_summary = TRUE** argument generates additional summaries from the bootstrap data, which improves accuracy and variance estimation

**sleuth_fit <- sleuth_fit(sleuth_prep, ~condition, 'full')** --> A full statistical model is fitted to the data, with **~condition** specifying that gene or transcript expression is modeled based on the **condition variable** (e.g., experimental groups)

**sleuth_fit <- sleuth_fit(sleuth_fit, ~1, 'reduced')** --> A reduced statistical model is created, where gene or transcript expression is not associated with condition (using only an intercept, ~1). This model is used for comparison with the full model

**sleuth_lrt <- sleuth_lrt(sleuth_fit, 'reduced', 'full')** --> A likelihood ratio test is performed to compare the reduced model ('reduced') with the full model ('full'). This test identifies genes or transcripts whose expression is significantly affected by the condition variable

**sleuth_table2 <- sleuth_results(sleuth_fit, 'reduced:full', 'lrt', show_all = FALSE)** --> The results of the likelihood ratio test are extracted. The show_all = FALSE argument ensures that only the most relevant or significant results are included

**sleuth_significant2 <- dplyr::filter(sleuth_table2, qval <= 0.05)** --> The significant results are filtered using the qval column, which contains adjusted p-values (to account for multiple testing). Only results with qval ≤ 0.05 are retained

**plot_pca(sleuth_fit, color_by = 'condition')** --> A Principal Component Analysis (PCA) plot is generated to visualize the overall variance and clustering of samples based on the condition variable
