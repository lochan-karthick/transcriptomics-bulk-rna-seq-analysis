# **RNA-seq Differential Expression Analysis: Bulk RNA-seq Analysis Project**

This repository contains the coursework project for RNA-seq analysis focusing on **differential expression** analysis using **DESeq2**. The project involves processing raw RNA sequencing data, performing quality control, identifying differentially expressed genes, and visualizing results using various statistical methods.

---

## **Project Goals**

### **1. Data Preprocessing and Quality Control**
- Convert raw sequencing data (in **fastq** format) into gene-based and transcript-based count tables.
- Tools used for preprocessing include:
  - **Scythe** and **Sickle** for adapter trimming.
  - **Hisat2** for alignment.
  - **Stringtie** for transcript assembly.
  - **prepDE.py** for generating count tables.
- Perform **quality control** to ensure the data is ready for analysis.

### **2. Differential Expression Analysis**
- Conduct **differential expression analysis** using the **DESeq2** R package.
- Compare gene expression across experimental groups **A**, **B**, and **C**, and identify significantly different genes between conditions (**B vs. A**, **C vs. A**).

### **3. Visualization of Results**
- Generate visualizations to better understand the data:
  - **Dispersion plots** to assess data variability.
  - **PCA plots** for visualizing sample relationships based on gene expression.
  - **MA (MvA) plots** for comparing expression levels across conditions.

### **4. Statistical Testing**
- Use **DESeq2** to test for differentially expressed genes and calculate **log2 fold changes** (LFC) for multiple contrasts.
- **Apply shrinkage** to LFC values for more accurate and stable estimates.
- **Batch correction** to eliminate unwanted technical variation using batch correction methods.
- Visualize corrected data with **PCA**.

### **5. Documentation and Reporting**
- Document the **data processing pipeline**, including the tools used and settings.
- Include a discussion of the **statistical results** and interpretations.
- Compare **log2 normalized counts** and **rlog-transformed data**, explaining their effects on the analysis and variance.

---

## **Data Used**

The data used for this analysis comes from a **Bulk RNA-seq experiment**. The dataset consists of RNA sequencing data from samples categorized into three experimental conditions:

- **Condition A**: Control or baseline condition.
- **Condition B**: Experimental condition 1.
- **Condition C**: Experimental condition 2.

### **Data Description:**
- **Sample Data**: The dataset includes RNA-seq data for multiple samples under each experimental condition (A, B, and C).
- **Raw Data**: The data is provided in **fastq** format, which was processed and aligned to a reference genome to generate count matrices.
- **Processed Data**: After quality control and alignment, the raw counts were used for differential expression analysis.

The **gene count matrix** was generated to measure the expression levels of genes across all samples. The following files are key outputs from this analysis:
- **`gene_count_matrix.csv`**: Contains the gene expression data in the form of raw counts per gene and sample.
- **`Significant.Genes.BvsA.LFC1.csv`**: Contains the differentially expressed genes for the **B vs. A** comparison.
- **`Significant.Genes.CvsA.LFC1.csv`**: Contains the differentially expressed genes for the **C vs. A** comparison.

---

## **Key Findings from the Report**

### **1. Data Quality and Preprocessing**
- **Contamination Rate**: The preprocessing pipeline effectively removed contamination, with contamination rates in the range of **10^-6** for all samples.
- **Alignment**: Over **85% of the reads** aligned uniquely for all samples, indicating good quality control. All samples achieved **100% alignment** for the reads.
- **Trimming**: Adapter trimming was performed using **Scythe**, and base quality trimming was done with **Sickle**. The quality of reads after trimming was optimal for downstream analysis.

### **2. Differential Expression Analysis with DESeq2**
- **Dispersion Estimation**: Differential expression analysis highlighted that genes with **low mean counts** had higher dispersion. **Shrinkage** was applied to reduce this variance, making the results more stable.
- **Log Fold Change (LFC) Shrinkage**: The use of **LFC shrinkage** prevented noise from lowly expressed genes, resulting in more reliable estimates of gene expression differences.
- **PCA and MA Plots**: The **PCA plot** showed clear separation between groups, particularly with group B standing out, suggesting strong biological differences. **MA plots** illustrated the number of **differentially expressed genes** with both **LFC=0** and **LFC<1**.

### **3. Transformation and Batch Correction**
- **Rlog Transformation**: The **rlog transformation** improved the homoscedasticity of the data compared to **log2 transformation**, making the analysis more robust.
- **Batch Correction**: Batch correction improved sample separation in **PCA plots**, particularly between groups A and C, where previously the overlap was high.

### **4. Statistical Testing and Null Hypothesis**
- **LFC = 0 vs. LFC < 1**: The analysis compared two null hypotheses:
  - **LFC = 0**: Genes with non-zero fold changes were considered significant.
  - **LFC < 1**: Only genes with **LFC > 1** were considered significant to reduce noise.
- **Shrinkage vs. Unshrunken Data**: **Shrinkage** narrowed the range of LFC estimates, reducing false positives and focusing on truly significant results, particularly useful with small sample sizes.

---

## **Files**

- **`transcriptomics_ass1_partb.Rmd`**: R Markdown file containing the complete **analysis workflow** for RNA-seq, including preprocessing, differential expression, and visualizations.
- **`Significant.Genes.BvsA.LFC1.csv`**: Differentially expressed genes for **B vs. A** comparison.
- **`Significant.Genes.CvsA.LFC1.csv`**: Differentially expressed genes for **C vs. A** comparison.
- **`gene_count_matrix.csv`**: Raw gene count data.
- **`transcriptomics_ass1_partb.html`**: HTML output generated from the R Markdown file.
- **`transcriptomics_ass1a.sh`**: Shell script for preprocessing steps and aligning raw data.

---

## **How to Run the Analysis**

1. **Clone the repository** to your local machine:
   ```bash
   git clone https://github.com/lochan-karthick/your-repository-name.git
