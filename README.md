# project_gurcan
## Project Overview
In this project, we observed the development of AML in the mouse model to gain more information about the transformation of AML. In order to understand the mechanisms underlying the development of the disease and to determine target molecules in AML treatment, proteomics analysis was performed on biological samples obtained from the AML disease (M) and disease-free control (C) groups. The obtained data were deeply analyzed using Python and Jupyter.
At first, packages needed were installed, and data was read and processed. Distribution of the data was checked. According to log transformed values histogram t was observed that the data was normally distributed. In the PCA analysis, The first PC is the highest one but the first 3 PC’s, which have 50-80%, can explain the variability of most of our data. Control (C) and AML (M) conditions showed a high level of separation, which means difference between groups is high. AML 3 has the highest PC and Control 1 and 2 show the lowest PC. Mean values ​​of the groups were calculated. Then, Log2 values ​​were computed to see the fold change. To determine statistical significance, pvalue was calculated using ttest. To get FDR values, FDR adjustment was applied. According to FDR values, 686 proteins with p<0.05 are significantly disregulated. 630 upregulated and 56 downregulated. 
Heatmap was created for disregulated proteins. Since many proteins were upregulated and downregulated. Trashold value = 2 was selected for Log2FC. A heatmap was created for the top 55 proteins that were upregulated and downregulated. 
Then, a volcano plot was created for the upregulated and downregulated proteins. NS (the green dots) refers the proteins with non-significant changes. fdr < 0.05 (blue dots) indicates significantly up- or down-regulated proteins in AML.fdr < 0.05 and Log2FC < -2 & Log2FC > 2 (the red dots) show significantly up-regulated and down-regulated proteins in AML conditions, compared to control, according to threshold value.Since there were too many upregulated proteins, their labels were overlapping each other. Therefore, labels were added for a certain number of proteins. In addition, an attempt was made to add an outline to the dots to determine which dots the labels belong to. Therefore, some corrections were made. However, it did not reach the desired level completely. After the error specified in the last code was fixed, a run was performed, but the process took a long time and did not produce results. 
In conclusion, upregulated and downregulated proteins were determined in the analyzed data due to the AML effect.

## Installation
To run this project, clone the repository and navigate to the project directory.
git clone https://github.com/gurcantunali/project_gurcan.git
cd project_gurcan

## Usage
1. Place the input file `Proteomics_Data.csv` in the same directory as the script or update the `input_file` path in the script to the location of your CSV file.
2. Run the script:
   python Proteomics-Project.py

## Dependencies
- PyDESeq2 (DeseqDataSet and DeseqStats): These classes are central to performing differential expression analysis with PyDESeq2. DeseqDataSet sets up the data structure, while DeseqStats runs the statistical tests.
- pandas: For data manipulation and processing
- numpy: For numerical operations
- matplotlib and seaborn: For data visualization
- scipy: For statistical analysis
- StandardScaler and PCA: StandardScaler scales data for analysis, while PCA reduces the dimensionality, which is especially useful for visualization and initial exploratory analysis.
- Multipletest: For adjusting p-values for multiple hypothesis testing. 
- Openpyxl: For reading and writing Excel files

## Data Requirements
The primary dataset needed for this analysis is Proteomics_Data.csv. The columns in the file should include:
- Accession: Protein identifiers
- Description: Description of each protein
- Abundance M (1,2,3): Abundance values for condition M
- Abundance C (1,2,3): Abundance values for condition C

## Expected Output
Summary Table: A summary of differentially expressed proteins with key statistics, such as fold change, fdr and p-values for significance.
Distribution of Data: A plot (such as a histogram or box plot) to visualize the overall distribution of protein abundances and assess data normalization.
Principal Component Analysis (PCA): Display a PCA plot to visualize sample separation across conditions.
Fold Change Calculation: The ratio of mean abundance between conditions, used to quantify changes in protein expression.
Statistical Significance: Statistical tests (e.g., t-tests or DESeq2) to assess whether differences in protein abundance between conditions are significant.
Visualizations: Graphs, such as heatmap and/or volcano plot, to highlight significant proteins, expression patterns, and group differences across conditions. 


