# project_gurcan
## Project Overview
In this project, proteomics data from AML (disease, M) and control (disease-free, C) samples were analyzed to explore the mechanisms of AML development and identify potential treatment targets. Using Python and Jupyter, necessary Python packages were installed, and proteomics data was imported, processed, and checked for quality. Log transformation revealed a normal distribution, verified using histograms. Principal Component Analysis (PCA) was conducted to reduce dimensionality and visualize differences between groups. The first three PCs explained 50-80% of data variability, showing clear separation between AML and control samples, with AML 3 showing the most variability and Control 1 & 2 showing the least. Group means were calculated to summarize protein expression levels. Log2 fold changes (Log2FC) were computed to quantify protein expression differences. A t-test was applied to identify statistically significant changes, followed by FDR adjustment to control for false discoveries. Proteins with FDR < 0.05 were deemed significant. 686 significantly dysregulated proteins (FDR < 0.05)(630 upregulated and 56 downregulated). Protein expression patterns visualized in heatmaps. The top 25 upregulated and downregulated proteins (Log2FC > 2 or < -2) were highlighted. Upregulated (red), downregulated (red), and non-significant proteins (green) were showed on volcano plot. Top 25 upregulated and downregulated proteins were labeled for emphasis. In conclusion, many proteins showed significant dysregulation in AML compared to controls, with a threshold of Log2FC > 2 or < -2 identifying the most notable changes. These proteins may provide insights into the molecular underpinnings of AML and serve as potential targets for therapeutic intervention.

## Installation
To run this project, clone the repository and navigate to the project directory.
- git clone https://github.com/gurcantunali/project_gurcan.git
- cd project_gurcan

## Usage
- Place the input file `Proteomics_Data_Dummy.csv` in the same directory as the script or update the `input_file` path in the script to the location of your CSV file.
- Run the script:
  python Proteomics_Dummy-Project.py

## Channels
environment.yml
conda environment: proteomics
- conda-forge
- defaults
- anaconda

## Dependencies
- python=3.10: Programming language and version used in proteomics analysis.
- PyDESeq2 (DeseqDataSet and DeseqStats): These classes are central to performing differential expression analysis with PyDESeq2. DeseqDataSet sets up the data structure, while DeseqStats runs the statistical tests.
- pandas: For data manipulation and processing
- numpy: For numerical operations
- matplotlib and seaborn: For data visualization
- scipy: For statistical analysis
- StandardScaler and PCA: StandardScaler scales data for analysis, while PCA reduces the dimensionality, which is especially useful for visualization and initial exploratory analysis.
- Multipletest: For adjusting p-values for multiple hypothesis testing. 
- Openpyxl: For reading and writing Excel files
- adjustText: To optimize overlapping text placement

## Data Requirements
The primary dataset needed for this analysis is Proteomics_Data.csv. The columns in the file should include:
- Accession: Protein identifiers
- Description: Description of each protein
- Abundance M (1,2,3): Abundance values for condition M
- Abundance C (1,2,3): Abundance values for condition C

## Expected Output
- Summary Table: A summary of differentially expressed proteins with key statistics, such as fold change, fdr and p-values for significance.
- Distribution of Data: A plot (such as a histogram or box plot) to visualize the overall distribution of protein abundances and assess data normalization.
- Principal Component Analysis (PCA): Display a PCA plot to visualize sample separation across conditions.
- Fold Change Calculation: The ratio of mean abundance between conditions, used to quantify changes in protein expression.
- Statistical Significance: Statistical tests (e.g., t-tests or DESeq2) to assess whether differences in protein abundance between conditions are significant.
- Visualizations: Graphs, such as heatmap and/or volcano plot, to highlight significant proteins, expression patterns, and group differences across conditions. 


