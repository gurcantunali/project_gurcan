#!/usr/bin/env python
# coding: utf-8

# In[1]:


from pydeseq2.dds import DeseqDataSet
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy import stats
from pydeseq2.ds import DeseqStats
from statsmodels.stats.multitest import multipletests
import openpyxl
from adjustText import adjust_text



# In[2]:


counts = pd.read_csv('Proteomics_Data_Dummy.csv')
counts


# In[3]:


# Remove rows that contain any NaN values
counts = counts.dropna()
counts


# In[4]:


counts = counts.set_index(['Accession', 'Description'])
counts


# In[5]:


# Plot a boxplot for all proteins
counts.boxplot(figsize=(10, 6))  # You can adjust the figsize for larger datasets

# Add labels and title
plt.title('Proteomics Data')
plt.xlabel('Conditions')
plt.ylabel('Expression')

# Show the plot
plt.show()


# In[6]:


unlisted_data = counts.values.flatten()

# Create a histogram
plt.figure(figsize=(10, 6))
plt.hist(unlisted_data, bins=30, edgecolor='black', alpha=0.7)  # Adjust bins as needed
plt.title('Histogram of Unlisted Values from Counts')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.grid(axis='y')

# Show the plot
plt.show()


# In[7]:


log_transformed_data = np.log2(counts + 1)

# Flatten the log-transformed data
unlisted_log_data = log_transformed_data.values.flatten()  # Use .values to get a NumPy array and flatten it

# Create a histogram of the log-transformed data
plt.figure(figsize=(10, 6))
plt.hist(unlisted_log_data, bins=30, edgecolor='black', alpha=0.7)  # Adjust bins as needed
plt.title('Histogram of Log-Transformed Values from Counts')
plt.xlabel('Log2 Transformed Value')
plt.ylabel('Frequency')
plt.grid(axis='y')

# Show the plot
plt.show()


# In[8]:


counts_log = np.log2(counts + 1)
counts_log = counts_log.T
counts_log.T


# In[9]:


counts_log


# In[10]:


scaler = StandardScaler()
counts_scaled = scaler.fit_transform(counts_log)

# Perform PCA
pca = PCA()
principal_components = pca.fit_transform(counts_scaled)

# Create a DataFrame for the principal components
pca_df = pd.DataFrame(data=principal_components, 
                      columns=[f'PC{i+1}' for i in range(principal_components.shape[1])])

# Explained variance
explained_variance = pca.explained_variance_ratio_

# Display the explained variance for each principal component
explained_variance


# In[11]:


plt.figure(figsize=(8,6))
plt.bar(range(1, len(explained_variance) + 1), explained_variance, alpha=0.7, color='b')
plt.title('Bar Plot: Explained Variance by Principal Components')
plt.xlabel('Principal Component')
plt.ylabel('Explained Variance Ratio')
plt.xticks(range(1, len(explained_variance) + 1))  # Set x-axis ticks to match component numbers
plt.grid(True)
plt.show()


# In[12]:


col1 = "blue"  # Group 1 color
col2 = "red"   # Group 2 color

# Prepare the scatter plot for the first 6 rows (samples)
plt.figure(figsize=(8,6))

# Plot the first 3 samples (rows 1-3) from Group 1
plt.scatter(principal_components[0:3, 0], principal_components[0:3, 1], color=col1, label='Group 1')

# Plot the next 3 samples (rows 4-6) from Group 2
plt.scatter(principal_components[3:6, 0], principal_components[3:6, 1], color=col2, label='Group 2')

# Add labels for Group 1 (first 3 rows)
for i in range(3):
    plt.text(principal_components[i, 0], principal_components[i, 1], counts_log.index[i], 
             color=col1, fontsize=12, ha='right')

# Add labels for Group 2 (next 3 rows)
for i in range(3, 6):
    plt.text(principal_components[i, 0], principal_components[i, 1], counts_log.index[i], 
             color=col2, fontsize=12, ha='left')

# Labels for axes
plt.xlabel("PC1")
plt.ylabel("PC2")

# Add legend and show plot
plt.legend()
plt.grid(True)
plt.title("PCA Plot of First 6 Samples (Rows)")
plt.show()


# In[13]:


# Calculate the mean across columns 1 to 3 (M) and columns 4 to 6 (Hypoxia)
mean_M = counts.iloc[:, 0:3].mean(axis=1)
mean_C = counts.iloc[:, 3:6].mean(axis=1)

# Create a new DataFrame with Accession, counts, mean_C, and mean_M
counts_out = pd.concat([
    counts,               # Adding the entire exprdata DataFrame
    mean_M.rename('mean_M'),  # Adding mean_M column
    mean_C.rename('mean_C')   # Adding mean_C column
], axis=1)

# Display the first few rows (equivalent to head() in R)
counts_out.head()


# In[14]:


# Add a new column for log2 Fold Change (log2FC)
counts_out['log2FC'] = np.log2(counts_out['mean_M'] / counts_out['mean_C'])

# Display the first few rows
counts_out.head()


# In[15]:


# Define the Fold Change threshold
FC = 2  # Change this value as needed

# Filter the DataFrame for log2FC > FC
upregulated = counts_out[counts_out['log2FC'] > FC]

# Filter the DataFrame for log2FC < -FC
downregulated = counts_out[counts_out['log2FC'] < -FC]

# Display the results
print("Upregulated genes (log2FC > FC):")
print(upregulated)

print("\nDownregulated genes (log2FC < -FC):")
print(downregulated)


# In[16]:


# Define a function to perform t-test and return p-value
def get_pvalue(row):
    # Perform t-test on the first three and last three elements
    return stats.ttest_ind(row[:3], row[3:]).pvalue

# Apply the function to each row in exprdata_log
p_values = counts_log.T.apply(get_pvalue, axis=1)

# Add p-values to exprdata_out
counts_out['pvalue'] = p_values

# Display the updated exprdata_out DataFrame
counts_out



# In[17]:


# Apply FDR adjustment using Benjamini-Hochberg method
counts_out['fdr'] = multipletests(counts_out['pvalue'], method='fdr_bh')[1]

# Display the updated DataFrame
counts_out


# In[18]:


# Define the threshold
alfa = 0.05  # Change this value as needed

# Count how many p-values and fdr values are below the threshold
count_pvalue_below_alfa = (counts_out['pvalue'] < alfa).sum()
count_fdr_below_alfa = (counts_out['fdr'] < alfa).sum()

# Display the results
print(f"Number of p-values below {alfa}: {count_pvalue_below_alfa}")
print(f"Number of FDR values below {alfa}: {count_fdr_below_alfa}")


# In[19]:


filtered_set = counts_log.T[counts_out['fdr'] < 0.05]
filtered_set


# In[20]:


allcounts_out = counts_out.reset_index()
allcounts_out


# In[21]:


# Set row names based on the 'Accession' column in exprdata_out
filtered_set.index = allcounts_out[allcounts_out['fdr'] < 0.05]['Accession']
filtered_set


# In[22]:


set_out = filtered_set.reset_index()
set_out


# In[23]:


# Convert filtered_set to a NumPy matrix (heatdata2)
heatdata2 = set_out.values

# Show the matrix
heatdata2


# In[24]:


# Filter out non-numeric columns from filtered_set
heatdata2_numeric = set_out.select_dtypes(include=[float, int]).to_numpy()

# Create the heatmap with numeric data
plt.figure(figsize=(7, 5))  # Adjust the figure size
sns.heatmap(heatdata2_numeric, cmap='coolwarm', cbar=True)

# Display the plot
plt.show()


# In[25]:


# Filter the data for fdr < 0.05 and log2FC conditions
filtered_set_upregulated = counts_log.T[(counts_out['fdr'] < 0.05) & (counts_out['log2FC'] > 0)]
filtered_set_downregulated = counts_log.T[(counts_out['fdr'] < 0.05) & (counts_out['log2FC'] < 0)]

# Combine the upregulated and downregulated sets
filtered_setdiff = pd.concat([filtered_set_upregulated, filtered_set_downregulated], axis=0)
filtered_setdiff


# In[26]:


# Filter the data for fdr < 0.05
filtered_counts_out = counts_out[counts_out['fdr'] < 0.05]

# Select most upregulated (top 25) and downregulated (top 25) proteins
top_25_upregulated = filtered_counts_out.nlargest(25, 'log2FC')
top_25_downregulated = filtered_counts_out.nsmallest(25, 'log2FC')

# Keep the counts_log.T for the selected proteins
upregulated_proteins25 = counts_log.T.loc[top_25_upregulated.index]
downregulated_proteins25 = counts_log.T.loc[top_25_downregulated.index]

# Combine the top 25 upregulated and downregulated proteins
combined_updownproteins25 = pd.concat([upregulated_proteins25, downregulated_proteins25])


# In[27]:


combined_updownproteins25


# In[28]:


# Filter out non-numeric columns from filtered_set
heatdata3_numeric = combined_updownproteins25.select_dtypes(include=[float, int]).to_numpy()

# Create the heatmap with numeric data
plt.figure(figsize=(7, 5))  # Adjust the figure size
sns.heatmap(heatdata3_numeric, cmap='coolwarm', cbar=True)

# Display the plot
plt.show()


# In[29]:


# Extract Accession labels if using MultiIndex
accession_labels = combined_updownproteins25.index.get_level_values(0)  # Change '0' if accessions are on a different level

# Create the heatmap with Accession labels on the y-axis and column names on the x-axis
plt.figure(figsize=(8, 15))  # Adjust the figure size
heatmap = sns.heatmap(
    heatdata3_numeric, 
    cmap='coolwarm', 
    cbar=True, 
    yticklabels=accession_labels,
    xticklabels=combined_updownproteins25.columns,  # Use column names for x-axis labels
    cbar_kws={"pad": 0.2}  # Adjust the pad value to increase space between heatmap and colorbar
)

# Move y-tick labels to the right
heatmap.yaxis.tick_right()
heatmap.yaxis.set_ticks_position('right')

# Set y-tick label rotation to horizontal
plt.yticks(rotation=0)  # Set rotation to 0 degrees for horizontal labels

# Set plot title
plt.title('Top 25 Upregulated and Downregulated Proteins', fontsize=14)

# Display the heatmap
plt.show()


# In[30]:


combined25 = combined_updownproteins25.reset_index()
combined25


# In[31]:


# Extract protein names from the combined50 DataFrame for y-axis labels
Description = combined25['Description'].values  # Get the values as a numpy array

# Create the heatmap with protein names on the y-axis and column names on the x-axis
plt.figure(figsize=(8, 18))  # Adjust the figure size
heatmap = sns.heatmap(
    heatdata3_numeric, 
    cmap='coolwarm', 
    cbar=True, 
    yticklabels=Description,  # Use protein numbers for y-axis labels
    xticklabels=combined_updownproteins25.columns,  # Use column names for x-axis labels
    cbar_kws={"pad": 0.2}  # Adjust the pad value to increase space between heatmap and colorbar
)

# Move y-tick labels to the right
heatmap.yaxis.tick_right()
heatmap.yaxis.set_ticks_position('right')

# Set y-tick label rotation to horizontal
plt.yticks(rotation=0)  # Set rotation to 0 degrees for horizontal labels

# Set plot title
plt.title('Top 25 Upregulated and Downregulated Proteins', fontsize=14)

# Display the heatmap
plt.show()


# In[32]:




# Define thresholds
pCutoff = 0.05
FCcutoff = 2.0  # Log2 fold-change cutoff

# Create significance flags
allcounts_out['significant'] = (allcounts_out['fdr'] < pCutoff) & (abs(allcounts_out['log2FC']) > FCcutoff)
allcounts_out['pval_significant'] = (allcounts_out['fdr'] < pCutoff) & (abs(allcounts_out['log2FC']) <= FCcutoff)

# Merge dataframes on the 'Accession' column
merged_data = allcounts_out.merge(combined25, on='Accession', how='left')

# Assign colors based on significance
colors = np.where(merged_data['significant'], 'red',
         np.where(merged_data['pval_significant'], 'blue', 'green'))

# Create the volcano plot
plt.figure(figsize=(12, 10))
plt.scatter(merged_data['log2FC'], -np.log10(merged_data['fdr']), c=colors, alpha=0.75)

# Add significance threshold lines
plt.axhline(y=-np.log10(pCutoff), color='grey', linestyle='--', linewidth=1)
plt.axvline(x=FCcutoff, color='grey', linestyle='--', linewidth=1)
plt.axvline(x=-FCcutoff, color='grey', linestyle='--', linewidth=1)

# Collect text labels to adjust later
texts = []

for i, row in merged_data.iterrows():
    protein_name = row['Description_x'] if pd.notna(row['Description_x']) else row['Description_x']
    
    # Check if the protein is in the list and if it's significant
    if protein_name in combined25['Description'].values or protein_name in combined25['Description'].values:
        texts.append(plt.text(row['log2FC'], -np.log10(row['fdr']), protein_name, fontsize=8))

# Adjust text to reduce overlap
adjust_text(
    texts, 
    arrowprops=dict(arrowstyle='->', color='black', lw=0.5),  # Optional: Add arrows for clarity
    expand_points=(1.2, 1.2),  # Adjusts spacing around points
    expand_text=(1.2, 1.2),  # Adjusts spacing around text
    force_text=0.7,  # Reduces overlap force
    force_points=0.3  # Adjusts distance of text from points
)


# Set axis labels
plt.xlabel('log2 Fold Change', fontsize=11)
plt.ylabel('-log10(FDR)', fontsize=11)

# Set x and y axis limits
plt.xlim([-8, 8])  # Adjust limits as needed
plt.ylim([0, 3.5])  # Adjust limits as needed

# Add plot title
plt.title('Volcano Plot', fontsize=14)

# Add custom legend
green_patch = plt.scatter([], [], color='green', label='NS')  # Non-significant
blue_patch = plt.scatter([], [], color='blue', label='p-value < 0.05')  # p-value < 0.05, but |log2FC| <= 2
red_patch = plt.scatter([], [], color='red', label='p-value < 0.05 and FC > 2')  # p-value < 0.05 and |log2FC| > 2

plt.legend(handles=[green_patch, blue_patch, red_patch], loc='upper left')

# Show plot
plt.show()


# In[33]:


# Define thresholds
pCutoff = 0.05
FCcutoff = 2.0  # Log2 fold-change cutoff

# Create significance flags
allcounts_out['significant'] = (allcounts_out['fdr'] < pCutoff) & (abs(allcounts_out['log2FC']) > FCcutoff)
allcounts_out['pval_significant'] = (allcounts_out['fdr'] < pCutoff) & (abs(allcounts_out['log2FC']) <= FCcutoff)

# Merge dataframes on the 'Accession' column
merged_data = allcounts_out.merge(combined25, on='Accession', how='left')

# Assign colors based on significance
colors = np.where(merged_data['significant'], 'red',
         np.where(merged_data['pval_significant'], 'blue', 'green'))

# Create the volcano plot
plt.figure(figsize=(12, 10))

# Plot all points with no fill (transparent) and outlined
plt.scatter(merged_data['log2FC'], -np.log10(merged_data['fdr']), c=colors, alpha=0.75, facecolor='none', edgecolor='gray', linewidth=0.5)

# Add significance threshold lines
plt.axhline(y=-np.log10(pCutoff), color='grey', linestyle='--', linewidth=1)
plt.axvline(x=FCcutoff, color='grey', linestyle='--', linewidth=1)
plt.axvline(x=-FCcutoff, color='grey', linestyle='--', linewidth=1)

# Collect text labels to adjust later
texts = []

# Plot outlined points for labeled proteins
for i, row in merged_data.iterrows():
    protein_name = row['Description_x'] if pd.notna(row['Description_x']) else row['Description_x']
    
    # Check if the protein is in the list and if it's significant
    if protein_name in combined25['Description'].values or protein_name in combined25['Description'].values:
        # Plot the labeled protein with filled color and outline
        color = 'red' if row['significant'] else 'blue'  # Choose color based on significance
        plt.scatter(row['log2FC'], -np.log10(row['fdr']),
                    c=color, alpha=0.75, edgecolor='black', linewidth=1)  # Filled color for label
        
        # Add the label
        texts.append(plt.text(row['log2FC'], -np.log10(row['fdr']), protein_name, fontsize=8, ha='right'))

# Adjust text to reduce overlap
adjust_text(
    texts, 
    arrowprops=dict(arrowstyle='->', color='black', lw=0.5),  # Optional: Add arrows for clarity
    expand_points=(1.2, 1.2),  # Adjusts spacing around points
    expand_text=(1.2, 1.2),  # Adjusts spacing around text
    force_text=0.7,  # Reduces overlap force
    force_points=0.3  # Adjusts distance of text from points
)

# Set axis labels
plt.xlabel('log2 Fold Change', fontsize=11)
plt.ylabel('-log10(FDR)', fontsize=11)

# Set x and y axis limits
plt.xlim([-8, 8])  # Adjust limits as needed
plt.ylim([0, 3.5])  # Adjust limits as needed

# Add plot title
plt.title('Volcano Plot', fontsize=14)

# Add custom legend
green_patch = plt.scatter([], [], color='green', label='NS')  # Non-significant
blue_patch = plt.scatter([], [], color='blue', label='p-value < 0.05')  # p-value < 0.05, but |log2FC| <= 2
red_patch = plt.scatter([], [], color='red', label='p-value < 0.05 and FC > 2')  # p-value < 0.05 and |log2FC| > 2

plt.legend(handles=[green_patch, blue_patch, red_patch], loc='upper left')

# Show plot
plt.show()





