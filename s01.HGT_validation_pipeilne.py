import pandas as pd
table_file = "/Users/mac/Desktop/Fu Group/1--Project1/2--DAG3_HGT/1--Recent_HGT/drep/FU/FU_Cdb.csv"

result_df = pd.read_csv(table_file, sep=",")

result_df['genome'] = result_df['genome'].str.replace('.fa$', '', regex=True)

genome_to_cluster = dict(zip(result_df['genome'], result_df['primary_cluster']))

HGT_file = "/Users/mac/Desktop/Fu Group/1--Project1/2--DAG3_HGT/1--Recent_HGT/HGT_events/HGT_events.csv"

df_abun = pd.read_csv(HGT_file, sep = ",")

#print(Cluster_group)
df_abun['MAG1_Group'] = df_abun['MAG 1'].map(genome_to_cluster)
df_abun['MAG2_Group'] = df_abun['MAG 2'].map(genome_to_cluster)
#df_abun.to_csv("/Users/mac/Desktop/Fu Group/1--Project1/2--DAG3_HGT/1--Recent_HGT/HGT_events/HGT_events_group.csv")

df_mag_abundance = pd.read_csv("/Users/mac/Desktop/Fu Group/1--Project1/2--DAG3_HGT/1--Recent_HGT/Species_abun/species_median.csv",index_col=0)
df_mag_existence = df_mag_abundance >= 0.2

sample_columns = df_mag_abundance.columns

def calculate_mag_existence_for_sample(row, sample):
    mag1 = row['MAG1_Group']
    mag2 = row['MAG2_Group']
    
    mag1_exists = df_mag_existence.loc[mag1, sample]
    mag2_exists = df_mag_existence.loc[mag2, sample]

    if mag1_exists and mag2_exists:
        return '2_2'
    elif mag1_exists:
        return '2_0'
    elif mag2_exists:
        return '0_2'
    else:
        return '0_0'

result_data = {'HGT_ID': df_abun['HGT_ID']}

for sample in sample_columns:
    result_data[sample] = df_abun.apply(lambda row: calculate_mag_existence_for_sample(row, sample), axis=1)

df_mags_coverage = pd.DataFrame(result_data)    
#df_mags_coverage
#df_mags_coverage.to_csv("/Users/mac/Desktop/Fu Group/1--Project1/2--DAG3_HGT/1--Recent_HGT/HGT_events/species_Exist.csv")

table_file2 = "/Users/mac/Desktop/Fu Group/1--Project1/2--DAG3_HGT/1--Recent_HGT/HGT_events/HGT_split_loose.csv"

merged_loose = pd.read_csv(table_file2, sep=",")
merged_loose


#result_df.set_index('HGT_ID', inplace=True)
#merged_loose.set_index('HGT_ID', inplace=True)

#validation_df.reset_index(inplace=True)
#result_df.reset_index(inplace=True)
#validation_df.to_csv("/Users/mac/Desktop/Fu Group/1--Project1/2--DAG3_HGT/1--Recent_HGT/HGT_events/split_loose_validate.csv")
#validation_df

#### HGT region coverage caculation
df_hgt_coverage = pd.read_csv("/Users/mac/Desktop/Fu Group/1--Project1/2--DAG3_HGT/1--Recent_HGT/HGT_events/HGT_cover_fraction.csv")

def check_values(value):
    num1, num2 = map(int, value.split('_'))
    if num1 >= 90 or num2 >= 90:
        return '2_2'
    else:
        return '0_0'

df_hgt_coverage.iloc[:, 1:] = df_hgt_coverage.iloc[:, 1:].applymap(check_values)

#df_hgt_coverage
def validate_values(loose_val, result_val):
    if loose_val == '2_2':
        if result_val == '2_2':
            return loose_val
        elif result_val == '2_0':
            return result_val
        elif result_val == '0_2':
            return result_val
        else:
            return '0_0'
    elif loose_val == '2_0':
        return loose_val if result_val.startswith('2_') else '0_0'
    elif loose_val == '0_2':
        return loose_val if result_val.endswith('_2') else '0_0'
    else:
        return loose_val


validation_df = merged_loose.copy()

for col in merged_loose.columns:
    validation_df[col] = [validate_values(loose_val, df_hgt_coverage.at[hgt_id, col]) 
                          for hgt_id, loose_val in merged_loose[col].items()]



def validate_values1(loose_val, result_val):
    if loose_val == '2_2':
        if result_val == '2_2':
            return loose_val
        elif result_val == '2_0':
            return "2_NA"
        elif result_val == '0_2':
            return "NA_2"
        else:
            return 'NA_NA'
    elif loose_val == '2_0':
        if result_val == '2_2':
            return "2_0"
        elif result_val == '2_0':
            return "2_NA"
        elif result_val == '0_2':
            return "NA_0"
        else:
            return 'NA_NA'
    elif loose_val == '0_2':
        if result_val == '2_2':
            return "0_2"
        elif result_val == '2_0':
            return "0_NA"
        elif result_val == '0_2':
            return "NA_2"
        else:
            return 'NA_NA'
    else:
        if result_val == '2_2':
            return "0_0"
        elif result_val == '2_0':
            return "0_NA"
        elif result_val == '0_2':
            return "NA_0"
        elif result_val == '0_0':
            return 'NA_NA'
        else:
            return loose_val

validation_df1 = validation_df.copy()

for col in validation_df.columns:
    validation_df1[col] = [validate_values1(loose_val, df_mags_coverage.at[hgt_id, col])
                          for hgt_id, loose_val in validation_df[col].items()]

## new species abun revise
data_cols = validation_df1.columns[validation_df1.columns != 'HGT_ID']


unique_values = validation_df1[data_cols].stack().unique()


counts_df = pd.DataFrame()

for value in unique_values:
    counts_df[value] = validation_df1[data_cols].apply(lambda row: sum(row.values == value), axis=1)

# Check if 'HGT_ID' exists in validation_df1 before inserting
if 'HGT_ID' in validation_df1.columns:
    counts_df.insert(0, 'HGT_ID', validation_df1['HGT_ID'])

#counts_df.to_csv("/Users/mac/Desktop/Fu Group/1--Project1/2--DAG3_HGT/1--Recent_HGT/read_split/count_loose_validate.csv")
condition1 = (counts_df['2_NA'] + counts_df['2_0']) > 0
condition2 = (counts_df['NA_2'] + counts_df['0_2']) > 0
condition3 = counts_df['2_2'] > 0

filtered_df = counts_df.loc[(condition1 & condition2) | condition3]

hgt_list = filtered_df['HGT_ID'].tolist()
validation_df2 = validation_df1[validation_df1["HGT_ID"].isin(hgt_list)]
df_abun1 = df_abun[df_abun["HGT_ID"].isin(hgt_list)]

validation_df2.to_csv("/Users/mac/Desktop/Fu Group/1--Project1/2--DAG3_HGT/1--Recent_HGT/HGT_events/HGT_Present_Table_Filter_New.csv")
df_abun1.to_csv("/Users/mac/Desktop/Fu Group/1--Project1/2--DAG3_HGT/1--Recent_HGT/HGT_events/HGT_events_Filter.csv")

HGT_genome_pair = df_abun1

HGT_genome_pair['mag_pair'] = HGT_genome_pair.apply(lambda row: '-'.join(sorted([row['MAG 1'], row['MAG 2']])), axis=1)

all_conserve = pd.read_csv("/Users/mac/Desktop/Fu Group/1--Project1/2--DAG3_HGT/1--Recent_HGT/ribo_protein/result/all_genes_long_100.csv")

average_per_genome_pair = all_conserve.groupby('genome_pair')['number'].mean().reset_index(name='average_number')

genomes_Pair_HGT_filter = HGT_genome_pair[~HGT_genome_pair['mag_pair'].isin(average_per_genome_pair['genome_pair'])]

genomes_Pair_HGT_filter.to_csv("/Users/mac/Desktop/Fu Group/1--Project1/2--DAG3_HGT/1--Recent_HGT/HGT_events/HGT_events_Filter_conserve_genomes_pair.csv", index=False)

#species_mag_number_filter = all_conserve.groupby('species_pair')['genome_pair'].nunique().reset_index(name='number_mags_pair')

#species_mag_number_filter.to_csv("/Users/mac/Desktop/Fu Group/1--Project1/2--DAG3_HGT/1--Recent_HGT/Tree_HGT_MAGs/species_pair_filter_number.csv", index=False)

