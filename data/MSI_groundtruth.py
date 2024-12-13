#%%
import pandas as pd
import numpy as np

# get PCR data
pre_pcr_df = pd.read_csv(".../MSI_PCR.csv")
pcr_df= pd.DataFrame()
pcr_df["Patient ID"] = pre_pcr_df["Tumor"]
pcr_df["tissue"] = pre_pcr_df["Project"]
pcr_df["MSI status"] = pre_pcr_df["MSI_Class"]
pcr_df['PCR'] = pcr_df['MSI status'].apply(lambda x: np.nan if pd.isna(x) else ('MSIH' if x ==1 else 'nonMSIH'))
pcr_df.drop(columns=["MSI status"],inplace=True)

# get NGS tools data
pre_clinical_data_df = pd.read_csv(".../MSI_NGS.tsv", sep='\t', header=0)
clinical_data_df = pd.DataFrame()
clinical_data_df["Patient ID"] = pre_clinical_data_df["Patient ID"]
clinical_data_df["MSI MANTIS Score"] = pre_clinical_data_df["MSI MANTIS Score"]
clinical_data_df["MSIsensor Score"] = pre_clinical_data_df["MSIsensor Score"]
clinical_data_df["Mutation Count"] = pre_clinical_data_df["Mutation Count"]
clinical_data_df["TMB (nonsynonymous)"] = pre_clinical_data_df["TMB (nonsynonymous)"]
clinical_data_df["Tissue"] = pre_clinical_data_df["TCGA PanCanAtlas Cancer Type Acronym"]

# apply thresholds
clinical_data_df['MANTIS'] = clinical_data_df['MSI MANTIS Score'].apply(lambda x: np.nan if pd.isna(x) else ('MSIH' if x >= 0.4 else 'nonMSIH'))
clinical_data_df['MANTIS'] = clinical_data_df[['MANTIS']].copy()
clinical_data_df['MSISENS'] = clinical_data_df['MSIsensor Score'].apply(lambda x: np.nan if pd.isna(x) else ('MSIH' if x >= 3.5 else 'nonMSIH'))
clinical_data_df['MSISENS'] = clinical_data_df[['MSISENS']].copy()

merged_df = pd.merge(pcr_df, clinical_data_df, on='Patient ID', how='outer')

tissues = ["UCEC","STAD","COAD","READ"]
merged_df = merged_df[merged_df["Tissue"].isin(tissues)]

# get consensus
condition = (~merged_df['MANTIS'].isna()) & (merged_df['MANTIS'] != 'empty') & \
            (~merged_df['MSISENS'].isna()) & (merged_df['MSISENS'] != 'empty') & \
            (~merged_df['PCR'].isna()) & (merged_df['PCR'] != 'empty')
merged_df = merged_df[condition]

consesnus_df = merged_df[(merged_df['MANTIS'] == merged_df['MSISENS']) & (merged_df['MSISENS'] == merged_df['PCR'])]
consesnus_df.rename(columns={'Patient ID': 'bcr_patient_barcode'}, inplace=True)

consesnus_df = consesnus_df.sort_values(by=["bcr_patient_barcode"])
consesnus_df.reset_index(inplace=True, drop = True)

consesnus_df["site"] = [x[5:7] for x in consesnus_df["bcr_patient_barcode"].tolist()]

# get source sites
tss = pd.read_csv(".../TCGA_TSS.csv",sep=';')
tss = tss.rename(columns={"TSS Code":"site"})
matched_tss = consesnus_df.merge(tss, on="site", how="left")

# remove the fixed test sites from the dataset
consesnus_df=matched_tss
fixed_test_sites = ['ILSBIO','Roswell Park','University of North Carolina','Peter MacCallum Cancer Center']

fixed_test_df = consesnus_df[consesnus_df['Source Site'].isin(fixed_test_sites)]
remaining_df = consesnus_df[~consesnus_df['Source Site'].isin(fixed_test_sites)]

remaining_sites = remaining_df['Source Site'].unique()
np.random.shuffle(remaining_sites)

folds = np.array_split(remaining_sites, 5)

# function to assign sites to a fold
def assign_fold(fold_idx):
    print()
    train_sites = np.concatenate([folds[i] for i in range(5) if i != fold_idx])
    
    df = remaining_df.copy()
    df['split'] = df['Source Site'].apply(
        lambda x: 'train' if x in train_sites else 'validation')

    test_df = fixed_test_df.copy()
    test_df['split'] = 'test'
    df = pd.concat([df, test_df])
    return df

fold_dfs = [assign_fold(i) for i in range(5)]

for idx, fold_df in enumerate(fold_dfs):
    print(f"Fold {idx + 1} sizes:")
    print("Train set size:", len(fold_df[fold_df['split'] == 'train']))
    print("Validation set size:", len(fold_df[fold_df['split'] == 'validation']))
    print("Test set size:", len(fold_df[fold_df['split'] == 'test']))
