#%%
import pandas as pd
import numpy as np

# get HRD scores
pre_df_hrd = pd.read_excel(".../HRD_scar_score.xlsx")
df_hrd = pd.DataFrame()
df_hrd["Patient ID"] = pre_df_hrd["PATIENT"]
df_hrd["Study ID"] = pre_df_hrd["Study ID"]
df_hrd["Site"] = pre_df_hrd["SITE"]
df_hrd["HRD_TAI_paper"] = pre_df_hrd["HRD_TAI_paper"]
df_hrd["HRD_LST_paper"] = pre_df_hrd["HRD_LST_paper"]
df_hrd["HRD_LOH_paper"] = pre_df_hrd["HRD_LOH_paper"]
df_hrd["HRD_Sum_paper"] = pre_df_hrd["HRD_Sum_paper"]
df_hrd["HRD_binary_paper"] = pre_df_hrd["HRD_binary_paper"]

df_hrd = df_hrd[df_hrd['HRD_binary_paper'].notna()]

study_id = df_hrd["Study ID"].tolist()
tissues_raw = [x.split("_")[0].upper() for x in study_id]
df_hrd["tissue"] = tissues_raw

# binarize at thresholds
cutoffs = {"BRCA":45,"OV":54,"PAAD":37,"PRAD":21}
tissues = ["BRCA","OV","PAAD","PRAD"]
df_hrd = df_hrd[df_hrd['tissue'].isin(tissues)]
df_hrd['cutoff'] = df_hrd['tissue'].map(cutoffs)

if df_hrd['cutoff'].isnull().any():
    missing_types = df_hrd[df_hrd['cutoff'].isnull()]['tissue'].unique()
    print(f"Warning: No cutoff defined for cancer types: {missing_types}")

df_hrd['HRD_status'] = np.where(df_hrd['HRD_Sum_paper'] > df_hrd['cutoff'], 'HRD_positive', 'HRD_negative')

# get source sites
tss = pd.read_csv(".../TCGA_TSS.csv",sep=';')
tss = tss.rename(columns={"TSS Code":"Site"})
matched_tss = df_hrd.merge(tss, on="Site", how="left")

# remove the fixed test sites from the dataset
consesnus_df=matched_tss
fixed_test_sites = ['BC Cancer Agency',
 'Baylor College of Medicine',
 "Brigham and Women's Hospital",
 'CHI-Penrose Colorado',
 'Emory University',
 'Essen',
 'Fondazione-Besta',
 'Fred Hutchinson',
 'Fundacio Clinic per a la Recerca Biomedica',
 'Garvan Institute of Medical Research',
 'Memorial Sloan Kettering Cancer Center',
 'Gundersen Lutheran Health System',
 'Heidelberg',
 'Hospital Louis Pradel',
 'ILSBIO',
 'NCI',
 'Roswell Park',
 'The University of New South Wales',
 'Translational Genomics Research Institute',
 'University of Liverpool',
 'University of North Carolina']

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
        lambda x: 'train' if x in train_sites else 'validation'
    )
    test_df = fixed_test_df.copy()
    test_df['split'] = 'test'
    df = pd.concat([df, test_df])
    return df

fold_dfs = [assign_fold(i) for i in range(5)]