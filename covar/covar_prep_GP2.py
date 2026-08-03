# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import pathlib

BASE_PATH = "/config/workspace/gp2_tier2_eu_release11"
TABULAR_PATH = f"{BASE_PATH}/clinical_data/master_key_release11_final_vwb.csv"
WGS_PATH = f"{BASE_PATH}/wgs/deepvariant_joint_calling"

df = pd.read_csv(TABULAR_PATH, low_memory=False)

# -------------------------------------------------------------------------
# STEP 1. Define Ancestry, case/ctrl Status and Age
# -------------------------------------------------------------------------

# Ancestry label 
df["ancestry"] = df["wgs_label"]

# Phenotype recode: 2 = PD case, 1 = control, NaN = ignore; Excluding genetically-enriched removes 399 cases and 321 ctrls (amongst Europeans with WGS data)
is_case = (
    df["baseline_GP2_phenotype"].isin(["PD", "Affected_PD"]) &
    ~df["study_type"].eq("Genetically Enriched")
)
is_ctrl = (
    df["baseline_GP2_phenotype"].isin(["Control", "Population Control", "Unaffected"]) &
    ~df["study_type"].eq("Genetically Enriched")
)

df["status"] = np.nan
df.loc[is_case, "status"] = 2
df.loc[is_ctrl, "status"] = 1

print(df.groupby("ancestry")["status"].value_counts(dropna=False))


# Age hierarchy: age_of_onset > age_at_diagnosis > age_at_death > age_at_last_follow_up > age_at_sample_collection
age_cols = ["age_of_onset", "age_at_diagnosis", "age_at_death", "age_at_last_follow_up", "age_at_sample_collection"]
df["AGE"] = np.nan
df["AGE_SOURCE"] = None

for col in age_cols:
    missing = df["AGE"].isna() & df[col].notna()
    df.loc[missing, "AGE"] = df.loc[missing, col]
    df.loc[missing, "AGE_SOURCE"] = col


# -------------------------------------------------------------------------
# STEP 2. Subset & Clean relevant columns
# -------------------------------------------------------------------------

# Subset for WGS = 1 and phenotype in [1,2]
WGS_df = df[(df["wgs"] == 1) & (df["status"].notna())].copy()


# Reorder (also removes non-mentioned) & Rename columns
WGS_df["FID"] = WGS_df["GP2ID"]

WGS_df = WGS_df.reindex(
    columns=[
        "FID",
        "GP2ID",
        "ancestry",
        "amppd_wgs",
        "AGE_SOURCE",
        "AGE",
        "biological_sex_for_qc",
        "status",
    ]
)

WGS_df.rename(
    columns={
        "GP2ID": "IID",
        "ancestry": "Ancestry",
        "amppd_wgs": "AMP_PD",
        "biological_sex_for_qc": "SEX",
        "status": "STATUS"
    },
    inplace=True,
)

# Recode sex
sex_mapping = {"Female": 2, "Male": 1}
WGS_df["SEX"] = WGS_df["SEX"].map(sex_mapping).fillna(0).astype(int) # 8 samples with unknown sex (which will be labelled as 0)
print(WGS_df["SEX"].value_counts(dropna=False))


# Convert Status type to integer (1/2)
WGS_df["STATUS"] = WGS_df["STATUS"].astype(int)
print(WGS_df["STATUS"].value_counts(dropna=False))


# Recode AMP_PD
WGS_df["AMP_PD"] = WGS_df["AMP_PD"].fillna(0).astype(int)


# Recode AGE_SOURCE
print(WGS_df.groupby("STATUS")["AGE_SOURCE"].value_counts(dropna=False)) # Pre-WGS prune stats for Age column
age_mapping = {"age_of_onset":1, "age_at_diagnosis":2, "age_at_death":3, "age_at_last_follow_up":4, "age_at_sample_collection":5} #s.t. SKAT doesn't throw non-numerical covar error
WGS_df["AGE_SOURCE"] = WGS_df["AGE_SOURCE"].map(age_mapping)


# -------------------------------------------------------------------------
# STEP 3. Create ancestry-specific covariate files
# -------------------------------------------------------------------------
ancestries = ["AFR", "AJ", "AMR", "CAS", "EAS", "EUR", "MDE", "SAS", "CAH"]

for ancestry in ancestries:
    print("────────────────────────────────────────────")
    print(f"=== PROCESSING: {ancestry} ===")

    # --- 3.1 Subset ancestry
    ancestry_df = WGS_df[WGS_df["Ancestry"] == ancestry].copy()
    ancestry_df.reset_index(drop=True, inplace=True)
    ancestry_df.drop(columns=["Ancestry"], inplace=True)
    print(f"Total WGS samples in {ancestry}: {len(ancestry_df):,}")

    
    # --- 3.2 Get PCs
    pcs_path = f"{WGS_PATH}/pcs/{ancestry}/{ancestry}_release11.eigenvec"
    pcs_df = pd.read_csv(pcs_path, sep="\t", usecols=range(1, 12))
    print(f"PCs loaded: {pcs_df.shape}")

    cov_df = pd.merge(ancestry_df, pcs_df, on="IID", how="inner") #Merge intersection of dfs
    print(f"Merged covariate table: {cov_df.shape}")


# --- 3.3 Relatedness filter (Greedy Algorithm: resolve relatedness conflicts by removing the least amount of people)
    rel_path = f"{WGS_PATH}/related_samples/{ancestry}/{ancestry}_release11.related"
    related_df = pd.read_csv(rel_path, sep=",", comment='#', header=None)
    
    samples_to_remove = set()
    
    # Use sorted tuples instead of sets to ensure deterministic iteration
    pairs = related_df.iloc[:, [1, 3]].values.astype(str)
    related_pairs = [tuple(sorted(pair)) for pair in pairs if pair[0] != pair[1]]
    print(f"Total related pairs found: {len(related_pairs)}")

    # Iteratively remove the sample with the most conflicts
    while related_pairs:
        # Count how many pairs each sample is currently involved in
        counts = {}
        for pair in related_pairs:
            for sample in pair:
                counts[sample] = counts.get(sample, 0) + 1
        
        # Break ties deterministically by sample ID
        worst_offender = max(counts, key=lambda s: (counts[s], s))
        samples_to_remove.add(worst_offender)

        # Remove any pair containing this offender
        related_pairs = [p for p in related_pairs if worst_offender not in p]
        
    cov_df = cov_df[~cov_df["IID"].isin(samples_to_remove)]
    print(f"Optimized removal list: {len(samples_to_remove)} individuals removed.")

    # --- 3.4 Get Stats
    print(f"\nMissing:\n{cov_df.isnull().sum()}")


    # --- 3.5 Save to file
    cov_df["STATUS"] = cov_df.pop("STATUS") # Move Status column to end
    cov_df.to_csv(f"{ancestry}_covar_GP2.txt", sep="\t", index=False, na_rep="NA")
    
