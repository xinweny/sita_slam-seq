#### Packages ####
import pandas as pd
import glob, os
import re

#### Load data ####
os.chdir("/Users/pomato/mrc/project/sita_slam-seq/")

tc_dir = "./processed/PROJ1624_TC"
files = glob.glob(f"{tc_dir}/*.csv")

tc_rates = []

# Calculate T > C conversion rate per sample
for file in files:
    # Get sample name
    sample = os.path.basename(file).replace("_R1_trimmed.fq_slamdunk_mapped_filtered_mutationrates_utr.csv", "")
    
    # Format grouping
    group = re.sub(r'^T[0-9]_', '', re.sub(r'_rep[0-9]$', '', sample))

    # Calculate overall T > C conversion rate
    df = pd.read_csv(file, sep="\t", header=2)
    tc = df['T_C'].sum() / (df['T_C'].sum() + df['T_T'].sum() + df['T_A'].sum() + df['T_G'].sum()) * 100
    
    tc_rates.append([sample, group, tc])

tc_df = pd.DataFrame(tc_rates, columns=['sample', 'group', 'TC'])

# Save output
tc_df.to_csv("./processed/PROJ1624_TC.txt", sep="\t", index=False)