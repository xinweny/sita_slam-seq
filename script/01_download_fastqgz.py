#### Packages ####
import argparse
import os
import re
import pandas as pd
import glob
from requests import get
from bs4 import BeautifulSoup

#### Functions ####
def scrape_sample_names(gse):
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}"
    response = get(url)

    html_soup = BeautifulSoup(response.text, 'html.parser')
    tds = [td.text for td in html_soup.findAll("td")]
    tables = html_soup.findAll('table')

    if 'NIH grant(s)' in tds:
        start = 19
        end = 21
    else:
        start = 18
        end = 20

    sample_table = tables[start:end]

    gsm_sample = {}

    for table in sample_table:
        for i, row in enumerate(table.findAll("tr")):
            cells = row.findAll("td")
            text = [cell.text for cell in cells]
            name = re.sub(r'[\/:]', "", text[0])

            gsm_sample[name] = text[1]

    return gsm_sample

def aspera_download(srr, paired=False):
    if paired == True:
        for i in [1, 2]:
            sys_value = os.system(f"ascp -QT -l 1000m -P33001 -i ~/.aspera/cli/etc/asperaweb_id_dsa.openssh \
                        era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/{srr[0:6]}/00{srr[-1]}/{srr}/{srr}_{i}.fastq.gz \
                        ./{srr}_{i}.fq.gz")

            if sys_value == 256:
                print("Retrying new ftp...")
                sys_value = os.system(f"ascp -QT -l 1000m -P33001 -i ~/.aspera/cli/etc/asperaweb_id_dsa.openssh \
                era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/{srr[0:6]}/0{srr[-2:]}/{srr}/{srr}_{i}.fastq.gz \
                ./{srr}_{i}.fq.gz")

            if sys_value == 256:
                print("Retrying new ftp...")
                os.system(f"ascp -QT -l 1000m -P33001 -i ~/.aspera/cli/etc/asperaweb_id_dsa.openssh \
                era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/{srr[0:6]}/{srr}/{srr}_{i}.fastq.gz \
                ./{srr}_{i}.fq.gz")
    else:
        sys_value = os.system(f"ascp -QT -l 1000m -P33001 -i ~/.aspera/cli/etc/asperaweb_id_dsa.openssh \
                    era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/{srr[0:6]}/00{srr[-1]}/{srr}/{srr}.fastq.gz \
                    ./{srr}.fq.gz")

        if sys_value == 256:
            print("Retrying new ftp...")
            os.system(f"ascp -QT -l 1000m -P33001 -i ~/.aspera/cli/etc/asperaweb_id_dsa.openssh \
                        era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/{srr[0:6]}/0{srr[-2:]}/{srr}/{srr}.fastq.gz \
                        ./{srr}.fq.gz")

        if sys_value == 256:
            print("Retrying new ftp...")
            os.system(f"ascp -QT -l 1000m -P33001 -i ~/.aspera/cli/etc/asperaweb_id_dsa.openssh \
            era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/{srr[0:6]}/{srr}/{srr}.fastq.gz \
            ./{srr}.fq.gz")

def samplename_from_metadata(gsm, metadata_df, columns):
    info = [str(metadata_df.loc[metadata_df['Sample Name'] == gsm, col].iloc[0]) for col in columns]
    name = re.sub(r'[\/:]', "", '_'.join(info))

    return name

def is_paired(metadata, col_name, element):
    return metadata.loc[metadata[col_name] == element, 'LibraryLayout'].iloc[0] == "PAIRED"

# Main
def main():
    # Parser
    parser = argparse.ArgumentParser(description='Download fq.gz files from SRA')
    parser.add_argument('-m', required=True, help='Path to SRA run table metadata.')
    parser.add_argument('-c', nargs='?', const='', required=True, help='Columns in metadata table for naming, separated by |. If given sample name on GSE page is preferred, leave as an empty string.')

    args = parser.parse_args()
    metadata_path = args.m
    columns = (args.c).split('|')

    if args.c == '':
        gsm_samples = scrape_sample_names(os.getcwd().split("/")[-1])

    metadata_df = pd.read_csv(metadata_path, header=0, sep=',')
    gsms = metadata_df['Sample Name'].unique()

    # Set up file architecture
    os.system("mkdir fastq")
    os.chdir("fastq")

    for i, gsm in enumerate(gsms):
        print(f"({i + 1}/{len(gsms)}) Processing {gsm}...")

        if len([file for file in glob.glob(f"{gsm}*.fq.gz")]) > 0:
            print(f"{gsm}*.fq.gz already downloaded. Skipping...")
            continue

        gsm_srrs = metadata_df.loc[metadata_df['Sample Name'] == gsm, 'Run']
        name = gsm_samples[gsm] if args.c == '' else samplename_from_metadata(gsm, metadata_df, columns)

        for srr in list(gsm_srrs):
            if len([file for file in glob.glob(f"{srr}*.fq.gz")]) > 0:
                print(f"{srr}*.fq.gz already downloaded. Skipping...")
                continue

            print(f"Downloading {srr}...")
            aspera_download(srr, paired=is_paired(metadata_df, 'Run', srr))

        if len(gsm_srrs) == 1:
            print(f"Renaming fastq files for {gsm}...")
            srr = gsm_srrs.iloc[0]

            if is_paired(metadata_df, 'Run', srr):
                for i in [1, 2]:
                    os.system(f"mv -v '{srr}_{i}.fq.gz' '{gsm}_{name}_{i}.fq.gz'")
            else:
                os.system(f"mv -v '{srr}.fq.gz' '{gsm}_{name}.fq.gz'")
        else:
            print(f"Merging technical runs and renaming fastq files for {gsm}...")
            if is_paired(metadata_df, 'Run', srr):
                for i in [1, 2]:
                    single_fqs = ' '.join(gsm_srrs + f"_{i}.fq.gz")
                    os.system(f"cat {single_fqs} > '{gsm}_{name}_{i}.fq.gz' && rm {single_fqs}")
            else:
                srr_fqs = ' '.join(gsm_srrs + '.fq.gz')
                os.system(f"cat {srr_fqs} > '{gsm}_{name}.fq.gz' && rm {srr_fqs}")

#### Execute code ####
if __name__ == "__main__":
    main()
