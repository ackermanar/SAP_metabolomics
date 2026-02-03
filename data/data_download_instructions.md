```
# Data Directory

This directory contains the data structure for the SAP Metabolomics project, but **NO DATA FILES ARE STORED IN THIS REPOSITORY**.

## Download Required Data

All data files must be downloaded from cloud storage using the following commands:

### Option 1: Download Individual Files

For specific files, use wget or curl:

**Using wget:**
```bash
wget https://sapmet.s3.us-east-2.amazonaws.com/SAP_Metabolomics/filename
```

**Using curl:**
```bash
curl -O https://sapmet.s3.us-east-2.amazonaws.com/SAP_Metabolomics/filename
```

Replace `filename` with the specific file you need from the directory structure below.

**For example, to download the SAP met data:**

```bash
wget https://sapmet.s3.us-east-2.amazonaws.com/SAP_Metabolomics/data/raw/untreated_CU_met/RAW_RESULTS_Sorghum_Grain_810_Samples_1M_cutoff_08082020.xlsx
```

### Option 2: Download Entire Dataset via AWS CLI (Recommended)

For the complete dataset, use AWS CLI:

```bash
aws s3 sync s3://sapmet/SAP_Metabolomics/ ./SAP_Metabolomics
```

**AWS CLI Setup Required:** Install and configure AWS CLI first. See instructions at:
https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-welcome.html

*Note: No AWS credentials needed for this public dataset - use `aws configure` with dummy values if prompted.*

### Option 3: Download Entire Dataset via Google Cloud CLI

Alternatively, use Google Cloud Storage with gsutil:

```bash
gsutil -m cp -r \
  "gs://sapmet/SAP_Metabolomics" \
  .
```

**Google Cloud CLI Setup Required:** Install and configure gcloud CLI first. See instructions at:
https://cloud.google.com/storage/docs/discover-object-storage-gcloud

*Note: No authentication required for this public dataset.*

## Directory Structure

After downloading, your data directory should look like this:

```
data/
├── raw/
│   ├── treated_TAMU_met/          # Fusarium-inoculated samples (TAMU IMAC)
│   └── untreated_CU_met/          # Control samples (Clemson MUAL)
└── processed/
    ├── Fig1/                      # Race-based pathway enrichment data
    ├── Fig2/                      # Panicle structure analysis
    ├── Fig3/                      # Grain color pathway analysis
    ├── Fig4/                      # Disease phenotype correlations
    ├── Fig5/                      # PCA and clustering data
    ├── Fig6/                      # Network-level pathway analysis
    ├── Fig7/                      # Volcano plot data
    ├── Fig8/                      # Correlation matrices
    └── Supplemental/              # Additional datasets
```

### Raw Dataset Details

- **treated_TAMU_met/**: *Fusarium verticillioides*-inoculated samples analyzed at Texas A&M University Integrated Metabolomic Analysis Core (IMAC)
- **untreated_CU_met/**: Control samples from uninoculated field trials analyzed at Clemson University Multi-User Analytics and Metabolomics Lab (MUAL)

## Essential Files for Basic Analysis

Download these key files first to run the main analyses:

1. `processed/Fig1/SAP19_MET_V4_ImputeChecks_NoBlanks.csv` - Main metabolomics dataset
2. `processed/Supplemental/SAP_accession_metadata.csv` - Sample metadata
3. `raw/field_data/SGM_DATA.csv` - Disease phenotype data

## File Listing

For a complete list of available files, contact the authors or check the cloud storage buckets directly.

⚠️ **Important**: The scripts expect data files to be in these exact directory locations. Make sure to download files to the correct subdirectories.

## Questions?

For data access issues, contact: Arlyn Ackerman (aja294@cornell.edu)
```
