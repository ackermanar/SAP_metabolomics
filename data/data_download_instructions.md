```
# Data Directory

This directory contains the data structure for the SAP Metabolomics project, but **NO DATA FILES ARE STORED IN THIS REPOSITORY**.

## Download Required Data

All data files must be downloaded from AWS S3. You have several options:

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

### Option 2: Download Entire Dataset (Recommended)

For the complete dataset, use AWS CLI:

```bash
aws s3 sync s3://sapmet/SAP_Metabolomics/ ./SAP_Metabolomics
```

**AWS CLI Setup Required:** Install and configure AWS CLI first. See instructions at:
https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-welcome.html

*Note: No AWS credentials needed for this public dataset - use `aws configure` with dummy values if prompted.*

## Directory Structure

After downloading, your data directory should look like this:

```
data/
├── raw/
│   ├── treated_TAMU_met/          # Fusarium-inoculated samples
│   └── untreated_CU_met/          # Control samples
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

## Essential Files for Basic Analysis

If downloading individual files, get these key datasets first:

1. `processed/Fig1/SAP19_MET_V4_ImputeChecks_NoBlanks.csv` - Main metabolomics dataset
2. `processed/Supplemental/SAP_accession_metadata.csv` - Sample metadata
3. `raw/field_data/SGM_DATA.csv` - Disease phenotype data

## File Listing

For a complete list of available files, contact the authors or check the AWS S3 bucket directly.

⚠️ **Important**: The scripts expect data files to be in these exact directory locations. Make sure to download files to the correct subdirectories.

## Questions?

For data access issues, contact: Arlyn Ackerman (aja294@cornell.edu)
```
