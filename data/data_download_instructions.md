```
# Data Directory

This directory contains the data structure for the SAP Metabolomics project, but **NO DATA FILES ARE STORED IN THIS REPOSITORY**.

## Download Required Data

All data files must be downloaded from AWS S3 using the following commands:

### Using wget:
```bash
wget https://sapmet.s3.us-east-2.amazonaws.com/SAP_Metabolomics/filename
```

### Using curl:
```bash
curl -O https://sapmet.s3.us-east-2.amazonaws.com/SAP_Metabolomics/filename
```

Replace `filename` with the specific file you need from the directory structure below.

## Directory Structure

After downloading, your data directory should look like this:

```
data/
├── raw/
│   ├── treated_TAMU_met/          # Treated Tx2911 and P850029, analyzed at TAMU
│   └── untreated_CU_met/          # Untreated SAP, analyzed at Clemson University
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

## File Listing

For a complete list of available files, contact the authors or check the AWS S3 bucket directly.

⚠️ **Important**: The scripts expect data files to be in these exact directory locations. Make sure to download files to the correct subdirectories.

## Questions?

For data access issues, contact: Arlyn Ackerman (aja294@cornell.edu)
```
