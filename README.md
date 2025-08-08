# SAP Metabolomics: Organizational Trends of the Grain Metabolome in the Sorghum Association Panel

## Project Overview

This repository contains code and analysis pipelines for **"The Organizational Trends of the Grain Metabolome in the Sorghum Association Panel"** by Ackerman et al. The study applies non-targeted LC-MS metabolomics to understand the biochemical architecture of sorghum grain across 407 diverse accessions from the Sorghum Association Panel (SAP), with special focus on metabolic responses to *Fusarium* grain mold infection.

### Research Highlights

- **Comprehensive metabolomic profiling** of 407 sorghum accessions representing global genetic diversity
- **Pathway enrichment analysis** across morphological subpopulations (race, panicle structure, grain color)
- **Stress response characterization** using *Fusarium verticillioides* grain mold as a model pathosystem
- **Network-level metabolic analysis** revealing constitutive vs. induced defense mechanisms
- **Biomarker identification** for disease resistance and metabolomic-assisted breeding applications

## Repository Structure

```
SAP_metabolomics/
├── scripts/
│   ├── enrichment_analysis.R          # Main mummichog pathway enrichment pipeline
│   ├── SAP__Met_Figures.R             # Publication figure generation pipeline
├── data/                              # DATA NOT INCLUDED - Download from AWS
│   ├── raw/
│   │   ├── treated_TAMU_met/          # Fusarium-inoculated Tx2911 and P850029 samples (TAMU analysis)
│   │   └── untreated_CU_met/          # Entire SAP (Clemson analysis)
│   └── processed/
│       ├── Fig1/                      # Race-based pathway enrichment data
│       ├── Fig2/                      # Panicle structure metabolomic profiles
│       ├── Fig3/                      # Grain color pathway analysis
│       ├── Fig4/                      # Disease phenotype correlations
│       ├── Fig5/                      # PCA and metabolomic clustering
│       ├── Fig6/                      # Network-level pathway analysis
│       ├── Fig7/                      # Volcano plot differential metabolites
│       ├── Fig8/                      # Partial correlation heatmaps
│       └── Supplemental/              # Additional datasets and analyses
└── README.md
```

## Data Access

⚠️ **Important**: Data files are **NOT** included in this GitHub repository due to size constraints. All datasets must be downloaded separately from AWS S3.

### Download Instructions

You have several options for downloading the required data:

#### Option 1: Download Individual Files

For specific files, use wget or curl:

```bash
# Using wget
wget https://sapmet.s3.us-east-2.amazonaws.com/SAP_Metabolomics/filename

# Using curl
curl -O https://sapmet.s3.us-east-2.amazonaws.com/SAP_Metabolomics/filename
```

Replace `filename` with the specific file you need.

#### Option 2: Download Entire Dataset (Recommended)

For the complete dataset, use AWS CLI:

```bash
aws s3 sync s3://sapmet/SAP_Metabolomics/ ./SAP_Metabolomics
```

**AWS CLI Setup Required**: Install and configure AWS CLI first. See instructions at:
https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-welcome.html

*Note: No AWS credentials needed for this public dataset - use `aws configure` with dummy values if prompted.*

#### Create Directory Structure

```bash
# Create directory structure
mkdir -p SAP_Metabolomics/data/raw/treated_TAMU_met
mkdir -p SAP_Metabolomics/data/raw/untreated_CU_met
mkdir -p SAP_Metabolomics/data/processed/Fig1
mkdir -p SAP_Metabolomics/data/processed/Fig2
mkdir -p SAP_Metabolomics/data/processed/Fig3
mkdir -p SAP_Metabolomics/data/processed/Fig4
mkdir -p SAP_Metabolomics/data/processed/Fig5
mkdir -p SAP_Metabolomics/data/processed/Fig6
mkdir -p SAP_Metabolomics/data/processed/Fig7
mkdir -p SAP_Metabolomics/data/processed/Fig8
mkdir -p SAP_Metabolomics/data/processed/Supplemental
```

### Data Organization

**Raw Datasets (`data/raw/`):**
- `treated_TAMU_met/`: *Fusarium verticillioides*-inoculated samples analyzed at Texas A&M University Integrated Metabolomic Analysis Core (IMAC) 
- `untreated_CU_met/`: Control samples from uninoculated field trials analyzed at Clemson University Multi-User Analytics and Metabolomics Lab (MUAL)

**Processed Data by Figure (`data/processed/`):**
- `Fig1/`: Pathway enrichment results across sorghum races and check lines (Tx2911, P850029)
- `Fig2/`: Metabolomic pathway analysis by panicle structure (compact, semi-compact, open)
- `Fig3/`: Grain color-associated metabolic pathway enrichment
- `Fig4/`: Disease phenotype data (PGMSR, fumonisin, FSDI) with metabolomic correlations
- `Fig5/`: Principal component analysis of treated vs. untreated metabolomes
- `Fig6/`: Network-level pathway fold-change analysis between genotypes
- `Fig7/`: Volcano plot data for differential metabolite abundance
- `Fig8/`: Partial correlation matrices with phenotype associations
- `Supplemental/`: Extended datasets, accession metadata, and additional analyses

### Essential Files for Basic Analysis

If downloading individual files, get these key datasets first:

```bash
# Main metabolomics dataset
wget https://sapmet.s3.us-east-2.amazonaws.com/SAP_Metabolomics/processed/Fig1/SAP19_MET_V4_ImputeChecks_NoBlanks.csv

# Sample metadata
wget https://sapmet.s3.us-east-2.amazonaws.com/SAP_Metabolomics/processed/Supplemental/SAP_accession_metadata.csv

# Disease phenotype data
wget https://sapmet.s3.us-east-2.amazonaws.com/SAP_Metabolomics/raw/field_data/SGM_DATA.csv
```

## Key Analysis Pipeline

### 1. Pathway Enrichment Analysis (`enrichment_analysis.R`)

The core analytical framework uses the **mummichog algorithm** through MetaboAnalystR to predict metabolic pathway activity without requiring upfront metabolite identification. This approach is particularly powerful for non-targeted metabolomics datasets where thousands of features need network-level interpretation.

**Pipeline Features:**
- **Parallel processing** for computational efficiency across multiple subpopulations
- **Zea mays KEGG reference** for pathway mapping (closest available plant model)
- **Statistical rigor** with t-tests between target subpopulations and remainder of SAP
- **Flexible grouping** easily adaptable for different categorical variables (race, panicle structure, grain color)

**Key Functions:**
```r
# Main analysis function with parallel processing
met_ttest <- function(peak_list, job)

# Safety wrapper for robust execution
met_ttest_safe <- safely(met_ttest, quiet = TRUE)

# Automated pipeline across all subpopulations
map2(met_list, jobs, ~met_ttest_safe(.x, .y))
```

### 2. Publication Figure Generation (`SAP__Met_Figures.R`)

Comprehensive visualization pipeline for creating all manuscript figures with publication-ready formatting and consistent aesthetics.

**Figure Generation Capabilities:**
- **Figure 1**: Race-based pathway enrichment bubble plots with sample size annotations
- **Figure 2**: Panicle structure metabolomic pathway analysis with color-coded significance
- **Figure 3**: Grain color pathway enrichment visualization
- **Figure 4**: Disease phenotype scatter plots (PGMSR vs. fumonisin with FSDI sizing)
- **Figure 5**: Principal Component Analysis plots with statistical ellipses
- **Figure 6**: Network-level pathway bar charts with fold-change directionality
- **Figure 7**: Enhanced volcano plots with metabolite class annotations and selective labeling
- **Figure 8**: Correlation heatmaps with phenotype annotations and hierarchical clustering
- **Supplementary Figure 1**: Multi-panel morphological trait distribution plots

**Key Visualization Features:**
```r
# Automated sample size annotation for bubble plots
merged_counts$Label <- paste(merged_counts$RACE, "\n (n = ", merged_counts$Count, ")", sep = "")

# Publication-ready color schemes
scale_color_gradient2(low = "lightblue1", mid = "lightskyblue", high = "dodgerblue4")

# Enhanced volcano plots with selective metabolite labeling
EnhancedVolcano(selectLab = c('Luteolin','Apigeninidin', 'Caffeic_Acid', 'Salicylic_acid'))

# Correlation heatmaps with custom annotations
pheatmap(annotation_col = metadata, annotation_row = rowData, cluster_cols = FALSE)
```

**Output Management:**
- **Multiple resolution options**: Small (7"×6") and large (14"×12") format exports
- **High-quality TIFF export**: 300 DPI publication-ready figures
- **Consistent theming**: Standardized fonts, colors, and layouts across all figures
- **Automated file naming**: Systematic figure file organization

### 3. Data Processing Workflow

1. **Quality Control**: Filter features, normalize data using QuantileNorm + LogNorm + AutoNorm
2. **Statistical Testing**: T-tests comparing each subpopulation against remainder of SAP
3. **Pathway Prediction**: Mummichog v2 algorithm with retention time integration
4. **Results Compilation**: Automated aggregation across all subpopulations

### 4. Metabolomic Network Analysis

- **Feature-Based Molecular Networking (FBMN)** for compound annotation
- **Network Annotation Propagation (NAP)** for unknown compound prediction
- **Partial correlation analysis** accounting for confounding variables
- **Pathway-level fold change analysis** for biological interpretation

## Dependencies

### R Packages
```r
# Core data manipulation and analysis
library(data.table)
library(tidyverse)
library(magrittr)

# Parallel processing
library(purrr)
library(furrr)
library(future)
library(parallel)

# Metabolomics analysis
library(MetaboAnalystR)
library(RJSONIO)
library(fitdistrplus)
library(memoise)

# Visualization and figure generation
library(pheatmap)
library(grid)
library(gridExtra)
library(ggforce)
library(EnhancedVolcano)

# Statistical analysis
library(ppcor)

# Color palettes and themes
library(wesanderson)
library(RColorBrewer)
library(viridis)
```

### External Tools
- **MS-DIAL**: Spectral feature alignment and processing
- **GNPS (Global Natural Products Social Molecular Networking)**: Compound annotation
- **MetaboAnalyst 5.0**: Web-based metabolomics analysis platform

## Usage

### Setup and Data Download

```bash
# Clone repository
git clone https://github.com/username/SAP_metabolomics.git
cd SAP_metabolomics

# Option A: Download entire dataset (recommended)
aws s3 sync s3://sapmet/SAP_Metabolomics/ ./SAP_Metabolomics

# Option B: Create structure and download individual files
mkdir -p data/raw/{treated_TAMU_met,untreated_CU_met}
mkdir -p data/processed/{Fig1,Fig2,Fig3,Fig4,Fig5,Fig6,Fig7,Fig8,Supplemental}

# Download key files individually
wget https://sapmet.s3.us-east-2.amazonaws.com/SAP_Metabolomics/processed/Fig1/SAP19_MET_V4_ImputeChecks_NoBlanks.csv \
     -P data/processed/Fig1/
```

### Generating Publication Figures

```r
# Load figure generation script
source("scripts/SAP__Met_Figures.R")

# Set working directory
setwd("SAP_Metabolomics")

# Generate individual figures (examples)
# Figure 1: Race-based pathway enrichment
# (Automatically reads data/processed/Fig1/final_peak_list_RACE.csv)

# Figure 4: Disease phenotype correlations
# (Reads data/raw/field_data/SGM_DATA.csv)

# Figure 7: Volcano plot with enhanced annotations
# (Uses data/processed/Fig7/Data/Tx2911upreg.csv)

# All figures are automatically saved as high-resolution TIFF files
```

### Basic Pathway Enrichment Analysis

```r
# Set up parallel processing
plan(multisession, workers = (detectCores()-2))

# Load and format data
met <- fread("data/processed/Fig1/SAP19_MET_V4_ImputeChecks_NoBlanks.csv")

# Define subpopulations (minimum n=5)
sum <- met %>% 
  group_by(RACE) %>% 
  summarize(n()) %>%
  filter(ncount >= 5 & !is.na(RACE))

# Execute parallel analysis
peak_lists <- map2(met_list, jobs, ~met_ttest_safe(.x, .y))

# Compile results
final_results <- peak_lists %>% 
  future_map("result") %>% 
  compact() %>% 
  rbindlist()
```

### Customizing Visualizations

```r
# Modify color schemes for different analyses
scale_color_manual(values = c("Tx2911" = "darkorange2", "P850029" = "dodgerblue2"))

# Adjust figure dimensions for different layouts
ggsave("custom_figure.tiff", device = "tiff", width = 10, height = 8, dpi = 300)

# Customize heatmap annotations
ann_colors <- list(
    Label = c("P85" = "dodgerblue2", "Tx2911" = "darkorange2"),
    P.value = c("Significant" = "chartreuse2", "Insignificant" = "coral2")
)
```

## Data Requirements

### Input Format
- **Rows**: Individual plant samples with metadata
- **Columns**: Metabolomic features (m/z/retention time) plus categorical variables
- **Grouping Variables**: Categorical columns for subpopulation analysis (RACE, PANICLE_STRUCTURE, GRAIN_COLOR, etc.)
- **Feature Format**: "m.z/retention_time" for mummichog compatibility

### Sample Metadata Requirements
- Biological replicates across genotypes
- Standardized categorical classifications
- Quality control samples interspersed (recommended every 10-15 injections)

## Statistical Considerations

### Power and Replication
- **Minimum group size**: n≥5 for pathway enrichment analysis
- **Biological replication**: Multiple plants per genotype recommended
- **Technical replication**: QC samples for measurement stability assessment

### Multiple Testing
- **Mummichog approach**: Inherently accounts for multiple testing through pathway-level analysis
- **P-value threshold**: 0.05 with 5ppm mass tolerance
- **Effect size**: log₁₀ EASE scores for pathway enrichment magnitude

## Publication Details

**Title**: The Organizational Trends of the Grain Metabolome in the Sorghum Association Panel

**Authors**: Arlyn Ackerman¹, Maria A Conti², Andrew Disharoon², Anthony Wenndt³, William Caughman², Richard Boyles²

**Affiliations**:
- ¹Breeding Insight, Cornell University, Ithaca, NY, USA
- ²Plant and Environmental Sciences, PeeDee Research and Education Center, Florence, South Carolina, USA  
- ³Global Alliance for Improved Nutrition (GAIN), Washington, DC, USA

**Funding**: Foundation for Food and Agriculture Research, project CA21-SS-0000000061

## Citation

If you use this code or methodology, please cite:

```
Ackerman, A., Conti, M.A., Disharoon, A., Wenndt, A., Caughman, W., & Boyles, R. (2025). 
The Organizational Trends of the Grain Metabolome in the Sorghum Association Panel. 
[Journal Name], [Volume], [Pages]. doi: [DOI]
```

## Contributing

We welcome contributions to improve the analytical pipeline. Please:

1. Fork the repository
2. Create a feature branch
3. Make your changes with clear documentation
4. Submit a pull request with detailed description

## Support

For questions about the methodology or code implementation:

- **Primary Contact**: Arlyn Ackerman (aja294@cornell.edu)
- **Issues**: Submit via GitHub Issues
- **Methodology Questions**: See `docs/methodology.md` for detailed protocols

---
