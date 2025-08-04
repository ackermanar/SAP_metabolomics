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
│   ├── data_preprocessing.R            # Data cleaning and normalization
│   ├── statistical_analysis.R         # T-tests and comparative statistics
│   └── visualization.R                # Figure generation scripts
├── data/
│   ├── SAP19_MET_V4_ImputeChecks_NoBlanks.csv    # Primary metabolomics dataset
│   ├── phenotype_data/                 # Morphological and disease phenotypes
│   └── processed/                      # Intermediate analysis files
├── outputs/
│   ├── pathway_enrichment/             # Mummichog analysis results
│   ├── figures/                        # Publication-ready plots
│   └── supplementary/                  # Additional analysis outputs
├── docs/
│   ├── methodology.md                  # Detailed analytical protocols
│   └── data_dictionary.md              # Variable definitions and units
└── README.md
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

### 2. Data Processing Workflow

1. **Quality Control**: Filter features, normalize data using QuantileNorm + LogNorm + AutoNorm
2. **Statistical Testing**: T-tests comparing each subpopulation against remainder of SAP
3. **Pathway Prediction**: Mummichog v2 algorithm with retention time integration
4. **Results Compilation**: Automated aggregation across all subpopulations

### 3. Metabolomic Network Analysis

- **Feature-Based Molecular Networking (FBMN)** for compound annotation
- **Network Annotation Propagation (NAP)** for unknown compound prediction
- **Partial correlation analysis** accounting for confounding variables
- **Pathway-level fold change analysis** for biological interpretation

## Dependencies

### R Packages
```r
# Core data manipulation and analysis
library(magrittr)
library(data.table)
library(tidyverse)

# Parallel processing
library(purrr)
library(furrr)
library(future)
library(parallel)

# Metabolomics analysis
library(MetaboAnalystR)
library(RJSONIO)

# Statistical modeling
library(fitdistrplus)
library(memoise)
```

### External Tools
- **MS-DIAL**: Spectral feature alignment and processing
- **GNPS (Global Natural Products Social Molecular Networking)**: Compound annotation
- **MetaboAnalyst 5.0**: Web-based metabolomics analysis platform

## Usage

### Basic Pathway Enrichment Analysis

```r
# Set up parallel processing
plan(multisession, workers = (detectCores()-2))

# Load and format data
met <- fread("SAP19_MET_V4_ImputeChecks_NoBlanks.csv")

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

### Customizing for Different Grouping Variables

Simply modify the `group_by()` statement to analyze different categorical variables:

```r
# For panicle structure analysis
sum <- met %>% group_by(PANICLE_STRUCTURE) %>% summarize(n())

# For grain color analysis  
sum <- met %>% group_by(GRAIN_COLOR) %>% summarize(n())

# For combined variables
sum <- met %>% group_by(RACE, PANICLE_STRUCTURE) %>% summarize(n())
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

## License

This project is licensed under the MIT License - see the LICENSE file for details.

