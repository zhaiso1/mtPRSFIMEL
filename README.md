# mtPRSFIMEL
Improving multi-trait polygenic risk score prediction using fine-mapping and ensemble learning

## Overview
mtPRSFIMEL package (Zhai et al., 2025) implements a novel multi-trait polygenic risk score (mtPRS) method, mtPRS-FIMEL. This method tackles key challenges in existing multi-trait PRS approaches by integrating information from three levels:
- variant-level: causal variants identified through fine-mapping analysis;
- trait-level: data from multiple genetically correlated traits;
- method-level: an ensemble learning framework that optimally combines complementary PRS methods.

## Installation

mtPRSFIMEL R package requires R version >= 4.0.3.

```
library(devtools)
devtools::install_github("yaowuliu/ACAT") # Need to install ACAT first through the Github
devtools::install_github("zhaiso1/mtPRSFIMEL")
```

## Usage

### Step 1: PolyPred

Update effect size estimates in disease GWAS summary statistics by conducting PolyPred with fine-mapping. Details can be found in **PolyPred** folder (**run_polypred_sim_example.sh**).

### Step 2: mtPRS-FIMEL

Run mtPRS-FIMEL with refined effect size estimates. Details can be found in **doc** folder for mtPRSFIMEL user manual (**mtPRSFIMEL_0.1.0.pdf**), and **vignettes** folder for a demo illustrating how to use our softwares (**README.Rmd**).

## References

Weissbrod, O., Kanai, M., Shi, H., et al., 2022. Leveraging fine-mapping and multipopulation training data to improve cross-population polygenic risk scores. Nature genetics, 54(4), 450-458.

Zhai, S., Guo, B., Wu, B., Mehrotra, D.V., and Shen, J., 2023. Integrating multiple traits for improving polygenic risk prediction in disease and pharmacogenomics GWAS. Briefings in Bioinformatics, 24(4), bbad181.

Zhai, S., Guo, B., and Shen, J., 2025. Improving multi-trait polygenic risk score prediction using fine-mapping and ensemble learning.
