Welcome to the GitHub repository for the following publication: [The energetic and allosteric landscape for KRAS inhibition (Weng C, Faure AJ & Lehner B, 2022)](https://www.biorxiv.org/content/10.1101/2022.12.06.519122v1)

Here you'll find an R package with all scripts to reproduce the figures and results from the computational analyses described in the paper.

# Table Of Contents

* **1. [Required Software](#required-software)**
* **2. [Required Data](#required-data)**
* **3. [Installation Instructions](#installation-instructions)**
* **4. [Usage](#usage)**

# Required Software

To run the krasddpcams pipeline you will need the following software and associated packages:

* **[_R_](https://www.r-project.org/) >=v3.6.1** (ggplot2, ggpubr, ggrepel, ROCR, bio3d, GGally, plot3D, Cairo, ggstatsplot, openxlsx, data.table, dplyr)

# Required Data

Fitness scores, inferred free energy changes and required miscellaneous files should be downloaded from **[here](https://crgcnag-my.sharepoint.com/:f:/g/personal/cweng_crg_es/EliX349TTkpIoMomBwphyRMBYI17nEt4XZ45XcTvWtpuyw)** and unzipped in your project directory (see 'base_dir' option) i.e. where output files should be written.

The following software is optional:

* **[DiMSum](https://github.com/lehner-lab/DiMSum) v1.2.9** (pipeline for pre-processing deep mutational scanning data i.e. FASTQ to fitness)
* **[MoCHI](https://github.com/lehner-lab/MoCHI)** (tool to fit mechanistic models to deep mutational scanning data i.e. fitness to free energy changes)

Configuration files and additional scripts used for running DiMSum and MoCHI are also available **[here](https://crgcnag-my.sharepoint.com/:f:/g/personal/cweng_crg_es/EliX349TTkpIoMomBwphyRMBYI17nEt4XZ45XcTvWtpuyw)**.

# Installation Instructions

Make sure you have git and conda installed and then run (expected install time <1h):

```
# Install dependencies (preferably in a fresh conda environment)
conda install -c conda-forge r-base>=3.6.1 r-ggplot2 r-ggpubr r-ggrepel r-rocr r-bio3d r-ggally r-plot3d r-cairo r-ggstatsplot r-openxlsx r-data.table r-dplyr

# Open an R session and install the krasddpcams R package
if(!require(devtools)) install.packages("devtools")
devtools::install_github("lehner-lab/krasddpcams")
```

# Usage

The top-level function **krasddpcams()** is the recommended entry point to the pipeline and by default reproduces the figures and results from the computational analyses described in the following publication: [The energetic and allosteric landscape for KRAS inhibition (Weng C, Faure AJ & Lehner B, 2022)](https://www.biorxiv.org/content/10.1101/2022.12.06.519122v1). See [Required Data](#required-data) for instructions on how to obtain all required data and miscellaneous files before running the pipeline. Expected run time <20min.

```
library(krasddpcams)
krasddpcams(base_dir = "MY_PROJECT_DIRECTORY")
```
