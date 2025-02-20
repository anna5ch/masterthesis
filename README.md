# Benchmarking Novel Isoform Discovery Tools: A Pipeline for Short- and Long-Read RNA-seq Data

This repository contains the benchmarking pipeline for various isoform discovery tools for long-read and short-read RNA-seq data as part of my master thesis. 
This thesis evaluates the performance of different methods in detecting novel isoform candidates, isolating the discovery step without filtering. 

## Overview
The tools included in this benchmarking are:
- Long-read RNA-seq tools: IsoSeq, Bambu, FLAIR, FLAMES, IsoQuant, StringTie Long
- Short-read RNA-seq tools: Scallop, StringTie2

## Key Findings
### Reference-based tools
Bambu demonstrated the best overall performance in runs with reference annotation, together with a logisitic regression approach as a filtering method.
### De novo tools 
Scallop performed best in detecting novel isoforms without available reference annotations in combination with CPM or predicted probability (logistic regression approach).

## Installation & Requirements
We used several conda environments for different tools. Each required conda environment has a corresponding YAML file listing all necessary packages and dependencies. You can create the environments using: \
`conda env create -f environment_name.yaml` \
Some tools may require specific installation steps. Please refer to their official documentation:
- [IsoSeq](https://github.com/PacificBiosciences/IsoSeq)
- [Bambu](https://github.com/GoekeLab/bambu)
- [IsoQuant](https://github.com/ablab/IsoQuant)
- [Scallop](https://github.com/Kingsford-Group/scallop)
- [StringTie2](https://github.com/gpertea/stringtie)
- [FLAIR](https://github.com/BrooksLabUCSC/flair)
- [FLAMES](https://github.com/mritchielab/FLAMES)


## Repository structure
```
📂 workflow
├── 📁 rules # Snakemake workflow rules
├── 📁 scripts # Analysis and evaluation scripts 
├── 📁 envs # Conda environment YAML files 
└── 📄 Snakefile # Main Snakemake workflow file
```

## Contact
For any questions, feel free to open an issue.
