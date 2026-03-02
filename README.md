# dna_methylation_smk  [![zread](https://img.shields.io/badge/Ask_Zread-_.svg?style=flat&color=00b0aa&labelColor=000000&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB3aWR0aD0iMTYiIGhlaWdodD0iMTYiIHZpZXdCb3g9IjAgMCAxNiAxNiIgZmlsbD0ibm9uZSIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KPHBhdGggZD0iTTQuOTYxNTYgMS42MDAxSDIuMjQxNTZDMS44ODgxIDEuNjAwMSAxLjYwMTU2IDEuODg2NjQgMS42MDE1NiAyLjI0MDFWNC45NjAxQzEuNjAxNTYgNS4zMTM1NiAxLjg4ODEgNS42MDAxIDIuMjQxNTYgNS42MDAxSDQuOTYxNTZDNS4zMTUwMiA1LjYwMDEgNS42MDE1NiA1LjMxMzU2IDUuNjAxNTYgNC45NjAxVjIuMjQwMUM1LjYwMTU2IDEuODg2NjQgNS4zMTUwMiAxLjYwMDEgNC45NjE1NiAxLjYwMDFaIiBmaWxsPSIjZmZmIi8%2BCjxwYXRoIGQ9Ik00Ljk2MTU2IDEwLjM5OTlIMi4yNDE1NkMxLjg4ODEgMTAuMzk5OSAxLjYwMTU2IDEwLjY4NjQgMS42MDE1NiAxMS4wMzk5VjEzLjc1OTlDMS42MDE1NiAxNC4xMTM0IDEuODg4MSAxNC4zOTk5IDIuMjQxNTYgMTQuMzk5OUg0Ljk2MTU2QzUuMzE1MDIgMTQuMzk5OSA1LjYwMTU2IDE0LjExMzQgNS42MDE1NiAxMy43NTk5VjExLjAzOTlDNS42MDE1NiAxMC42ODY0IDUuMzE1MDIgMTAuMzk5OSA0Ljk2MTU2IDEwLjM5OTlaIiBmaWxsPSIjZmZmIi8%2BCjxwYXRoIGQ9Ik0xMy43NTg0IDEuNjAwMUgxMS4wMzg0QzEwLjY4NSAxLjYwMDEgMTAuMzk4NCAxLjg4NjY0IDEwLjM5ODQgMi4yNDAxVjQuOTYwMUMxMC4zOTg0IDUuMzEzNTYgMTAuNjg1IDUuNjAwMSAxMS4wMzg0IDUuNjAwMUgxMy43NTg0QzE0LjExMTkgNS42MDAxIDE0LjM5ODQgNS4zMTM1NiAxNC4zOTg0IDQuOTYwMVYyLjI0MDFDMTQuMzk4NCAxLjg4NjY0IDE0LjExMTkgMS42MDAxIDEzLjc1ODQgMS42MDAxWiIgZmlsbD0iI2ZmZiIvPgo8cGF0aCBkPSJNNCAxMkwxMiA0TDQgMTJaIiBmaWxsPSIjZmZmIi8%2BCjxwYXRoIGQ9Ik00IDEyTDEyIDQiIHN0cm9rZT0iI2ZmZiIgc3Ryb2tlLXdpZHRoPSIxLjUiIHN0cm9rZS1saW5lY2FwPSJyb3VuZCIvPgo8L3N2Zz4K&logoColor=ffffff)](https://zread.ai/YuanfengZhang/dna_methylation_smk)
Snakemake pipelines to run NGS-based methylation analysis and benchmark. [check the DeepWiki explaination](https://deepwiki.com/YuanfengZhang/dna_methylation_smk)
The python script and .ipynb files to reproduce the benchmarking results, including evaluation and statistical analysis, are in the `benchmark` folder.

The scripts for data visualization are listed in `benchmark/figures` folder.

Following is a simplified schematic diagram involving the dna_methylation_smk and MethCali:
![Schematic diagram](MISC/overview.png)

Here is the way to calculate the RMSE and SpearmanR:
![RMSE](MISC/RMSE.png)
![SpearmanR](MISC/SpearmanR.png)

Please read the `benchmark/README.md` for more details.

## Requirements
While the snakemake and python scripts are compatible with most operating systems, there are many bioinformatic tools evaluated here only work on x64 Linux. Please use a x64 Linux server / container with at least 128GB RAM and 8 CPU cores to run the pipelines.

## Installation
1. Install [miniforge](https://github.com/conda-forge/miniforge) following the instructions on the github page.
2. Create a conda environment with [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) following the instructions.
3. Clone this git repo with submodules:
    ```bash
    git clone https://github.com/YuanfengZhang/dna_methylation_smk --recurse-submodules
    ```
4. Install and compile the third-party tools: abismal, BatMeth2, bsgenova, FAME, hisat-3n, Whisper, BioSeqZip, and rastair.
5. Modify the `config/runtime_config.yaml` to set the paths of the files and parameters required by third-party tools.
6. Fill the sample sheet csv file according to the format described in `utils/sample_sheet_parser.py` and `config/sample_sheet.csv`.
7. Run the snakemake pipelines, e.g.:
    ```bash
    snakemake --snakefile fq2bedgraph.smk --config sample_sheet=config/sample_sheet.csv --cores 32 --printshellcmds --use-conda --keep-incomplete --keep-going --nolock --rerun-triggers mtime --benchmark-extended
    ```
8. Use the scripts in `benchmark` to reproduce the benchmarking results.
