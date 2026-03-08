# MethFlow  [![zread](https://img.shields.io/badge/Ask_Zread-_.svg?style=flat&color=00b0aa&labelColor=000000&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB3aWR0aD0iMTYiIGhlaWdodD0iMTYiIHZpZXdCb3g9IjAgMCAxNiAxNiIgZmlsbD0ibm9uZSIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KPHBhdGggZD0iTTQuOTYxNTYgMS42MDAxSDIuMjQxNTZDMS44ODgxIDEuNjAwMSAxLjYwMTU2IDEuODg2NjQgMS42MDE1NiAyLjI0MDFWNC45NjAxQzEuNjAxNTYgNS4zMTM1NiAxLjg4ODEgNS42MDAxIDIuMjQxNTYgNS42MDAxSDQuOTYxNTZDNS4zMTUwMiA1LjYwMDEgNS42MDE1NiA1LjMxMzU2IDUuNjAxNTYgNC45NjAxVjIuMjQwMUM1LjYwMTU2IDEuODg2NjQgNS4zMTUwMiAxLjYwMDEgNC45NjE1NiAxLjYwMDFaIiBmaWxsPSIjZmZmIi8%2BCjxwYXRoIGQ9Ik00Ljk2MTU2IDEwLjM5OTlIMi4yNDE1NkMxLjg4ODEgMTAuMzk5OSAxLjYwMTU2IDEwLjY4NjQgMS42MDE1NiAxMS4wMzk5VjEzLjc1OTlDMS42MDE1NiAxNC4xMTM0IDEuODg4MSAxNC4zOTk5IDIuMjQxNTYgMTQuMzk5OUg0Ljk2MTU2QzUuMzE1MDIgMTQuMzk5OSA1LjYwMTU2IDE0LjExMzQgNS42MDE1NiAxMy43NTk5VjExLjAzOTlDNS42MDE1NiAxMC42ODY0IDUuMzE1MDIgMTAuMzk5OSA0Ljk2MTU2IDEwLjM5OTlaIiBmaWxsPSIjZmZmIi8%2BCjxwYXRoIGQ9Ik0xMy43NTg0IDEuNjAwMUgxMS4wMzg0QzEwLjY4NSAxLjYwMDEgMTAuMzk4NCAxLjg4NjY0IDEwLjM5ODQgMi4yNDAxVjQuOTYwMUMxMC4zOTg0IDUuMzEzNTYgMTAuNjg1IDUuNjAwMSAxMS4wMzg0IDUuNjAwMUgxMy43NTg0QzE0LjExMTkgNS42MDAxIDE0LjM5ODQgNS4zMTM1NiAxNC4zOTg0IDQuOTYwMVYyLjI0MDFDMTQuMzk4NCAxLjg4NjY0IDE0LjExMTkgMS42MDAxIDEzLjc1ODQgMS42MDAxWiIgZmlsbD0iI2ZmZiIvPgo8cGF0aCBkPSJNNCAxMkwxMiA0TDQgMTJaIiBmaWxsPSIjZmZmIi8%2BCjxwYXRoIGQ9Ik00IDEyTDEyIDQiIHN0cm9rZT0iI2ZmZiIgc3Ryb2tlLXdpZHRoPSIxLjUiIHN0cm9rZS1saW5lY2FwPSJyb3VuZCIvPgo8L3N2Zz4K&logoColor=ffffff)](https://zread.ai/YuanfengZhang/MethFlow)
Snakemake pipelines to run NGS-based methylation analysis and benchmark. [check the DeepWiki explaination](https://deepwiki.com/YuanfengZhang/MethFlow)
The python script and .ipynb files to reproduce the benchmarking results, including evaluation and statistical analysis, are in the `benchmark` folder.

The scripts for data visualization are listed in `benchmark/figures` folder.

Following is a simplified schematic diagram involving the MethFlow and MethCali:
![Schematic diagram](MISC/overview.png)

Here is the way to calculate the RMSE and SpearmanR:
![RMSE](MISC/RMSE.png)
![SpearmanR](MISC/SpearmanR.png)

Please read the `benchmark/README.md` for more details.

## Requirements
While the snakemake and python scripts are compatible with most operating systems, there are many bioinformatic tools evaluated here only work on x64 Linux. Please use a x64 Linux server / container with at least 128GB RAM and 8 CPU cores to run the pipelines.

## Installation

### Step 1: Clone the Repository

Clone this repository with all submodules:

```bash
git clone https://github.com/YuanfengZhang/MethFlow --recurse-submodules
cd MethFlow
```

### Step 2: Build the Docker Image

Build the Docker image from the provided Dockerfile:

```bash
docker build -t MethFlow:latest .
```

### Step 3: Build Genome Index Files

Use the `utils/build_index_en.sh` script inside the Docker container to build indices for the alignment tools. You need to provide your own reference genome file:

```bash
docker run -it --rm \
    -v /path/to/your/reference:/ref \
    -v /path/to/output:/output \
    MethFlow:latest \
    /opt/MethFlow/utils/build_index_en.sh \
    -r /ref/genome.fa -o /output -t bwa-meme,bwa-meth,bismark-bowtie2 -@ 16
```

**Note:** The `-t` flag specifies which tools to build indices for. Supported tools include: `bwa-meme`, `bwa-meth`, `bwa-mem`, `astair`, `batmeth2`, `bismark-bowtie2`, `bismark-hisat2`, `biscuit`, `bowtie2`, `bsbolt`, `bsmapz`, `fame`, `gatk`, `gem3`, `hisat2`, `hisat-3n`, `strobealign`, `whisper`.

### Step 4: Download Required Resources

Download additional required resources using the provided script:

```bash
bash utils/download.sh
```

This will download dbSNP and GENCODE annotation files to the `resources/` directory.

### Step 5: Configure Runtime Parameters

Copy and modify the Docker runtime configuration file:

```bash
cp config/runtime_config_docker.yaml config/my_runtime_config.yaml
```

Edit `config/my_runtime_config.yaml` to set:
- `input_dir`: Path to your input FASTQ files
- Reference genome paths under the `ref:` section (update to match your genome indices)
- `snp_vcf` and `dbsnp_file`: Path to these two files downloaded by `utils/download.sh` if you want to use `gatk-calibrate`.
- Any other tool-specific parameters

### Step 6: Prepare Sample Sheet

Create a sample sheet CSV file according to the format described in `utils/sample_sheet_parser.py`. See `utils/conda_trigger.csv` for an example. See `utils/sample_sheet_parser.py` for detailed usage.

### Step 7: Run the Pipeline

Execute the Snakemake pipeline using Docker:

```bash
docker run -it --rm \
    -v /folder/of/config:/data
    -v /path/to/your/input:/input \
    -v /path/to/output:/output \
    MethFlow:latest \
    snakemake --snakefile fq2bedgraph.smk \
        --config sample_sheet=/data/sample_sheet.csv \
        --cores 32 --use-conda --printshellcmds
```

**Note:** Mount your data directory to `/data` in the container and ensure your runtime config file paths are set accordingly. The `/data` folder should contain `runtime_config.yaml` and `sample_sheet.csv`.
