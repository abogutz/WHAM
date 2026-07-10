# WHAM

**W**GBS **H**eterogeneity **A**nalysis at the **M**olecular level

A toolkit for analysing the distribution of DNA methylation across individual
reads from Bismark-aligned WGBS data. It produces UCSC genome-browser tracks for
three analyses:

1. **Lollipop** visualization of per-read methylation calls
2. **Heatmap** of read-level methylation distribution
3. **Diptest** modality test for non-unimodal (heterogeneous) methylation

## Installation

Install all dependencies into an isolated environment using conda:

```bash
git clone <this-repo-url> WHAM
cd WHAM
conda env create -f environment.yml
conda activate wham
```

> **Apple Silicon (osx-arm64):** some UCSC tools (`ucsc-*`) are only built for
> `osx-64` on bioconda. If conda cannot find them, create the env under the
> `osx-64` subdir (runs via Rosetta 2) using this trick:
> ```bash
> CONDA_SUBDIR=osx-64 conda env create -f environment.yml
> ```

### Dependencies

Installed automatically by `environment.yml`:

- `samtools`, `bedtools`, `gawk`
- UCSC `bedToBigBed`, `bedGraphToBigWig`
- `R` with `diptest`, `foreach`, `doParallel`
- Optional (only for `-p` traditional tracks): `bismark`, `deeptools` (`bamCoverage`)

Verify everything is available:

```bash
./WHAM.sh -d
```

## Usage

```bash
./WHAM.sh -i <bismark.bam> -z <genome.sizes> -o <output_dir>
```

The input must be a **Bismark-aligned** BAM (reads carry `XM:Z:` methylation and
`XG:Z:` genome-strand tags). Both single-end and paired-end data are supported;
paired mates are automatically combined. Output is written to
`<output_dir>/Track_Hub/`, ready to load at in the UCSC Genome Browser (Track Hubs).

### Example

```bash
./WHAM.sh -i Rat_BNxWKY_Epiblast_rep1.bam -z rn7.sizes -o results
```

### Options

| Flag | Description | Default |
|------|-------------|---------|
| `-i` | Input Bismark `.bam` file (**required**) | |
| `-z` | Chromosome sizes file | mm10 |
| `-o` | Output directory for `Track_Hub` | current directory |
| `-s` | Scratch directory | a temp folder created beside the input BAM |
| `-t` | Number of threads | SLURM allocation, else all cores |
| `-q` | Minimum mapping quality | 40 |
| `-C` | Minimum CpGs per read | 4 |
| `-G` | Genome bin size for Diptest | 100 |
| `-H` | Genome bin size for Heatmap | 25 |
| `-B` | Number of methylation bins for Heatmap | 5 |
| `-R` | Max reads for color scaling in Heatmap | 10 |
| `-c` | Number of color bins for Heatmap | 5 |
| `-p` | Also create traditional methylation/coverage tracks | off |
| `-l` | Skip lollipop track | on |
| `-m` | Skip heatmap track | on |
| `-D` | Skip diptest track | on |
| `-b` | Reference BED file (runs ONLY diptest over these regions) | |
| `-d` | Check dependencies and exit | |
| `-h` | Full help | |