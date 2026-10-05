<h1 align="center">BTC-CCS</h1>
<p align="center"><b>Consensus cancer subtypes for biliary tract cancer</b></p>

<p align="center">
  <a href="https://www.biorxiv.org/content/10.64898/2025.12.12.693962v3"><img alt="bioRxiv preprint" src="https://img.shields.io/badge/preprint-bioRxiv-b31b1b"></a>
  <a href="https://pancurx-oicr.github.io/btc-ccs/"><img alt="Analysis walkthrough" src="https://img.shields.io/badge/walkthrough-R%20Markdown-2c7fb8"></a>
  <img alt="R" src="https://img.shields.io/badge/R-%E2%89%A54.0-276DC3?logo=r&logoColor=white">
  <a href="https://ega-archive.org/studies/EGAS50000000972"><img alt="EGA: EGAS50000000972" src="https://img.shields.io/badge/EGA-EGAS50000000972-555"></a>
  <a href="LICENSE"><img alt="License: all rights reserved" src="https://img.shields.io/badge/license-all%20rights%20reserved-lightgrey"></a>
</p>

This repository accompanies the manuscript
**"Integrated whole-genome and transcriptome sequencing reveals divergent
evolutionary processes across biliary tract cancer subtypes"** (under review;
[preprint on bioRxiv](https://www.biorxiv.org/content/10.64898/2025.12.12.693962v3)).

It contains two things:

| | What it is | Start here |
|---|---|---|
| 🔬 **Single-sample classifiers** | Command-line tools that assign one tumour (or many) to a consensus subtype, from RNA expression (**eCCS**) or from copy number (**gCCS**) | [Classify your own samples](#classify-your-own-samples) |
| 📖 **Manuscript analyses** | The code and step-by-step walkthroughs behind every result in the paper | [Reproduce the manuscript](#reproduce-the-manuscript) |

---

## Classify your own samples

Each classifier is a single R script that can be run from any directory. The
trained models and reference files in `data/` are used by default. The
bundled example inputs let you try both in under a minute. Add `--help` to
either script to see all options.

### Install

With [conda](https://github.com/conda-forge/miniforge) (recommended):

```bash
git clone https://github.com/PanCuRx-OICR/btc-ccs.git
cd btc-ccs
./install.sh             # creates the 'btc-ccs' env and tests it
conda activate btc-ccs
```

`install.sh` builds the environment from [`environment.yml`](environment.yml),
then runs [`test_install.sh`](test_install.sh). The test classifies the
bundled examples and checks the results against the expected outputs. On
Apple Silicon Macs it builds an Intel environment that runs under
Rosetta 2, because Bioconductor's CNTools has no native conda build there.

<details>
<summary>Without conda</summary>

In R:

```r
install.packages(c("data.table", "optparse", "stringr", "dplyr", "tibble", "BiocManager"))
BiocManager::install("CNTools")
```

Then check the installation with `./test_install.sh`. eCCS alone only needs
`data.table`, `optparse` and `stringr`.
</details>

### eCCS: expression subtype (RNA)

eCCS is a top-scoring-pair classifier. It compares expression within 16
gene pairs, and each pair votes for **eCCS-A** or **eCCS-B**. Because it uses
only within-sample comparisons, it needs no normalisation against a cohort.

```bash
Rscript call_btc_eCCS.R -i data/tpm.txt -o eCCS.txt
```

**Input** (`-i`): a tab-separated expression table with one row per gene and
one column per sample. The first column holds HGNC gene symbols (any column
name). The classifier was trained on TPM. Because it relies on rankings, it
has also behaved well on other quantifications. If a gene is missing, the
script tries its known aliases. If a value is missing (`NA`), that gene pair
is skipped for that sample, with a warning.

```
gene_list   sample_1   sample_2
OR4F5       0.1        0.0
KCNN4       225.5      12.3
...
```

**Output** (`-o`):

| Column | Meaning |
|---|---|
| `sample_id` | Column name from the input |
| `rna_class` | `eCCS-A`, `eCCS-B`, or `NA` if fewer than 11 of 16 pairs agree |
| `confidence.polarized` | Fraction of pairs voting eCCS-A (1 = all A, 0 = all B) |

### gCCS: genomic subtype (copy number)

gCCS is a logistic-regression classifier. It uses median arm-level copy
number on 1q, 4p, 13q, 18q and 19p, plus the sample's mean variant allele
frequency.

```bash
Rscript call_btc_gCCS.R -s data/example.seg -v data/mean_vaf.txt -o gCCS.txt
```

For hg19 segments, add `-b data/genebed.hg19.txt`.

**Inputs:**

| Flag | File | Format |
|---|---|---|
| `-s` | Copy-number segments | Tab-separated, columns `ID chrom loc.start loc.end num.mark seg.mean` (see [`data/example.seg`](data/example.seg)). Chromosomes may be named `chr1` or `1`. |
| `-v` | Mean VAF | Tab-separated, columns `sample mean.vaf`; `sample` must match the segment file's `ID` (see [`data/mean_vaf.txt`](data/mean_vaf.txt)) |
| `-b` | Gene coordinates | Default [`data/genebed.hg38.txt`](data/genebed.hg38.txt); must match the genome build of your segments |
| `-c` | Gene → chromosome arm map | Default [`data/cytoband.txt`](data/cytoband.txt) |
| `-m` | Trained model | Default [`data/cnv.model.rds`](data/cnv.model.rds) |

**Output** (`-o`):

| Column | Meaning |
|---|---|
| `sample` | Sample ID |
| `glm_labels` | `gCCS-A` if `prob` ≥ 0.5, otherwise `gCCS-B` |
| `prob` | Predicted probability of gCCS-A |

Example outputs for the bundled inputs are in
[`data/eCCS.txt`](data/eCCS.txt) and [`data/gCCS.txt`](data/gCCS.txt).

---

## Reproduce the manuscript

### Browse the walkthrough

Every analysis in the paper is a rendered R Markdown page, linked in order on
the **[project site](https://pancurx-oicr.github.io/btc-ccs/)**:

| Part | Topic |
|---|---|
| 0 | Cohort setup: tumour specimens and case histories |
| 1 | Consensus subtypes: signature network, classifier training and robustness, anatomical location, survival (internal and external cohorts) |
| 2 | Expression: differential expression, NMF, immune signatures, and re-analysis of three public single-cell datasets |
| 3 | Genetic drivers: small mutations, structural-variant and copy-number recurrence, co-occurrence |
| 4 | Mutational landscape: MSI / HRD / TMB, subtype differences, and evolutionary timing |

The `.Rmd` sources and their rendered `.html` pages are in [`docs/`](docs).

### Start from raw data

Raw FASTQ files are on the European Genome-phenome Archive under
[EGAS50000000972](https://ega-archive.org/studies/EGAS50000000972). Access
requires a Data Access Agreement. Variants were called with the pipeline
described in
[Chan-Seng-Yue *et al.* 2020](https://www.nature.com/articles/s41588-019-0566-9).
The HPC analyses are launched by
[`bin/0_evolution.collations.sh`](bin/0_evolution.collations.sh), which calls
the other shell scripts in `bin/`.

## Citation

If you use the classifiers or code, please cite the manuscript. Until it is
published, cite the
[bioRxiv preprint](https://www.biorxiv.org/content/10.64898/2025.12.12.693962v3).

## License

Copyright © 2026 Felix Beaudry. All rights reserved.

This is unpublished, proprietary work, provided **for viewing only**. No
part of it may be copied, modified, distributed or used without prior
written permission from the author. See [LICENSE](LICENSE).
