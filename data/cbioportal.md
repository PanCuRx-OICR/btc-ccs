# External datasets (`source_data/cbioportal/`)

The manuscript walkthroughs in [`docs/`](../docs) validate the classifiers and
findings on published biliary tract cancer cohorts. These data are not
included in this repository. To re-run those steps, download them into
`source_data/cbioportal/` inside your data directory (`data.dir` in each
`.Rmd`), arranged as below.

## Expected layout

```
source_data/cbioportal/
├── gene_info.txt
├── chol_tcga_gdc/
│   ├── data_cna_hg38.seg
│   ├── data_mutations.txt
│   └── data_mrna_seq_tpm.txt
├── chol_icgc_2017/
│   ├── data_clinical_sample.txt
│   ├── data_clinical_patient.txt
│   ├── data_mutations.txt
│   ├── GSE89747_CCA_batch01_illumina_Gene_expression_noNorm_noBKGD.txt
│   ├── GSE89748_CCA_batch02_illumina_Gene_expression_noNorm_noBKGD.txt
│   ├── 21598290cd170368-sup-181437_2_supp_4118307_1ryl1l.xlsx
│   ├── 21598290cd170368-sup-181437_2_supp_4118310_krylkl.xlsx
│   └── GSE89803_RAW/
│       ├── jusakul_methylation.txt
│       ├── GPL16304_Gene_features_PlatformTable.txt
│       └── *.idat   (methylation array files)
├── hcc_msk_2024/
│   ├── data_clinical_patient.txt
│   ├── data_clinical_sample.txt
│   ├── data_mutations.txt
│   └── data_cna_hg19.seg
├── ihch_msk_2021/
│   ├── data_clinical_patient.txt
│   ├── data_mutations.txt
│   └── data_cna_hg19.seg
└── ihch_mskcc_2020/
    └── data_clinical_patient.txt
```

## Studies

| Folder | Paper | cBioPortal study | Used in |
|---|---|---|---|
| `chol_tcga_gdc/` | [Farshidfar et al. 2017](http://dx.doi.org/10.1016/j.celrep.2017.02.033) | [chol_tcga_gdc](https://www.cbioportal.org/study/summary?id=chol_tcga_gdc) | [1.2.2](../docs/1.2.2_CNV_classifier.Rmd) |
| `chol_icgc_2017/` | [Jusakul et al. 2017](https://doi.org/10.1158/2159-8290.CD-17-0368) | [chol_icgc_2017](https://www.cbioportal.org/study/summary?id=chol_icgc_2017) | [1.3.1](../docs/1.3.1_external_survival.Rmd), [2.0.2](../docs/2.0.2.external_expression.Rmd) |
| `hcc_msk_2024/` | [Song Y et al. 2024](https://pubmed.ncbi.nlm.nih.gov/38864854/) | [hcc_msk_2024](https://www.cbioportal.org/study/summary?id=hcc_msk_2024) | [1.3.1](../docs/1.3.1_external_survival.Rmd) |
| `ihch_msk_2021/` | [Boerner et al. 2021](https://pubmed.ncbi.nlm.nih.gov/33765338/) | [ihch_msk_2021](https://www.cbioportal.org/study/summary?id=ihch_msk_2021) | [1.3.1](../docs/1.3.1_external_survival.Rmd) |
| `ihch_mskcc_2020/` | [Jolissaint et al. 2020](https://pubmed.ncbi.nlm.nih.gov/33963001/) | [ihch_mskcc_2020](https://www.cbioportal.org/study/summary?id=ihch_mskcc_2020) | [1.3.1](../docs/1.3.1_external_survival.Rmd) |

## Files

The `data_*` files come from each study's cBioPortal download (the study
links above). The other files come from elsewhere:

| File | What it is | Source | Used in |
|---|---|---|---|
| `gene_info.txt` | Tab-separated gene ID lookup with two columns, `ENTREZID` and `Hugo_Symbol`, used to match TCGA expression to gene symbols | | 1.2.2 |
| `chol_icgc_2017/GSE89747_…_noBKGD.txt` | Jusakul expression, batch 1 | [GEO GSE89747](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89747) | 1.3.1 |
| `chol_icgc_2017/GSE89748_…_noBKGD.txt` | Jusakul expression, batch 2 | [GEO GSE89748](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89748) | 1.3.1 |
| `chol_icgc_2017/21598290cd170368-…_4118307_1ryl1l.xlsx` | Jusakul et al. 2017 supplementary table (read as S3A) | | 2.0.2 |
| `chol_icgc_2017/21598290cd170368-…_4118310_krylkl.xlsx` | Jusakul et al. 2017 supplementary table (read as S1A) | | 2.0.2 |
| `chol_icgc_2017/GSE89803_RAW/*.idat` | Jusakul methylation arrays (raw); must be unzipped from `.idat.gz` | [GEO GSE89803](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89803) | 2.0.2 |
| `chol_icgc_2017/GSE89803_RAW/GPL16304_Gene_features_PlatformTable.txt` | Methylation array probe annotation | [GEO GPL16304](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16304) | 2.0.2 |
| `chol_icgc_2017/GSE89803_RAW/jusakul_methylation.txt` | Metadata for the 142 methylation arrays (tab-separated, no header). Columns 2–4 (GEO sample ID, chip ID, array position) are joined with `_` to match the `.idat` names. | Provided in this repo as [`jusakul_methylation.txt`](jusakul_methylation.txt); copy it into `GSE89803_RAW/` | 2.0.2 |
