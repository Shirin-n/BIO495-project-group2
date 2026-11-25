# Coral NLR and microbiome analysis

This project tests whether coral species with more predicted NOD-like receptor (NLR) genes have more diverse tissue microbiomes, and whether species with similar candidate counts have similar microbial communities. The complete workflow is in **`coral_nlr_microbiome_analysis.ipynb`**, with English explanations, executable code, tables and figures.

The primary analysis includes **8 species and 25 tissue samples**. Neither primary test supports the predicted positive association. This is not evidence that immune genes have no biological effect or that host species have identical microbiomes.

## Start here

This repository includes the executed notebook, curated results and figures, versioned reference snapshots, and the BUSCO OCR record. The saved results are ready to inspect. The download and notebook steps below recreate the source dataset and derived files needed for a full run.

To rerun, first create the `bio495` environment using `requirement.txt` as described below. Then run these commands from the project root:

```bash
conda activate bio495
python scripts/download_project_data.py
python scripts/verify_downloaded_data.py
jupyter lab coral_nlr_microbiome_analysis.ipynb
```

The download step retrieves 255 files (about 786 MiB) from the fixed source commits, verifies sizes and Git blob checksums, and exits unsuccessfully if any file fails validation. It uses the system `curl` command. It requires space for missing downloads plus a 3 GiB working reserve; reuse of existing downloads does not count their size as new space. The small reference files under `data/reference/` are already included, so no live reference API fetch is needed.

Select **Python (bio495)** and choose **Restart Kernel and Run All Cells**. The kernel working directory must be the project root. Run cells in order: they rebuild processed datasets and screening intermediates before the later analyses consume them. Outputs and summary reports are regenerated. A rerun may replace files in `results/`; the distributed results represent a verified run.

A clean download and complete notebook run were verified in the `bio495` environment on Apple Silicon using Rosetta. The two primary test results matched the earlier run exactly.

## Project structure

```text
BIO495-project-group2/
├── README.md
├── requirement.txt
├── coral_nlr_microbiome_analysis.ipynb
├── config/                           Fixed settings, mappings and source metadata
├── data/
│   ├── reference/                    Exact HMM and reference protein snapshots
│   └── interim/proteins/busco_ocr.json  Required checked OCR record
├── scripts/                          Download, analysis, summary and test code
├── results/                          Selected results, figures and validation records
└── docs/                             Fixed protocol and decision record
```

The notebook calls the reusable scripts. Share this directory as a unit. The notebook alone does not contain the supporting scripts, reference snapshots or original source datasets; the download step restores the omitted source datasets.

## Folder roles and runtime dependencies

“Required” below means required during an end-to-end run, including verification. The raw, processed and report folders are absent from this shared copy and are created again during download or execution. Some folders contain both analytical inputs and audit records; they cannot safely be classified as disposable checking material as a whole.

| Folder | Contents and purpose | Relationship to running the notebook |
| --- | --- | --- |
| `config/` | Frozen analysis settings, species and reference manifests, sample mappings, source versions and checksums | **Required input.** The notebook and scripts read settings and manifests, regenerate selected preparation manifests, and verify frozen representations. Download receipts are provenance records, but the folder as a whole is essential. |
| `data/raw/` | Original protein/GTF annotations, microbiome BIOM/metadata and source host tree | **Downloaded input.** Omitted from the shared copy; restored by the download script and then read unchanged. |
| `data/reference/` | Pfam profiles, reference protein source records and FASTA files | **Required input.** Used for NLR identification and controls; some derived reference files are regenerated. |
| `data/interim/` | BUSCO OCR transcription, temporary protein queries and HMMER/BLAST intermediates | **Required OCR record plus generated cache.** Only `proteins/busco_ocr.json` is included. It is read for genome QC and is not recreated by the notebook. NLR intermediate folders are regenerated during execution. |
| `data/processed/` | Cleaned metadata, cohort count matrices, NLR exposures and diversity metrics | **Generated outputs reused as inputs.** Omitted from the shared copy; earlier cells produce the files later cells read. Cohort representations are also checked against the frozen manifest. Run in order if regenerating them. |
| `scripts/` | NLR screening, diversity calculations, statistics, figures and test code | **Required executable code.** The notebook calls the analysis, summary and test scripts. Download/verification helpers and `ocr_busco.swift` support initial preparation and are not called in a normal notebook run. `__pycache__/` is a disposable Python cache. |
| `results/` | QC tables, candidate evidence, summary metrics, figures, permutation distributions and validation records | **Generated outputs and verification inputs.** Later cells display and validate earlier outputs. Existing screening outputs may be reused after hash checks. It is not merely an optional inspection folder. |
| `docs/` | Fixed protocol and decision record | **Required verification input.** `analysis_protocol_v1.txt` and `analysis_decisions_v1.tsv` are checksum-checked. The literature PDF is omitted; source citations remain in this README and the notebook. |
| `reports/` | Four existing technical summaries from the NLR, microbiome and association stages | **Generated summaries, also checksum-verified.** Omitted from the shared copy; summary scripts recreate these files; they are not source data and are not a separate final project report. The directory must be writable and is created by setup if missing. |

Within `results/`, `qc/` stores input/eligibility audits, `tables/` stores coverage summaries, and `figures/` contains the notebook figures. `nlr_pilot_v1/` and `nlr_full_v1/` retain candidate summaries and selected screening evidence; `microbiome_metrics_v1/` retains depth, ordination and eligibility diagnostics; `association_v1/` retains tests, permutation nulls and the independent audit. `statistics/` is currently an empty reserved directory and is not used by the notebook.

The repository includes result summaries, figures and validation records. A full notebook run recreates the detailed NLR and feature taxonomy tables, local execution logs, and other generated files listed in the validation records.

## Notebook contents

| Section | Content |
| --- | --- |
| 1 Environment and frozen plan | Environment checks, metadata integrity and protocol settings |
| 2 Metadata cleaning | Missing values, species reconciliation, tissue and health eligibility, tree mapping |
| 3 Genome quality | Annotation statistics, source BUSCO summaries and reference resources |
| 4 Microbiome preparation | BIOM validation, taxonomy filtering, mitochondrial QC, depth thresholds and cohort exports |
| 5 NLR identification | Three-species pilot, controls, nine-species screening and candidate review |
| 6 Diversity and composition | Standardized alpha diversity, Bray–Curtis distances, PCoA and depth diagnostics |
| 7 Association analysis | Primary Spearman and Mantel tests, all planned sensitivities and species influence |
| 8 Interpretation and verification | Main conclusion, limitations, independent audits and compact results |

The notebook contains data preparation, cleaning, preprocessing and exploratory figures in one place. It reuses the published OTU table; it does not reprocess raw 16S reads or run QIIME 2. Stage-specific flags such as `association_tests_run: false` describe a preparation stage, not the completion status of the whole workflow.

## Environment and dependencies

**`requirement.txt` is the sole maintained environment specification.** It contains 128 pinned Python distributions and an embedded comment block with 264 exact conda package builds, including Python **3.11.18**, HMMER **3.4**, BLAST+ **2.17.0** and native libraries. Comments are ignored by pip, so the Python list remains a standard requirements file.

### Recreate the complete environment

From the project root, use bash or zsh and an existing conda installation. Run this only if an environment named `bio495` does not already exist:

```bash
conda create --name bio495 --platform osx-64 --override-channels -c conda-forge -c bioconda --file <(sed -n 's/^# conda-lock: //p' requirement.txt)
conda activate bio495
python -m ipykernel install --user --name bio495 --display-name "Python (bio495)"
python -m pip check
```

The shell exposes the embedded specification to conda without creating another maintained environment file. This includes all Python packages and sequence tools; do not run a second pip installation over the recreated environment.

The exact specification targets **macOS osx-64**. On Apple Silicon it runs under Rosetta. Recreating the environment requires access to the package archives in the lock. The verified build matched all Conda and Python pins and passed `pip check`.

### Python packages only

For a compatible environment where Python, sequence tools and native dependencies are already prepared:

```bash
python -m pip install -r requirement.txt
```

This command installs the Python layer only. Pip does not install HMMER, BLAST or Python itself; use the complete recreation command above for the full project environment. All package versions and commands are recorded in the same `requirement.txt`.

## Rerunning and validation

All notebook data, configuration, script and output paths are relative to the project root through `ROOT = Path(".")`. Start Jupyter from the repository root and use that folder as the kernel working directory. There is no hard-coded user directory or environment-variable override. When moving the project, move the entire folder and start Jupyter from its new root.

The Python executable and installed tool locations come dynamically from the active conda environment; these are system executable locations rather than project data paths. Supporting scripts derive the project root from their own location, so they do not depend on the original Desktop path. Local execution logs record the actual command paths used during each run.

For a headless notebook run from the project directory:

```bash
jupyter nbconvert --to notebook --execute coral_nlr_microbiome_analysis.ipynb \
  --ExecutePreprocessor.kernel_name=bio495 \
  --ExecutePreprocessor.timeout=1800 \
  --output=coral_nlr_microbiome_analysis.ipynb --output-dir=.
```

Optionally generate a local HTML reading copy after execution (excluded from the shared source files):

```bash
jupyter nbconvert --to html coral_nlr_microbiome_analysis.ipynb --output-dir=.
```

The notebook runs the stage checks automatically. Individual checks can also be run from the project root after the notebook has regenerated their required inputs:

```bash
python scripts/test_nlr_pilot.py
python scripts/test_microbiome_metrics.py
python scripts/test_association_tests.py
python scripts/audit_association_results.py
```

Key validation records are `results/nlr_full_v1/validation.json`, `results/microbiome_metrics_v1/validation.json`, `results/association_v1/independent_audit.json`. HMMER cache reuse makes a repeat run faster than a fresh screening run; elapsed time depends on the computer and cache state.

## Data provenance and fixed methods

- **Genome annotations:** [Metz coral annotation release v1.0.0](https://doi.org/10.5281/zenodo.17713987), published 25 November 2025; repository commit `5d18908cfb89af83181d5205f003c05d9494c0cc`. Forty annotation sets were collected, nine overlapping species were screened, and eight entered the primary association tests.
- **Microbiome:** [Pollock et al. 2018](https://doi.org/10.1038/s41467-018-07275-x) and the [GCMP Australia analysis repository](https://github.com/zaneveld/GCMP_Australia_Coevolution), commit `fba13cd74d717727d1990e5848d9c86695454a90`. The original BIOM has 129,305 OTUs and 692 samples; the metadata has 695 records.
- **Primary cohort:** coral tissue, strict recorded health criteria, source mitotype agreement, assigned Bacteria/Archaea and at least 1,000 filtered reads. Recorded colony identifiers do not establish biological independence.
- **NLR exposure:** unique GTF genes whose longest eligible protein contains accepted NACHT and downstream LRR domains. Model gathering thresholds and at least 50% HMM coverage are required. BLAST is supporting evidence, not experimental validation.
- **Diversity:** 100 subsamples without replacement at 1,000 reads, natural-log Shannon, with recorded-colony and species aggregation. Composition uses equally aggregated unrarefied relative abundances and Bray–Curtis distances.
- **Inference:** species-level, two-sided permutations. Both primary tests enumerate 40,320 label permutations and receive Holm correction as a two-test family. Sensitivity p values are exploratory and unadjusted.

The protocol is `docs/analysis_protocol_v1.txt`, with machine-readable settings in `config/analysis_protocol_v1.json`. The 65 frozen input representations are listed in `config/analysis_inputs_v1.tsv`. The plan was fixed locally after preparation and before identification and association testing; it was not externally preregistered. Method changes should be recorded as a versioned amendment.

## Main results and where to find them

| Primary test | Species | rho | Raw p | Holm p |
| --- | ---: | ---: | ---: | ---: |
| A Candidate count versus Shannon | 8 | -0.227549 | 0.585913 | 1.000000 |
| B Gene-count distance versus Bray–Curtis | 8 | -0.088469 | 0.617609 | 1.000000 |

Nine screened species contain **266 canonical architecture candidates**, including **238** in the eight primary species. Seventeen analysis scenarios produce two primary tests, 30 eligible sensitivity endpoint tests and two descriptive-only endpoint rows. None of the 30 eligible sensitivity tests has raw p < 0.05; these tests share data and are not independent replications.

| Output | Location |
| --- | --- |
| Cleaned primary counts and metadata | `data/processed/primary_strict_mito_agree/` |
| NLR counts and candidate sequences | `data/processed/nlr_v1/` |
| Species-level diversity and distance inputs | `data/processed/microbiome_metrics_v1/` |
| Primary test results | `results/association_v1/primary_results.tsv` |
| Every scheduled endpoint | `results/association_v1/all_results.tsv` |
| Species omission diagnostics | `results/association_v1/leave_one_species_out.tsv` |
| Primary association figure | `results/figures/association_primary_v1.png` |
| Sensitivity figure | `results/figures/association_sensitivity_v1.png` |
| Independent statistical audit | `results/association_v1/independent_audit.json` |

## Interpretation limits

Only eight species contribute to primary inference. Sampling is uneven, five species have fewer than three recorded colonies, and 1,000-read standardization does not establish diversity saturation. Predicted NLR architecture counts are not experimentally validated immune repertoires; 22 candidates have overlapping core-model assignments.

A suitable species-resolved host phylogeny is unavailable, so **PGLS was not performed**. Ancestry and geography were not adjusted; species exchangeability is an assumption of the permutation tests. The single-expedition analysis is a sensitivity check. No PERMANOVA or partial Mantel test was added. Failure to support A or B does not establish the original separate hypothesis of no host-species effect.

## Storage and handoff

The tracked files include the code, executed notebook, core outputs, fixed analysis settings, small reference snapshots, and BUSCO OCR record. `requirement.txt` is the environment specification.

See `.gitignore` for local-only files, including downloaded datasets, derived matrices, working caches, detailed evidence tables, machine-specific result logs, and course materials. The rules explicitly retain `data/interim/proteins/busco_ocr.json`.

The three auxiliary source-method files previously downloaded outside the main 255-file manifest are not needed by the runtime. Their source URLs remain recorded in source metadata or the pinned source repository. The fixed reference snapshots avoid relying on changing live model or protein records when reproducing v1.
