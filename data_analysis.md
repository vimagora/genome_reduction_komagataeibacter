# Data analysis

Step-by-step recipe of the analysis: what to run, in which order, what each
step reads and what it writes. Results and their interpretation are in the
article.

Every command is run from the repository root. All outputs are written under
`$PROJECT_ROOT`; paths below are relative to it.

## 0. Setup

**Inputs.** The pipeline starts from three files in `config/`:

| File | Content |
|---|---|
| `genomes.tsv` | One row per genome: sample name, NCBI assembly accession, taxonomy (`Phylum`, `Class`, `Order`, `Family`), BUSCO lineage, and `include` (`yes` for genomes used from step 5 on) |
| `region.tsv` | Target genome (sample name), chromosome accession, start and end of the region |
| `config.env` | Paths to the work directory, environments and databases, plus analysis settings |

Fill in the `<TODO>` accessions in `genomes.tsv` and `region.tsv`. Put your paths
and any other changed settings in `config/local.env` (not tracked by git); it is
loaded at the end of `config.env` and overrides its values.

**Software.** Create the conda environments:

```bash
set -a; source config/config.env; set +a
for env in ncbi-tools busco orthofinder gtdbtk ete3; do
    conda env create -f envs/$env.yaml -p "$ENV_ROOT/$env-env"
done
```

**Databases.** Set in `config/local.env`: GTDB-Tk reference data release 232
(`GTDBTK_DATA_PATH`), the NCBI taxonomy dump (`TAXONKIT_DB`), the NCBI `nr`
BLAST database with its taxonomy files (`BLASTP_DB`), and the command that makes
`blastp` available (`BLAST_SETUP`; on CSC Roihu, the `blast-plus` module). BUSCO
downloads its lineage datasets into `BUSCO_DOWNLOADS` on first use.

The NCBI taxonomy dump is downloaded into `TAXONKIT_DB` (where internet is
available); only the four `.dmp` files read by taxonkit are extracted:

```bash
mkdir -p "$TAXONKIT_DB"
wget -P "$TAXONKIT_DB" https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz \
                       https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz.md5
(cd "$TAXONKIT_DB" && md5sum -c taxdump.tar.gz.md5 \
    && tar -xzf taxdump.tar.gz names.dmp nodes.dmp delnodes.dmp merged.dmp)
```

The taxonomy is updated frequently; record the download date.

**Each session.** Load the configuration and create the SLURM log folders:

```bash
set -a; source config/config.env; set +a
PY="$ENV_ROOT/ete3-env/bin/python3"
mkdir -p logs/busco logs/blast logs/orthofinder logs/gtdbtk
```

SLURM scripts are submitted with `sbatch --account="$SLURM_ACCOUNT"`.

## 1. Genomes

Run `download_genomes.sh` to download the genome sequence, proteins and GFF3
annotation of every genome in `genomes.tsv` from NCBI (run it where internet is
available):

```bash
bash scripts/01_genomes/download_genomes.sh
```

Output: `genomes/fna/<sample>.fna`, `genomes/faa/<sample>.faa`,
`genomes/gff/<sample>.gff`. They are read by steps 2, 3, 5 and 7.

## 2. Genome quality

Run `busco_array.sh` to compute BUSCO completeness of each genome in `genomes.tsv`
(one array task per row), with the lineage from the `busco_lineage` column:

```bash
sbatch --account="$SLURM_ACCOUNT" --array=1-$(($(wc -l < "$GENOMES") - 1))%8 scripts/02_genome_qc/busco_array.sh
```

Output: `busco/<sample>/`. It is read by `generate_qc_report.py`, which combines
BUSCO scores, assembly statistics from `genomes/fna/` and the taxonomy in
`genomes.tsv` into one HTML report:

```bash
"$PY" scripts/02_genome_qc/generate_qc_report.py \
    --busco-dir "$PROJECT_ROOT/busco" --assembly-dir "$PROJECT_ROOT/genomes/fna" \
    --metadata "$GENOMES" --output "$PROJECT_ROOT/busco/genome_qc_report.html"
```

Output: `busco/genome_qc_report.html`. The `include` column of `genomes.tsv`
records which genomes are used from step 5 on.

## 3. Genes of the region

Run `extract_region_genes.py` to list the protein-coding genes of the region from
the target's NCBI annotation. A CDS belongs to the region when its midpoint lies
between the coordinates in `region.tsv`; transposases/IS (products matching
`IS_PATTERN`) and pseudogenes are excluded.

```bash
"$PY" scripts/03_region/extract_region_genes.py
```

Output: `region/region_cds.tsv` (every CDS in the region, with exclusion flags),
`region/genes_of_interest.txt` (protein IDs, genome order) and
`region/genes_of_interest.faa`. The protein list is read by steps 4 and 6;
the sequences by step 4.

## 4. Outgroup search

Run `split_query.sh` to split the genes of interest into `BLAST_CHUNKS` query files:

```bash
bash scripts/04_outgroup_search/split_query.sh
```

Output: `blast_results/query_chunks/chunk_<NNN>.faa`. They are read by
`blastp_nr.sh`, which searches each chunk (one array task per chunk) against
`BLASTP_DB`, restricted to the taxa in `BLASTP_TAXIDS`:

```bash
N=$(ls "$PROJECT_ROOT"/blast_results/query_chunks/chunk_*.faa | wc -l)
sbatch --account="$SLURM_ACCOUNT" --array=1-${N} scripts/04_outgroup_search/blastp_nr.sh
```

Output: `blast_results/chunk_results/chunk_<NNN>.tsv`. They are read by
`merge_blastp.sh`, which checks that every chunk finished and concatenates them:

```bash
bash scripts/04_outgroup_search/merge_blastp.sh
```

Output: `blast_results/blastp_nr.tsv`. It is read by `taxonomy_lookup.sh`, which
retrieves the lineage of every hit taxon with taxonkit:

```bash
bash scripts/04_outgroup_search/taxonomy_lookup.sh
```

Output: `blast_results/taxid_lineage_lookup.tsv` and `taxonkit_warnings.log`
(taxids of the BLAST database that were merged or deleted in the current
taxonomy). Both tables are read by `add_taxonomy.py`, which joins them and drops
hits whose taxid has no lineage:

```bash
"$PY" scripts/04_outgroup_search/add_taxonomy.py
```

Output: `blast_results/results_with_taxonomy.tsv`. It is read by
`taxonomy_coverage_check.py`, which counts the genes with hits outside the
target's family, order, class and phylum, and by
`outgroup_candidate_selection.py`, which ranks the species outside each of the
target's taxa by the number of genes they match (top `OUTGROUP_TOP_N`):

```bash
"$PY" scripts/04_outgroup_search/taxonomy_coverage_check.py
"$PY" scripts/04_outgroup_search/outgroup_candidate_selection.py
```

Both print their tables. The genomes chosen from these candidates are those
listed in `genomes.tsv`.

## 5. Phylogenomics

### 5.1 Inputs

Run `prepare_inputs.sh` to link the proteomes and genomes of the rows with
`include` = `yes`:

```bash
bash scripts/05_phylogenomics/prepare_inputs.sh
```

Output: `orthofinder/proteins/` (read by 5.2) and `gtdbtk/input/` (read by 5.3 and 5.5).

### 5.2 Orthogroups

Run `orthofinder.sh` to infer orthogroups and gene trees (DIAMOND search, MSA gene trees):

```bash
sbatch --account="$SLURM_ACCOUNT" scripts/05_phylogenomics/orthofinder.sh
```

Output: `orthofinder/run_<timestamp>/Results_<date>/`, linked as
`orthofinder/results`. It is read by step 6.

### 5.3 GTDB taxonomy

Run `gtdbtk_classify.sh` to classify the genomes with GTDB-Tk `classify_wf`:

```bash
sbatch --account="$SLURM_ACCOUNT" scripts/05_phylogenomics/gtdbtk_classify.sh
```

Output: `gtdbtk/results/classify/gtdbtk.bac120.summary.tsv`. It is read by 5.4.

### 5.4 Custom taxonomy

Run `make_custom_taxonomy.sh` to build the taxonomy file for the de novo tree;
genomes without a GTDB species get `s__<genus> sp`:

```bash
bash scripts/05_phylogenomics/make_custom_taxonomy.sh
```

Output: `gtdbtk/custom_taxonomy.tsv`. It is read by 5.5.

### 5.5 Species tree

Run `gtdbtk_de_novo.sh` to build the species tree with GTDB-Tk `de_novo_wf`,
without GTDB reference genomes, rooted on `GTDBTK_OUTGROUP`:

```bash
sbatch --account="$SLURM_ACCOUNT" scripts/05_phylogenomics/gtdbtk_de_novo.sh
```

Output: `gtdbtk/results/de_novo/gtdbtk.bac120.decorated.tree`. It is read by 5.6.

### 5.6 Tree conversion

Run `species_tree.sh` to write an iTOL version and a plain newick version of the tree:

```bash
bash scripts/05_phylogenomics/species_tree.sh
```

Output: `hgt_analysis/gtdbtk_itol.tree` and `hgt_analysis/species_tree.newick`.
The newick tree is used for HGT analysis as described in the article.

## 6. Horizontal gene transfer

### 6.1 Orthogroups of the region

Run `extract_hgt_orthogroups.sh` to find the orthogroup of each gene of interest
(searching the target's column of `Orthogroups.tsv`) and collect their gene trees
and gene counts:

```bash
bash scripts/06_hgt/extract_hgt_orthogroups.sh
```

Output: `hgt_analysis/gene_to_orthogroup.tsv`, `orthogroups_of_interest.txt`,
`gene_trees/` and `gene_count_matrix.tsv`.

The outputs are used for HGT analysis as described in the article.

## 7. Transposase enrichment

Run `permutation_test_transposases.py` to compare the number of transposases/IS
in the region with `N_PERMUTATIONS` random windows of the same size on the same
sequence (not overlapping the region; seed `PERMUTATION_SEED`), using the
target's genome and GFF3 annotation from step 1:

```bash
"$PY" scripts/07_statistics/permutation_test_transposases.py
```

Output: printed test statistics and `statistics/permutation_null_distribution.txt`.
This step only depends on step 1.
