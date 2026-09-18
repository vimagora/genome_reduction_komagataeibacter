# HGT analysis scripts — *Komagataeibacter rhaeticus* iGEM-hybrid

Analysis scripts for detecting horizontal gene transfer (HGT) in a 220 kb region
of the *K. rhaeticus* iGEM-hybrid genome.

Used for the analyses described in <TODO: paper reference / preprint URL when available>.
Contact the corresponding author for questions about the paper itself;
open an issue here for questions about the scripts.

## What's here

- `data_analysis.md` — the recipe: every step in order, the script that runs it,
  what it reads and what it writes.
- `config/` — the three inputs of the pipeline:
  - `genomes.tsv` — the genomes analysed (NCBI accessions and taxonomy)
  - `region.tsv` — the target genome and the coordinates of the region
  - `config.env` — paths and analysis settings
- `scripts/` — one folder per pipeline stage, numbered in running order.
- `envs/` — conda environment files.
- `deprecated/` — deprecated steps of the pipeline.

## Reproducing

The pipeline starts from the assembly of the *K. rhaeticus* iGEM-hybrid genome
and public NCBI records: it rotates and annotates the iGEM-hybrid assembly (dnaapler,
Bakta), downloads the reference genomes with their NCBI annotations, and every
other file is produced by the scripts. It assumes a SLURM cluster for the heavy
steps; the others run with `bash` or `python3`.

1. Put your paths, including the iGEM-hybrid assembly (`TARGET_ASSEMBLY`), in
   `config/local.env` (overrides `config/config.env`).
2. Create the conda environments from `envs/`.
3. Follow `data_analysis.md` from step 0.

| Tool | Version | Environment |
|---|---|---|
| dnaapler | 1.4.0 | `envs/dnaapler.yaml` |
| Bakta | 1.12.0 (database v6.0, full) | `envs/bakta.yaml` |
| NCBI Datasets CLI | not pinned | `envs/ncbi-tools.yaml` |
| BLAST+ | 2.17 | `envs/ncbi-tools.yaml` (or a cluster module, `BLAST_SETUP`) |
| TaxonKit | 0.9.0 | `envs/ncbi-tools.yaml` |
| BUSCO | 6.1.0 | `envs/busco.yaml` |
| OrthoFinder | 2.5.5 | `envs/orthofinder.yaml` |
| GTDB-Tk | 2.7.2 (reference data release 232) | `envs/gtdbtk.yaml` |
| ete3, pandas, matplotlib | not pinned (Python 3.10) | `envs/ete3.yaml` |

## Data availability

The reference genomes are public NCBI records listed in `config/genomes.tsv`. The
*K. rhaeticus* iGEM-hybrid genome and its sequencing reads are available under
<TODO: accession>; see also the Data Availability statement of the article.

## License

The code is released under the MIT License; see `LICENSE`.

## Deprecated analysis

`deprecated/discarded_analysis.md` describes an HGT inference with RANGER-DTL
that is not part of the article; it is kept for reference. The article's HGT
analysis used alternative methods, which are described in the article.
