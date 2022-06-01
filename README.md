This repository implements an example workflow for performing metapangenomics using amino acid k-mers, and compares this pipeline against existing metapangenome approaches.
The results of this workflow are written up at github.com/dib-lab/2021-paper-metapangenomes.

## Getting started

The workflow is written in snakemake with software managed by conda. The workflow can be re-run with the following (snakemake command assumes a slurm cluster, and executes as a dry run first).

```
conda env create --name orph --file environment.yml
conda activate orph

snakemake -j 16 --use-conda --rerun-incomplete --latency-wait 15 --resources mem_mb=200000 --cluster "sbatch -t 120 -J metap -p bmm -n 1 -N 1 -c {threads} --mem={resources.mem_mb}" -k -n
```
