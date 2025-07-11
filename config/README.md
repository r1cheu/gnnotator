## Workflow overview

This workflow is a workflow for `genome annotation`.
The workflow is built using [snakemake](https://snakemake.readthedocs.io/en/stable/) and consists of the following steps:

1. Download genome reference from NCBI
2. Validate downloaded genome (`python` script)
3. Simulate short read sequencing data on the fly (`dwgsim`)
4. Check quality of input read data (`FastQC`)
5. Collect statistics from tool output (`MultiQC`)

## Running the workflow

### Input data

Modify the samplesheet file `config/samples.tsv` and prepare the data for the workflow.

| id      |
| ------- |
| 1GS-002 |

Then, create two directocries data/assembly and data/rnaseq, note the tissue of RNA-seq should match that in config.yaml.
e.g.

```bash
data
├── assembly
│   └── 1GS-002.fa
└── rnaseq
    ├── 1GS-002_fringe_R1.fastq
    ├── 1GS-002_fringe_R2.fastq
    ├── 1GS-002_leaf_R1.fastq
    ├── 1GS-002_leaf_R2.fastq
    ├── 1GS-002_root_R2.fastq
    ├── 1GS-002_seedling_R1.fastq
    └── 1GS-002_seedling_R2.fastq
```

### Parameters

Change config.yaml to set the parameters for the workflow.

E.g. change the `rna_seq_tissue` to the tissue of RNA-seq data you have.

```yaml
rna_seq_tissue:
  - LEAF
  - ROOT
```

For parameters of pasa_config and evm_weights, it is recommended to use the default values, unless you know what you are doing.

For user who use slurm, change the slurm account in `slurm/config.yaml`
E.g.

```yaml
default-resources:
  slurm_account: "your account"
```

And see [Snakemake executor plugin: slurm](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) for documentation.
