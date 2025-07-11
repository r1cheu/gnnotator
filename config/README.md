## Workflow overview

This workflow is a workflow for `genome annotation`.
The workflow is built using [snakemake](https://snakemake.readthedocs.io/en/stable/) and consists of the following steps:

1. Annotate and mask interspersed repeats and low complexity DNA sequences (`RepeatMasker`) in assemblies
2. Conduct transcript-based gene predictions, including RNA-seq reads alignment (`hisat2`) and alignment assembly pipeline (`StringTie`, `Trinity`, `PASA`) 
3. Conduct ab initio gene predictions based on `Fgenesh`, `SNAP` and `Eviann`
4. Conduct homology-based annotation by `miniprot`
5. Combine ab intio gene predictions and protein and transcript alignments into weighted consensus gene structures by `EVM`

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
