# Snakemake workflow: `gnnotator`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/r1cheu/gnnotator/actions/workflows/main.yml/badge.svg?branch=main)](https://github.com/r1cheu/gnnotator/actions?query=branch%3Amain+workflow%3ATests)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with apptainer](https://img.shields.io/badge/run%20with-apptainer-52307c?labelColor=000000&logo=apptainer)](https://apptainer.org/)

[![workflow catalog](https://img.shields.io/badge/Snakemake%20workflow%20catalog-darkgreen)](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/r1cheu/gnnotator)

A Snakemake workflow for `genome annotation`

- [Snakemake workflow: `gnnotator`](#snakemake-workflow-gnnotator)
  - [Usage](#usage)
  - [Deployment options](#deployment-options)
  - [Authors](#authors)
  - [References](#references)
  - [TODO](#todo)

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/r1cheu/gnnotator).

Detailed information about input data and workflow configuration can also be found in the [`config/README.md`](config/README.md).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository or its DOI.

## Deployment options

To run the workflow from command line, change the working directory. This workflow need to run with **conda** and **apptainer** / **singularity**.

```bash
cd gnnotator
```

Provide your own input data in the `data/` directory. The default config file is located at `config/config.yml`.
Before running the complete workflow, you can perform a dry run using:

```bash
snakemake --dry-run
```

You can prepare the environments with:

```bash

snakemake --sdm conda apptainer --conda-create-envs-only

```

Run the workflow with 20 cores:

```bash
snakemake --cores 20 --sdm conda apptainer
```

It's recommended to run the workflow with slurm, and do not forget to change account in the `slurm/config.yaml` file if you are using slurm.

```bash
snakemake --sdm conda apptainer --profile slurm

```

Use the `--ri` to rerun the uncompleted jobs:

```bash
snakemake --sdm conda apptainer --profile slurm --ri

```

## Authors

- RuLei Chen
  - Develop the snakemake workflow
  - Center for Excellence in Molecular Plant Sciences
  - [ORCID](https://orcid.org/0009-0000-9645-0951)
  - [github](https://github.com/r1cheu)

- ZhouLin Gu
  - Design the original pipeline
  - Center for Excellence in Molecular Plant Sciences
  - [ORCID](https://orcid.org/0000-0001-7732-9515)

## References

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, 10, 33, **2021**. https://doi.org/10.12688/f1000research.29032.2.

## TODO
