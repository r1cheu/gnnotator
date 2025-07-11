# Snakemake workflow: `gnnotator`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/r1cheu/gnnotator/workflows/Tests/badge.svg?branch=main)](https://github.com/r1cheu/gnnotator/actions?query=branch%3Amain+workflow%3ATests)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with apptainer](https://img.shields.io/badge/run%20with-apptainer-52307c?labelColor=000000&logo=apptainer)](https://apptainer.org/)

[![workflow catalog](https://img.shields.io/badge/Snakemake%20workflow%20catalog-darkgreen)](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/r1cheu/gnnotator)

A Snakemake workflow for `genome annotation`

- [Snakemake workflow: `gnnotator`](#snakemake-workflow-name)
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

Adjust options in the default config file `config/config.yml`.
Before running the complete workflow, you can perform a dry run using:

```bash
snakemake --dry-run
```

You can prepare the environments with:

```bash

snakemake --cores 20 --sdm conda apptainer --conda-create-envs-only

```

To run the workflow with **apptainer** / **singularity**, add a link to a container registry in the `Snakefile`, for example `container: "oras://ghcr.io/<user>/<repository>:<version>"` for Github's container registry.
Run the workflow with:

```bash
snakemake --cores 20 --sdm conda apptainer
```

It's recommended to run the workflow with slurm

## Authors

- RuLei Chen
  - develop the snakemake workflow
  - Affiliation
  - ORCID profile
  - home page

- ZhouLin Gu
  - design the original pipeline
  - Affiliation
  - ORCID profile
  - home page

## References

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, 10, 33, **2021**. https://doi.org/10.12688/f1000research.29032.2.

## TODO
