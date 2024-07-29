# GSAP: Genome Sequence Assembly-and-Annotation Pipeline

This Python package automates the assembly, annotation, and variant analysis of genome sequence data from NGS short read sequences.

---

<p align="center">
  <a href="https://github.com/sungshic/gsap/actions/workflows/ci.yml?query=branch%3Amain">
    <img src="https://img.shields.io/github/actions/workflow/status/sungshic/gsap/ci.yml?branch=main&label=CI&logo=github&style=flat-square" alt="CI Status" >
  </a>
  <a href="https://gsap.readthedocs.io">
    <img src="https://img.shields.io/readthedocs/gsap.svg?logo=read-the-docs&logoColor=fff&style=flat-square" alt="Documentation Status">
  </a>
  <a href="https://codecov.io/gh/sungshic/gsap">
    <img src="https://img.shields.io/codecov/c/github/sungshic/gsap.svg?logo=codecov&logoColor=fff&style=flat-square" alt="Test coverage percentage">
  </a>
</p>
<p align="center">
  <a href="https://python-poetry.org/">
    <img src="https://img.shields.io/endpoint?url=https://python-poetry.org/badge/v0.json" alt="Poetry">
  </a>
  <a href="https://github.com/astral-sh/ruff">
    <img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json" alt="Ruff">
  </a>
  <a href="https://github.com/pre-commit/pre-commit">
    <img src="https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white&style=flat-square" alt="pre-commit">
  </a>
</p>
<p align="center">
  <a href="https://pypi.org/project/gsap/">
    <img src="https://img.shields.io/pypi/v/gsap.svg?logo=python&logoColor=fff&style=flat-square" alt="PyPI Version">
  </a>
  <img src="https://img.shields.io/pypi/pyversions/gsap.svg?style=flat-square&logo=python&amp;logoColor=fff" alt="Supported Python versions">
  <img src="https://img.shields.io/pypi/l/gsap.svg?style=flat-square" alt="License">
</p>

---

**Documentation**: <a href="https://gsap.readthedocs.io" target="_blank">https://gsap.readthedocs.io </a>

**Source Code**: <a href="https://github.com/sungshic/gsap" target="_blank">https://github.com/sungshic/gsap </a>

---

## Introduction

NGS devices, such as those made by Illumina, generate short read sequences that are collated into raw, unassembled sequence data formats like FASTQ. Converting this raw data into a format ready for high-level bioinformatic analyses involves multiple laborious steps that are prone to error if performed manually. This Python package implements a software pipeline designed to automate the assembly and annotation of NGS sequence data from short sequence reads. The pipeline also performs variant analysis, annotating the assembled sequence for SNPs and INDELs against a chosen reference genome.

## Overview of the Pipeline

The pipeline consists of four main parts: sequence assembly, sequence pre-processing, variant analysis, and sequence annotation.

### Sequence Assembly

![image of the pipeline overview](https://github.com/sungshic/gsap/blob/main/docs/_static/assets/pipeline_overview1.png?raw=true)

## Evolving Research Computing: From Legacy Code to Modern Software Engineering Practices in research

The initial code base for this project dates back to around 2015. The code here was manually forked and migrated from <a href="https://bitbucket.org/sungshic/genomeassemblypipeline" target="_blank">this Bitbucket repository</a>. It was originally a rudimentary implementation of a genome assembly pipeline written in Python 2.7, now refactored in Python 3.11 by adopting popular modern practices in software engineering, including meticulous documentation, containerisation, test-driven development, and continuous integration. The primary purpose of the initial code base was not to achieve high software engineering standards but to quickly assemble sequence samples of about 20 bacterial strains being studied at the Centre for Bacterial Cell Biology, Newcastle University.

Over the past decade, computer science has made inroads into many research fields traditionally not oriented toward computing. The fields of biological science, systems biology, and synthetic biology, in particular, have seen a significant adoption of research computing practices. Yet, software in these fields is often considered one-off disposable-ware, used primarily to support research output in the form of publications.

One of the primary purposes of the GSAP repository is to showcase the need for adopting agile software engineering practices in research computing. These practices make it easier to develop, maintain, and distribute software, ultimately fostering long-term computational reproducibility in scientific research.

## Strategic choice of GSAP's tech stack to maximise computational reproducibility

Making the code base of a research output open source is certainly an important first step in addressing the concern of reproducibility in scientific research and research computing. However, source code alone is often insufficient to reproduce results. In fact, there are three essential components to computational reproducibility: source code, data, and the runtime environment. Having this trio of assets ready to provision the execution of computational tasks is conceptually similar to providing Infrastructure as Code (IaC).

This repository exemplifies the open-source provisioning of a computerized scientific pipeline, based on the following choices in its tech stack:

Python:

- One of the most versatile languages for integrating libraries written in multiple programming languages.
- Strong and widespread community support for libraries in bioinformatics, AI, and backend microservice frameworks.

Open Science Framework:

- Serves as a repository to store large datasets in research.
- A non-profit initiative with ample storage and bandwidth quotas far exceeding those of GitHub and GitLab.

Docker:

- The de facto community standard in containerization.
- Widespread community adoption and support for providing IaC at scale.
- Compatible with massively scalable cloud architectures such as Kubernetes.

Best Practices in Software Engineering for Test-Driven Agile Development and the Adoption of Continuous Integration:

- Poetry for Python package management.
- Pytest for Python unit testing.
- Pre-commit for managing the linting, static types and error checks
- Codespell, Ruff & Mypy for applying various rule-based validity checks via pre-commit.
- GitHub Actions for build, test, and release workflows.

By integrating these technologies and best practices, this open source project can provide a complete and standardised runtime environment along with necessary data and code to ensure the general reproducibility of computational research.

## Installation

Install this via docker and docker-compose by:

```bash
docker compose up -d
```

It would launch a container in your local computer, by downloading the latest pre-built GSAP image from the Docker Hub.

Alternatively, this repo contains all the codes necessary to rebuild and run a Docker image locally by executing the following commands:

```bash
docker compose -f docker-compose.dev.yml build
docker compose -f docker-compose.dev.yml up -d
```

Please note, if running the build process on Apple silicon computers, the following general setting of Docker Desktop should be opted out (FYI: it is opted in by default as of Docker Desktop version 4.30.0 (149282)).

"Use Rosetta for x86_64/amd64 emulation on Apple Silicon"

Once a GSAP container is launched, please check if the container is running without error like so:

```bash
docker ps
```

Then enter the container's command line space by executing:

```bash
docker exec -it $(docker ps -aqf "ancestor=sungshic/gsap") /bin/bash
```

Or, by using the following, if the image was rebuilt locally:

```bash
docker exec -it $(docker ps -aqf "ancestor=gsap") /bin/bash
```

## Contributors ✨

Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)):

<!-- prettier-ignore-start -->
<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- markdownlint-disable -->
<!-- markdownlint-enable -->
<!-- ALL-CONTRIBUTORS-LIST:END -->
<!-- prettier-ignore-end -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!
