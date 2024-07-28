# GSAP: Genome Sequence Assembly-and-Annotation Pipeline

This Python package automates the assembly, annotation, and variant analysis of genome sequence data from NGS short read seqeunces.

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

###

## Installation

Install this via pip (or your favourite package manager):

`pip install gsap`

## Evolving Research Computing: From Legacy Code to Modern Software Engineering Practices in research

The initial code base for this project dates back to around 2015. The code here was manually forked and migrated from <a href="https://bitbucket.org/sungshic/genomeassemblypipeline" target="_blank">this Bitbucket repository</a>. It was originally a rudimentary implementation of a genome assembly pipeline written in Python 2.7, now refactored in Python 3.11 by adopting popular modern practices in software engineering, including meticulous documentation, containerisation, test-driven development, and continuous integration. The primary purpose of the initial code base was not to achieve high software engineering standards but to quickly assemble sequence samples of about 20 bacterial strains being studied at the Centre for Bacterial Cell Biology, Newcastle University.

Over the past decade, computer science has made inroads into many research fields traditionally not oriented toward computing. The fields of biological science, systems biology, and synthetic biology, in particular, have seen a significant adoption of research computing practices. Yet, software in these fields is often considered one-off disposable-ware, used primarily to support research output in the form of publications.

Another purpose of the GSAP repository is to showcase the need for adopting agile software engineering practices in research computing. These practices make it easier to develop, maintain, and distribute software, ultimately fostering long-term computational reproducibility in scientific research.

## Contributors ✨

Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)):

<!-- prettier-ignore-start -->
<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- markdownlint-disable -->
<!-- markdownlint-enable -->
<!-- ALL-CONTRIBUTORS-LIST:END -->
<!-- prettier-ignore-end -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!

## Credits

This package was created with
[Copier](https://copier.readthedocs.io/) and the
[browniebroke/pypackage-template](https://github.com/browniebroke/pypackage-template)
project template.
