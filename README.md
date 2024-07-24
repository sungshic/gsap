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

## The pipeline overview

![image of the pipeline overview](https://github.com/sungshic/gsap/blob/main/docs/_static/assets/pipeline_overview1.png?raw=true)

## Installation

Install this via pip (or your favourite package manager):

`pip install gsap`

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
