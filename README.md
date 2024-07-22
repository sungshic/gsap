# Genome Sequence Assembly Pipeline

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

NGS devices generate short sequence reads collated together into raw unassembled
sequence data formats, such as FASTQ. There are a series of steps involved in
converting the raw data into some form ready to be consumed by other high-
level bioinformatic analyses. This Python package is a software pipeline implementation that automates
the assembly and the annotation of NGS sequence data from short sequence reads. The pipeline also performs
variant analysis to annotate the assembled sequence about any SNPs and INDELs with respect to a reference genome of choice.

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
