(installation)=

# Installation

Start by first cloning the repository:

```bash
git clone https://github.com/sungshic/gsap
cd gsap
```

In order to execute the pipeline with a custom data set, a specific directory structure is needed, as expected by GSAP's docker container. These directories serve as the interface between the local computer and the GSAP's runtime environment:

```bash
mkdir -p data_input/sample
mkdir -p data_input/refseq
mkdir -p data_output
```

Using GSAP's unit testing data as example, the relevant directories can be populated with the data set like so:

```bash
curl -sSL https://osf.io/xr48c/download -o data_input/sample/SP2_S6_L001_R1_001.fastq.gz
curl -sSL https://osf.io/9gya4//download -o data_input/sample/SP2_S6_L001_R2_001.fastq.gz
curl -sSL https://osf.io/d4gjv/download -o data_input/refseq/AL009126_v11.tar.gz
tar -xzf data_input/refseq/AL009126_v11.tar.gz -C data_input/refseq/
```

Download and spawn the latest GSAP runtime environment:

```bash
docker pull sungshic/gsap:latest
docker compose up -d
```

While the latest pre-built GSAP image is available from the Docker Hub, this repo contains all the codes necessary to rebuild and run a Docker image locally by executing the following commands:

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

Next, see the {ref}`section about usage <usage>` to see how to use it.
