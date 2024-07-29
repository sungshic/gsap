FROM ubuntu:24.04 as builder

ENV DEBIAN_FRONTEND noninteractive

# Poetry's configuration:
ENV POETRY_NO_INTERACTION=1 \
    POETRY_VIRTUALENVS_CREATE=true \
    POETRY_CACHE_DIR='/var/cache/pypoetry' \
    POETRY_HOME='/usr/local' \
    POETRY_VERSION=1.8.3

# REPOS
RUN apt-get -y update
RUN apt-get install -y -q software-properties-common
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt-get -y update
RUN apt-get install -y -q build-essential curl git vim tmux unzip nano python3 python3-virtualenv python3-pip
#RUN pip install poetry
RUN curl -sSL https://install.python-poetry.org -o /tmp/install-poetry.py && python3 /tmp/install-poetry.py && rm /tmp/install-poetry.py
RUN git clone https://github.com/sungshic/gsap.git
WORKDIR /gsap
RUN mkdir -p tests/data/test_ref_data
RUN poetry config virtualenvs.in-project true --local
RUN cat requirements.txt | xargs -I % sh -c 'poetry add "%"'
RUN poetry install --without dev --no-interaction

FROM ubuntu:24.04 as gsap
MAINTAINER Sunny Park <sp2307@columbia.edu>

ENV DEBIAN_FRONTEND noninteractive


# REPOS
RUN apt-get -y update
RUN apt-get install -y -q software-properties-common
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt-get -y update

# INSTALL
RUN apt-get install -y -q build-essential curl git git-lfs vim tmux unzip nano python3 python3-virtualenv python3-pip x11-xserver-utils
# matplotlib dependencies
#RUN apt-get install -y -q pkg-config libfreetype6-dev libpng-dev
# scipy dependencies
#RUN apt-get install -y -q libatlas-base-dev libatlas3-base libopenblas-base libopenblas-dev gfortran

# gsap dependencies
COPY --from=builder /gsap /gsap
WORKDIR /gsap
ENV PATH="/gsap/.venv/bin:$PATH"
#RUN git clone https://github.com/sungshic/gsap.git
RUN apt-get install -y -q openjdk-21-jre npm
RUN update-alternatives --set java /usr/lib/jvm/java-21-openjdk-amd64/bin/java
#RUN curl -sSL https://install.python-poetry.org -o /tmp/install-poetry.py && python3 /tmp/install-poetry.py && rm /tmp/install-poetry.py
RUN curl -sSL https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Linux.tar.gz | tar -xzf - -C /gsap/src/gsap/toolset/
RUN curl -sSL https://github.com/broadinstitute/gatk/releases/download/4.6.0.0/gatk-4.6.0.0.zip -o gatk-4.6.0.0.zip && unzip -d /gsap/src/gsap/toolset/ gatk-4.6.0.0.zip && rm gatk-4.6.0.0.zip
RUN apt-get install -y -q seqtk samtools picard-tools bowtie2 bcftools fastqc trimmomatic

#RUN cat requirements.txt | xargs -I % sh -c 'poetry add "%"'
#RUN poetry install --no-interaction --no-ansi

#RUN mkdir -p /gsap/tests/data/sample1
#RUN curl -sSL 'https://osf.io/xr48c/download' -o /gsap/tests/data/sample1/SP2_S6_L001_R1_001.fastq.gz
#RUN curl -sSL 'https://osf.io/9gya4//download' -o /gsap/tests/data/sample1/SP2_S6_L001_R2_001.fastq.gz
#RUN curl -sSL 'https://osf.io/auqje/download' -o /gsap/tests/data/test_ref_data.tar.gz
#RUN tar -kxzf /gsap/tests/data/test_ref_data.tar.gz -C /gsap/tests/data/

# just to keep a container running indefinitely
ENTRYPOINT ["tail", "-f", "/dev/null"]
