FROM python:3.8-buster@sha256:7e7f4c5508b85268a93b573566c8eb321a6fdb466e3b60c663a42300c73a7400

LABEL maintainer="Mei Knudson"

# To prevent time zone prompt
ENV DEBIAN_FRONTEND=noninteractive
ENV SAMTOOLS_VERSION 1.9

# Install softwares from apt repo
RUN apt-get update && apt-get install -y \
    autoconf \
    binutils \
    build-essential \
    gcc \
    git \
    libcurl4-openssl-dev \
    libjpeg-dev \
    liblz4-dev \
    liblzma-dev \
    libncurses5-dev \
    libbz2-dev \
    python3-pip \
    wget \
    zlib1g-dev &&\
    rm -rf /var/lib/apt/lists/*

# Install samtools 1.9
RUN git clone --branch ${SAMTOOLS_VERSION} --single-branch https://github.com/samtools/samtools.git && \
    git clone --branch ${SAMTOOLS_VERSION} --single-branch https://github.com/samtools/htslib.git && \
    cd samtools && make && make install && cd ../ && rm -rf samtools* htslib*

# Install packages for python3 scripts
RUN python3 -m pip install numpy pysam TOBIAS

# Create and setup new user
ENV USER=tobias
WORKDIR /home/$USER
RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

ENV PYTHONPATH="/usr/local/python:$PYTHONPATH"

USER ${USER}
