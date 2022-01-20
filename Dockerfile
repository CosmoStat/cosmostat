FROM ubuntu

LABEL Description="CosmoStat"
WORKDIR /home

SHELL ["/bin/bash", "-c"]

ARG DEBIAN_FRONTEND=noninteractive
ARG CC=gcc-9
ARG CXX=g++-9

RUN apt-get update && \
    apt-get install -y autoconf automake libtool pkg-config libgl1-mesa-glx && \
    apt-get install -y gcc-9 g++-9 && \
    apt-get install -y cmake git pkg-config wget && \
    apt-get install -y libarmadillo-dev && \
    apt-get install -y libcfitsio-dev && \
    apt-get install -y libfftw3-dev && \
    apt-get install -y libgsl-dev && \
    apt-get install -y libsharp-dev && \
    apt-get install -y libhealpix-cxx-dev && \
    apt-get clean

RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /miniconda && \
    rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH /miniconda/bin:${PATH}

RUN conda create -n cosmostat python=3.9 -y && \
    conda install jupyter -y

ENV PATH /miniconda/envs/cosmostat/bin:${PATH}

RUN mkdir /cosmostat && \
    cd /cosmostat && \
    git clone https://github.com/sfarrens/cosmostat

RUN cd /cosmostat/cosmostat && \
    git checkout gcc_build && \
    source activate cosmostat && \
    python setup.py install

ENV PATH /cosmostat/cosmostat/build/bin:${PATH}
