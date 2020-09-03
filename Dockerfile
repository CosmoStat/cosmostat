FROM ubuntu

LABEL Description="CosmoStat"

ARG DEBIAN_FRONTEND=noninteractive
ARG CC=gcc-9
ARG CXX=g++-9

RUN apt-get update && \
    apt-get install -y autoconf automake libtool pkg-config && \
    apt-get install -y gcc-9 g++-9 && \
    apt-get install -y git cmake wget && \
    apt-get install -y libcfitsio-dev libhealpix-cxx-dev && \
    apt-get clean

RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    chmod +x miniconda.sh && \
    ./miniconda.sh -b -p $HOME/miniconda && \
    export PATH=$HOME/miniconda/bin:$PATH && \
    conda init bash && \
    conda create -n cosmostat python=3.8 -y


RUN cd home && \
    git clone https://github.com/CosmoStat/cosmostat && \
    cd cosmostat && \
    conda activate cosmostat && \
    python setup.py install
