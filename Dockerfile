FROM ubuntu

LABEL Description="CosmoStat"

SHELL ["/bin/bash", "-c"]
WORKDIR /workdir

ARG DEBIAN_FRONTEND=noninteractive
ARG CC=gcc-9
ARG CXX=g++-9

RUN apt-get update && \
    apt-get install -y autoconf automake libtool pkg-config libgl1-mesa-glx && \
    apt-get install -y gcc-9 g++-9 && \
    apt-get install -y cmake git wget && \
    apt-get install -y libarmadillo-dev && \
    apt-get install -y libcfitsio-dev && \
    apt-get install -y libfftw3-dev && \
    apt-get install -y libgsl-dev && \
    apt-get install -y libsharp-dev && \
    apt-get install -y libhealpix-cxx-dev && \
    apt-get install -y healpy-data && \
    apt-get clean

RUN wget https://github.com/catchorg/Catch2/archive/refs/tags/v3.1.0.tar.gz && \
    tar -xvf v3.1.0.tar.gz && \
    cd Catch2-3.1.0 && \
    cmake -Bbuild -H. -DBUILD_TESTING=OFF && \
    cmake --build build/ --target install

ENV HEALPIX /usr/share/healpy

RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /miniconda && \
    rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH /miniconda/bin:${PATH}

RUN conda create -n cosmostat python=3.10 pip -y

ENV PATH /miniconda/envs/cosmostat/bin:${PATH}

RUN python -m pip install jupyter

COPY . /home

RUN cd /home && \
    python -m pip install .

ENV LD_LIBRARY_PATH /miniconda/envs/cosmostat/lib

RUN echo -e '#!/bin/bash\njupyter notebook --port=8888 --no-browser --ip=0.0.0.0 --allow-root' > /usr/bin/notebook && chmod +x /usr/bin/notebook
