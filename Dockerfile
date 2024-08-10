FROM ubuntu

LABEL Description="CosmoStat"

SHELL ["/bin/bash", "-c"]
WORKDIR /workdir

ARG DEBIAN_FRONTEND=noninteractive
ARG CC=gcc-9
ARG CXX=g++-9

RUN apt-get update && \
    apt-get install -y autoconf automake libtool pkg-config && \
    apt-get install -y gcc-9 g++-9 && \
    apt-get install -y cmake git wget && \
    apt-get install -y libarmadillo-dev && \
    apt-get install -y libcfitsio-dev && \
    apt-get install -y libfftw3-dev && \
    apt-get install -y libgl1 && \
    apt-get install -y libgsl-dev && \
    apt-get install -y libsharp-dev && \
    apt-get install -y libhealpix-cxx-dev && \
    apt-get install -y healpy-data && \
    apt-get install -y python3 python3-pip && \
    apt-get clean

ENV HEALPIX=/usr/share/healpy

RUN PIP_VERSION=$(python3 -m pip --version | awk '{print $2}') && \
    PIP_OPTS="" && \
    if [[ $(echo -e "23.1\n$PIP_VERSION" | sort -V | head -n1) == "23.1" && "$PIP_VERSION" != "23.1" ]]; then \
        PIP_OPTS="--break-system-packages"; \
    fi && \
    echo "export PIP_OPTS=$PIP_OPTS" >> ~/.bashrc && \
    export PIP_OPTS=$PIP_OPTS

RUN echo $PIP_OPTS

RUN python3 -m pip install jupyter $PIP_OPTS

COPY . /home

RUN cd /home && \
    python3 -m pip install . $PIP_OPTS

RUN echo -e '#!/bin/bash\njupyter notebook --port=8888 --no-browser --ip=0.0.0.0 --allow-root' > /usr/bin/notebook && chmod +x /usr/bin/notebook
