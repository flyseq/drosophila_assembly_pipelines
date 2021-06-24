# This Dockerfile contains the environment used for Nanopore base
# calling, assembly, and polishing. BLAST is not included due to the
# disk space required for the nt database, and RepeatMasker is not
# included due to our usage of the RepBase RepeatMasker library

# Use CUDA 10.0 for tensorflow compatibility (in medaka)
FROM nvidia/cuda:10.0-cudnn7-devel-ubuntu16.04

# install libs/dependencies
# manual install of Cmake for compatibility with Racon
RUN apt-get update && apt-get install -y \
    bzip2 \
    curl \
    default-jdk \
    default-jre \
    gcc \
    git \
    g++ \
    libboost-all-dev \
    libbz2-dev \
    libcurl4-openssl-dev \
    libffi-dev \
    libhdf5-cpp-11 \
    liblzma-dev \
    libncurses5-dev \
    libssl-dev \
    libzmq5 \
    pigz \
    python3-all-dev \
    python3-pip \
    python-virtualenv \
    wget \
    unzip \
    zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*
RUN git clone --branch 'v3.17.2' \
    https://gitlab.kitware.com/cmake/cmake.git \
 && cd cmake \
 && ./bootstrap \
 && make -j \
 && make install \
 && cd / && rm -r cmake

# install Guppy 3.2.4 (not latest version - but if you change the
# version number in the link to match the version you want it 
# should work
RUN wget -q https://mirror.oxfordnanoportal.com/software/analysis/\
ont_guppy_3.2.4-1~xenial_amd64.deb \
 && dpkg -i --ignore-depends=nvidia-384,libcuda1-384 \
    ont_guppy_3.2.4-1~xenial_amd64.deb \
 && rm ont_guppy_3.2.4-1~xenial_amd64.deb

# install Flye 2.6
RUN pip3 install setuptools
RUN git clone --branch 2.6 https://github.com/fenderglass/Flye.git \
 && cd /Flye && python3 setup.py install \
 && rm -r /Flye

# install CUDA-enabled Racon
RUN git clone --branch '1.4.3' --recursive \
    https://github.com/lbcb-sci/racon.git racon \
 && cd /racon \
 && mkdir -p build \
 && cd build \
 && cmake -DCMAKE_BUILD_TYPE=Release -Dracon_enable_cuda=ON .. \
 && make -j \
 && cp /racon/build/bin/racon /usr/local/bin \
 && rm -r /racon

# install minimap2
RUN git clone --branch 'v2.17' https://github.com/lh3/minimap2 \
 && cd /minimap2 \
 && make -j \
 && cp minimap2 /usr/local/bin \
 && rm -r /minimap2

# install samtools
RUN wget -q https://github.com/samtools/samtools/releases/download/\
1.12/samtools-1.12.tar.bz2 \
 && tar -xjf samtools-1.12.tar.bz2 \
 && cd /samtools-1.12 \
 && ./configure \
 && make -j \
 && make install \
 && rm -r /samtools-1.12*

# install htslib
RUN wget -q https://github.com/samtools/htslib/releases/download/\
1.12/htslib-1.12.tar.bz2 \
 && tar -xjf htslib-1.12.tar.bz2 \
 && cd /htslib-1.12 \
 && ./configure \
 && make -j \
 && make install \
 && rm -r /htslib-1.12*

# install bcftools
RUN wget -q https://github.com/samtools/bcftools/releases/download/\
1.12/bcftools-1.12.tar.bz2 \
 && tar -xjf bcftools-1.12.tar.bz2 \
 && cd /bcftools-1.12 \
 && ./configure \
 && make -j \
 && make install \
 && rm -r /bcftools-1.12*

#install medaka 0.9.1 (not latest version)
RUN wget -q https://github.com/git-lfs/git-lfs/releases/download/\
v2.10.0/git-lfs-linux-amd64-v2.10.0.tar.gz \
 && tar -zxf git-lfs-linux-amd64-v2.10.0.tar.gz \
 && ./install.sh \
 && git lfs install \
 && rm install.sh git-lfs*
RUN git clone --branch 'v0.9.1' \
    https://github.com/nanoporetech/medaka.git \
 && cd medaka \
 && sed -i 's/tensorflow/tensorflow-gpu/' requirements.txt \
 && sed -i 's/biopython/biopython==1.76/' requirements.txt \
 && make -j install

#install Purge Haplotigs
RUN mkdir -p tools \
 && cd tools \
 && git clone --branch 'v1.1.1' \
    https://bitbucket.org/mroachawri/purge_haplotigs.git \
 && ln -s /tools/purge_haplotigs/bin/purge_haplotigs \
    /usr/local/bin/purge_haplotigs

#install npscarf
RUN cd tools \
 && wget -q https://github.com/mdcao/japsa/releases/download/v1.9-2b/\
JapsaRelease.zip \
 && unzip JapsaRelease.zip \
 && rm JapsaRelease.zip \
 && cd JapsaRelease \
 && sed -i 's/--mode install/--mode install --installDir \/tools\/japsa --server yes --xmx=7G/' install.sh \
 && yes "" | ./install.sh \

#install Pilon
RUN cd tools \
 && wget -q https://github.com/broadinstitute/pilon/releases/download/\
v1.23/pilon-1.23.jar

# install sambamba
RUN wget -q https://github.com/biod/sambamba/releases/download/\
v0.8.0/sambamba-0.8.0-linux-amd64-static.gz \
 && gunzip sambamba-0.8.0-linux-amd64-static.gz \
 && chmod +x sambamba-0.8.0-linux-amd64-static \
 && cp sambamba-0.8.0-linux-amd64-static /usr/local/bin/sambamba

# install Jellyfish
RUN wget -q https://github.com/gmarcais/Jellyfish/releases/download/\
v2.3.0/jellyfish-linux \
 && chmod +x jellyfish-linux \
 && mv jellyfish-linux /usr/local/bin/jellyfish