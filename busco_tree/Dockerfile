# Use NVIDIA image for GPU compatibility in ASTRAL, but...
# ASTRAL runs quickly enough to where we don't need it yet
#FROM nvidia/cuda:10.0-cudnn7-devel-ubuntu18.04
FROM ubuntu:18.04

#install libs/dependencies
#manual install of Cmake for compatibility with Racon
RUN apt-get update && apt-get install -y \
    curl \
    default-jdk \
    default-jre \
    gcc \
    git \
    g++ \
    make \
    pigz \
    wget \
    unzip \
    zip \
    zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*
RUN mkdir /tools

# Install MAFFT v7.453
RUN wget -q https://mafft.cbrc.jp/alignment/software/\
mafft-7.453-without-extensions-src.tgz \
 && tar zxf mafft-7.453-without-extensions-src.tgz \
 && rm mafft-7.453-without-extensions-src.tgz \
 && cd mafft-7.453-without-extensions/core \
 && make clean && make -j && make install \
 && cd / && rm -r mafft-7.453-without-extensions

# Install RAxML-NG
RUN wget -q https://github.com/amkozlov/raxml-ng/releases/download/\
0.9.0/raxml-ng_v0.9.0_linux_x86_64.zip \
 && unzip -p raxml-ng_v0.9.0_linux_x86_64.zip raxml-ng > raxml-ng \
 && rm raxml-ng_v0.9.0_linux_x86_64.zip \
 && chmod +x raxml-ng && mv raxml-ng /usr/local/bin

# Install ASTRAL-MP
RUN git clone --branch 'MP' https://github.com/smirarab/ASTRAL.git \
 && cd ASTRAL \
 && bash ./make.sh \
 && cp astral.5.14.7.jar /tools \
 && cp -r ./lib /tools

#udo apt install ocl-icd-libopencl1
#sudo apt install ocl-icd-opencl-dev

