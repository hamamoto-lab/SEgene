# base image
FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive


RUN apt-get update &&  \
    apt-get install -y \
        wget \ 
        build-essential \
        libcurl4-openssl-dev \
        fort77 \
        xorg-dev \
        liblzma-dev \
        libblas-dev \
        gfortran \
        gcc-multilib \
        gobjc++ \
        aptitude \
        libbz2-dev \
        libpcre2-dev \
        locales && \
    aptitude install libreadline-dev -y

RUN mkdir /tmp_work/
WORKDIR /tmp_work


RUN wget https://cran.r-project.org/src/base/R-4/R-4.2.2.tar.gz \
    && tar xvf R-4.2.2.tar.gz \
    && cd R-4.2.2 \
    && ./configure \
    && make -j12 \
    && make install \
    && cd ../../ \
    && rm -rf /tmp_work

RUN R -e \
    "install.packages(c('BiocManager', \
        'data.table', \
        'openxlsx', \
        'optparse', \
        'pbmcapply', \
        'stringr'), \
        repos='https://cran.r-project.org/'); \
        BiocManager:::install(c('GenomicRanges', 'rhdf5'))"


WORKDIR /opt
ENV JULIA_DEPOT_PATH=/opt/julia-1.8.3/depot
RUN  wget https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.3-linux-x86_64.tar.gz \
        && tar xvf julia-1.8.3-linux-x86_64.tar.gz \
        && chown -R root:root julia-1.8.3 \
        && ln -s /opt/julia-1.8.3/bin/julia /usr/local/bin/julia \
        && rm /opt/julia-1.8.3-linux-x86_64.tar.gz  \
        && mkdir -p /opt/julia-1.8.3/depot \
        && chmod -R 755 /opt/julia-1.8.3/depot

RUN julia -e 'using Pkg; Pkg.add("ArgParse")'
RUN julia -e 'using Pkg; Pkg.add("HDF5")'
RUN julia -e 'using Pkg; Pkg.add("RData")'
RUN julia -e 'using Pkg; Pkg.add("StatsBase")'

RUN mkdir /opt/P2GL

WORKDIR /opt/P2GL

RUN mkdir /opt/P2GL/output
RUN mkdir /data

COPY ../src/peakToGeneLinks* /opt/P2GL/

CMD ["bash"]