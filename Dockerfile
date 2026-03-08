FROM ubuntu:24.04

LABEL maintainer="Yuanfeng Zhang"
LABEL description="MethFlow: snakemake pipeline for NGS-based methylation analysis"

ENV DEBIAN_FORNTED=noninteractive
ENV CORES=64

SHELL ["/bin/bash", "-c"]

# Use Alibaba APT source for users in China
# Remove line 13-14 if you are outside of China
RUN sed -i 's@//.*archive.ubuntu.com@//mirrors.aliyun.com/@g' /etc/apt/sources.list.d/ubuntu.sources && \
    sed -i 's@//.*security.ubuntu.com@//mirrors.aliyun.com/@g' /etc/apt/sources.list.d/ubuntu.sources

RUN apt-get update -y && \
    apt-get install -y build-essential \
        aria2 git curl perl cmake ca-certificates \
        automake autoconf libtool pkg-config \
        zlib1g-dev libbz2-dev libz-dev libgsl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /opt && chmod 777 -R /opt

#! create a normal user and install other packages
#! as the normal user to prevent permission issues
RUN printf 'CREATE_MAIL_SPOOL=no' >> /etc/default/useradd \
    && mkdir -p /home/snake \
    && groupadd snake \
    && useradd snake -g snake -d /home/snake \
    && chown -R snake:snake /home/snake \
    && chown -R snake:snake /opt \
    && usermod -aG sudo snake

USER snake

RUN aria2c -c \
        "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" \
        -d /tmp -o miniforge.sh && \
    bash /tmp/miniforge.sh -b -f -p /opt/miniforge && \
    rm /tmp/miniforge.sh

# RUN bash resources/Miniforge3-Linux-x86_64.sh -b -f -p /opt/miniforge

ENV PATH=/opt/miniforge/bin:$PATH

RUN mamba shell init && \
    source ~/.bashrc

RUN mamba create -n snakemake \
        -c conda-forge -c bioconda snakemake htslib -y && \
    mamba create -n cpython \
        -c conda-forge python=3.12 pybind11 -y && \
    mamba create -n genomic_tools -c bioconda -c conda-forge \
        liblzma samtools -y && \
    mamba create -n reckoner -c conda-forge \
        python=3.9.22 pybind11=2.13.6 -y

SHELL ["conda", "run", "-n", "snakemake", "/bin/bash", "-c"]

# Install abismal

ENV CONDA_PREFIX=/opt/miniforge/envs/genomic_tools
COPY --chown=snake:snake . /opt/MethFlow
WORKDIR /opt/MethFlow
RUN tar -xvzf resources/abismal-3.2.4.tar.gz -C resources/  && \
    cd resources/abismal-3.2.4 && \
    mkdir build && cd build && ../configure \
      --prefix /opt/MethFlow/abismal \
      CPPFLAGS="-I${CONDA_PREFIX}/include" \
      LDFLAGS="-L${CONDA_PREFIX}/lib" \
      LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH && \
    make -j ${CORES} && make install && \
    rm ../../abismal-3.2.4.tar.gz && \
    rm -rf ../../abismal-3.2.4

# Install rastair
WORKDIR /opt/MethFlow/resources
RUN tar -xzvf rastair-v2.0.0-x86_64-unknown-linux-gnu.tar.gz && \
    mv rastair-v2.0.0-x86_64-unknown-linux-gnu rastair && \
    rm rastair-v2.0.0-x86_64-unknown-linux-gnu.tar.gz

ENV PATH="$(pwd)/rastair:${PATH}"

# Install BatMeth2
WORKDIR /opt/MethFlow/resources/BatMeth2
RUN ./configure --prefix=$(pwd) && \
    make -j ${CORES} && \
    make install

ENV PATH="$(pwd)/bin:${PATH}"

# Install FAME
WORKDIR /opt/MethFlow/resources/FAME
RUN make -j ${CORES}

ENV PATH="$(pwd):${PATH}"

# Install hisat-3n
WORKDIR /opt/MethFlow/resources/hisat-3n
RUN make -j ${CORES}

ENV PATH="$(pwd):${PATH}"

# Install hisat2
WORKDIR /opt/MethFlow/resources/hisat2
RUN make -j ${CORES}

ENV PATH="$(pwd):${PATH}"

# Install Whisper
WORKDIR /opt/MethFlow/resources/Whisper/src
RUN chmod +x whisper && \
    chmod +x whisper-index

# Install BioSeqZip
WORKDIR /opt/MethFlow/resources/BioSeqZip
RUN mkdir build && cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release && \
    make -j ${CORES}
ENV PATH="$(pwd)/apps:${PATH}"

# Install RECKONER2
WORKDIR /opt/MethFlow/resources/RECKONER
# Skip py_kmc_api (Python API incompatible with newer GCC/Python)
RUN if grep -q "py_kmc_api" kmc_dir/makefile; then \
        sed -i 's/ py_kmc_api$//' kmc_dir/makefile; \
    fi && \
    mamba run -n reckoner make -j ${CORES}
ENV PATH="$(pwd)/bin:${PATH}"

# Install GEM3 (gem3-mapper)
# ! Watch out! GEM3 `make -j ${cores}` will pop an error.
# gem3-mapper with CUDA has not been tested yet.
WORKDIR /opt/MethFlow/resources/gem3-mapper
RUN ./configure && \
    make HAVE_CUDA=0 all
ENV PATH="$(pwd)/bin:${PATH}"

# Use snakemake to install all the conda env
WORKDIR /opt/MethFlow
RUN mkdir -p input && \
    touch input/BS1_BL_1.R1.fq.gz input/BS1_BL_1.R2.fq.gz
RUN snakemake --snakefile fq2bedgraph.smk \
        --config sample_sheet=utils/conda_trigger.csv \
        --cores 1 --use-conda --conda-create-envs-only \
        --verbose --software-deployment-method conda

# Use your own runtime_config.yaml mount its folder as /data instead
RUN find . -name "*.smk" \
    -exec sed \
    -i 's|config/runtime_config.yaml|/data/runtime_config.yaml|g' {} +
