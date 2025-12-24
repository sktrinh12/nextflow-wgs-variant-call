FROM debian:bullseye-slim

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    wget \
    procps \
    bzip2 \
    ca-certificates \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Install Miniforge
ENV CONDA_DIR=/opt/conda
RUN wget --quiet https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O ~/miniforge.sh && \
    /bin/bash ~/miniforge.sh -b -p $CONDA_DIR && \
    rm ~/miniforge.sh && \
    $CONDA_DIR/bin/conda clean -afy

ENV PATH=$CONDA_DIR/bin:$PATH

RUN echo "channels:" > $CONDA_DIR/.condarc && \
    echo "  - conda-forge" >> $CONDA_DIR/.condarc && \
    echo "  - bioconda" >> $CONDA_DIR/.condarc && \
    echo "channel_priority: strict" >> $CONDA_DIR/.condarc

# Install bioinformatics tools via bioconda
RUN conda install -y \
    openjdk=17 \
    r-base \
    fastqc samtools bwa gatk4 seqkit awscli && \
    conda clean -afy

WORKDIR /data

CMD ["/bin/bash"]
