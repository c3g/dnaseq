FROM ubuntu:21.10

RUN apt-get update -y && apt-get upgrade -y && apt-get install -y build-essential wget

# Skewer
RUN apt-get install -y skewer

# BWA
RUN apt-get install -y bwa samtools bamtools tabix

# Sambamba
RUN mkdir -p /tmp/sambamba \
    && cd /tmp/sambamba \
    && wget --quiet https://github.com/biod/sambamba/releases/download/v0.8.1/sambamba-0.8.1-linux-amd64-static.gz \
    && gunzip *.gz \
    && chmod +x sambamba* \
    && mv sambamba* /usr/local/bin/sambamba \
    && rm -r /tmp/sambamba

# FastQC
RUN apt-get install -y fastqc

# GATK
RUN apt-get install -y python unzip \
    && mkdir -p /tmp/gatk \
    && cd /tmp/gatk \
    && wget --quiet https://github.com/broadinstitute/gatk/releases/download/4.2.2.0/gatk-4.2.2.0.zip \
    && unzip gatk-*.zip \
    && rm *.zip \
    && mv gatk* /opt/gatk
ENV PATH="/opt/gatk:$PATH"

# VT
RUN apt-get install -y zlib1g-dev \
    &&  mkdir -p /tmp/vt \
    && cd /tmp/vt \
    && wget --quiet https://github.com/atks/vt/archive/refs/tags/0.57.tar.gz \
    && tar -xzvf *.tar.gz \
    && rm  *.tar.gz \
    && cd /tmp/vt/vt* \
    && make \
    && mv vt /usr/local/bin \
    && rm -r /tmp/vt

# SnpEff
RUN apt-get install -y unzip \
    && cd /opt \
    && wget --quiet https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip \
    && unzip snpEff*.zip \
    && rm snpEff*.zip
ENV SNPEFF_HOME=/opt/snpEff

# VCFtools
RUN apt-get install -y pkg-config \
    && mkdir -p /opt/vcftools \
    && cd /opt/vcftools \
    && wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz \
    && tar -xzvf vcftools*.tar.gz \
    && ln -s vcftools-0.1.16 latest \
    && cd latest \
    && ./configure && make && make install

# Fastp
RUN apt-get install -y fastp
