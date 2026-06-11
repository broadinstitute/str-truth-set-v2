# Image for aligning PacBio Iso-Seq reads to hg38 with a current pbmm2.
# Reuses the existing long-reads toolchain (gcloud, htslib/samtools 1.20, mosdepth) but
# replaces the base pbmm2 1.10.0 with pbmm2 26.1.0 installed from bioconda via micromamba.
FROM pacificbiosciences/pbmm2:1.10.0

RUN apt update && apt-get install --no-install-recommends -y \
    ca-certificates \
    pkg-config \
    dpkg-dev \
    build-essential \
    libssl-dev \
    libncurses5-dev \
    libcurl4-openssl-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    gcc \
    g++ \
    autoconf \
    make \
    unzip \
    bzip2 \
    procps \
    curl \
    wget \
    git \
    jq \
    less \
    time

# install gcloud
ENV GCLOUD_SDK_VERISON="481.0.0"
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && tar xzf google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && /google-cloud-sdk/install.sh -q \
    && /google-cloud-sdk/bin/gcloud components update --quiet \
    && rm google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz

ENV PATH=/google-cloud-sdk/bin:$PATH

# install htslib
ENV SAMTOOLS_VERSION="1.20"
RUN wget https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar xjf htslib-${SAMTOOLS_VERSION}.tar.bz2 \
    && rm htslib-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd htslib-${SAMTOOLS_VERSION} \
    && ./configure \
    && make \
    && make install \
    && make clean

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd samtools-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& make clean

# install mosdepth
RUN wget https://github.com/brentp/mosdepth/releases/download/v0.3.5/mosdepth -O /usr/local/bin/mosdepth \
	&& chmod 777 /usr/local/bin/mosdepth

# upgrade pbmm2 (base image ships 1.10.0) to a current release via micromamba + bioconda.
# pbmm2 is installed into an isolated prefix so its bundled libs/tools don't shadow the
# hand-built samtools 1.20 above; a wrapper at /usr/local/bin/pbmm2 overrides the base binary
# while keeping the real executable's $ORIGIN-based RPATH intact.
ENV PBMM2_VERSION=26.1.0
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
RUN curl -L -o /usr/local/bin/micromamba \
        https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64 \
    && chmod +x /usr/local/bin/micromamba \
    && micromamba create -y -p /opt/pbmm2 -c bioconda -c conda-forge "pbmm2=${PBMM2_VERSION}" \
    && printf '#!/bin/bash\nexec /opt/pbmm2/bin/pbmm2 "$@"\n' > /usr/local/bin/pbmm2 \
    && chmod +x /usr/local/bin/pbmm2 \
    && micromamba clean -a -y \
    && pbmm2 --version

ENTRYPOINT ["/bin/bash"]
