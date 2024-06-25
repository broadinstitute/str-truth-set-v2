FROM python:3.9-slim-bullseye

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


ENV SAMTOOLS_VERSION="1.20"
RUN wget https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar xjf htslib-${SAMTOOLS_VERSION}.tar.bz2 \
    && rm htslib-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd htslib-${SAMTOOLS_VERSION} \
    && ./configure \
    && make \
    && make install \
    && make clean

# install gcloud
ENV GCLOUD_SDK_VERISON="481.0.0"
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && tar xzf google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && /google-cloud-sdk/install.sh -q \
    && /google-cloud-sdk/bin/gcloud components update --quiet \
    && rm google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz

ENV PATH=/google-cloud-sdk/bin:$PATH

# install bedtools
ENV BED_TOOLS_VERSION="2.31.0"
RUN wget https://github.com/arq5x/bedtools2/releases/download/v${BED_TOOLS_VERSION}/bedtools.static \
	&& chmod 777 bedtools.static \
	&& mv bedtools.static /usr/bin/bedtools

# install dipcall
RUN wget https://github.com/lh3/dipcall/releases/download/v0.3/dipcall-0.3_x64-linux.tar.bz2 \
 	&& tar -jxf dipcall-0.3_x64-linux.tar.bz2 \
	&& rm dipcall-0.3_x64-linux.tar.bz2


ENTRYPOINT ["/bin/bash"]