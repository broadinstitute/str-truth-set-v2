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

ENTRYPOINT ["/bin/bash"]