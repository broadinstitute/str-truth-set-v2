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
    time


# install gcloud
ENV GCLOUD_SDK_VERISON="481.0.0"
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && tar xzf google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && /google-cloud-sdk/install.sh -q \
    && /google-cloud-sdk/bin/gcloud components update --quiet

ENV PATH=/google-cloud-sdk/bin:$PATH


# install htslib  --  install in downstream images with:
#     COPY --from=weisburd/base-bam-tools:latest  /usr/local/bin/tabix /usr/local/bin/tabix
#     COPY --from=weisburd/base-bam-tools:latest  /usr/local/bin/htsfile /usr/local/bin/htsfile
#     COPY --from=weisburd/base-bam-tools:latest  /usr/local/bin/bgzip /usr/local/bin/bgzip
#     COPY --from=weisburd/base-bam-tools:latest  /usr/local/bin/samtools /usr/local/bin/samtools
ENV SAMTOOLS_VERSION="1.20"
RUN wget https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd htslib-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& make clean


# install samtools  --  install in downstream images with:
#     COPY --from=weisburd/base-bam-tools:latest /usr/local/bin/samtools /usr/local/bin/samtools
RUN wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd samtools-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& make clean

RUN python3 -m pip install --no-cache-dir --upgrade cmake pip


# install bwa  --  install in downstream images with:
#     COPY --from=weisburd/base-bam-tools:latest  /usr/local/bin/bwa /usr/local/bin/bwa
RUN wget https://github.com/lh3/bwa/archive/master.zip -O bwa.zip \
	&& unzip bwa.zip \
	&& cd bwa-master \
	&& make \
	&& cp bwa /usr/local/bin/ \
	&& cd .. \
	&& rm -rf bwa*
