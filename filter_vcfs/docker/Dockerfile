FROM python:3.10-slim-bullseye

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

ENV SAMTOOLS_VERSION="1.20"
RUN wget -q https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
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

RUN apt-get install -y --no-install-recommends bc dvipng texlive-latex-extra texlive-fonts-recommended cm-super

# set up "aliases" for MacOSX commands
RUN cp /bin/zcat /bin/gzcat

# install str-analysis
RUN python3 -m pip install --no-cache-dir --upgrade cmake pip

# pinned to a specific str-analysis commit for reproducible builds (previously tracked main HEAD via a randbyte cache-buster)
RUN python3 -m pip install --no-cache-dir git+https://github.com/broadinstitute/str-analysis@2d325969b3e53148dd1d66459d9ce51e19019642


# install str-truth-set repo
COPY str-truth-set/requirements.txt /
# step-pipeline is only used by the pipeline drivers that run locally (not by the analysis/plot scripts baked into
# this image) and has no PyPI release installable on this Python, so exclude it from the image install
RUN grep -vE '^step-pipeline([<>=!~ ]|$)' /requirements.txt > /requirements.image.txt \
    && python3 -m pip install --no-cache-dir -r /requirements.image.txt

COPY str-truth-set /str-truth-set
COPY filter_vcfs/scripts /filter_vcfs/scripts


ENTRYPOINT ["/bin/bash"]