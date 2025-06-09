# Define a base to inherit from so ARGs can be collected here
FROM alpine:latest as base
ARG NUM_JOBS=16
ARG HTSLIB_VERSION=1.21
ARG SAMTOOLS_VERSION=1.21

# All builds from source are installed to their own directory so the
# essential file from them can be retrieved separately as needed.

# Build htslib and samtools from source
FROM base as build_htslib
RUN apk update
WORKDIR /build

RUN apk add --no-cache \
    build-base \
    autoconf \
    automake \
    gsl-dev \
    zlib-dev \
    bzip2-dev \
    xz-dev \
    ncurses-dev \
    wget
ENV HTSLIB=htslib-${HTSLIB_VERSION}
ENV SAMTOOLS=samtools-${SAMTOOLS_VERSION}
RUN cd /build \
    && wget -nv https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/${HTSLIB}.tar.bz2 \
    && tar -xf ${HTSLIB}.tar.bz2 \
    && cd ${HTSLIB} \
    && autoreconf -i \
    && ./configure --prefix=$(pwd) \
    && make -j${NUM_JOBS} \
    && make install
RUN cd /build \
    && wget -nv https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/${SAMTOOLS}.tar.bz2 \
    && tar -xf ${SAMTOOLS}.tar.bz2 \
    && cd ${SAMTOOLS} \
    && ./configure --prefix=$(pwd) \
    && make -j${NUM_JOBS} \
    && make install

# Build dnmtools
FROM base as build_dnmtools
RUN apk update
WORKDIR /build
ENV HTSLIB=htslib-${HTSLIB_VERSION}
ENV SAMTOOLS=samtools-${SAMTOOLS_VERSION}
RUN apk add --no-cache \
    build-base \
    autoconf \
    automake \
    gsl-dev \
    zlib-dev \
    bzip2-dev \
    xz-dev \
    wget \
    gzip \
    bash

# Copying the install from within the build tree helps keep things
# smaller than installing earlier and copying the /usr or /usr/local
COPY --from=build_htslib /build/${HTSLIB}/lib /usr/lib
COPY --from=build_htslib /build/${HTSLIB}/include /usr/include
COPY --from=build_htslib /build/${SAMTOOLS}/bin /usr/bin

RUN mkdir /build/dnmtools
COPY . /build/dnmtools
RUN cd /build/dnmtools \
    && ./autogen.sh \
    && ./configure --prefix=$(pwd) \
    && make -j${NUM_JOBS} \
    && make -j${NUM_JOBS} check \
    && make -j${NUM_JOBS} distcheck \
    && make install

# Build a light-weight image just with binaries
FROM base
ENV HTSLIB=htslib-${HTSLIB_VERSION}
RUN apk update
WORKDIR /build
RUN apk add \
    gsl-dev \
    zlib-dev \
    bzip2-dev \
    xz-dev \
    libstdc++ \
    libgomp
COPY --from=build_htslib /build/${HTSLIB}/lib /usr/lib
COPY --from=build_dnmtools /build/dnmtools/dnmtools /usr/bin

ENTRYPOINT ["dnmtools"]
