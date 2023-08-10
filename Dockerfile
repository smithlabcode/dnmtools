

# Build htslib from source
FROM alpine:latest as build_sam_hts
RUN apk update 
WORKDIR /app

RUN apk add build-base 
RUN apk add autoconf 
RUN apk add automake 
RUN apk add gsl-dev zlib-dev bzip2-dev xz-dev ncurses-dev
RUN apk add wget

ENV HTS=1.18
ENV HTSDIR=htslib-${HTS}

RUN mkdir /app/htslib &&  cd /app/htslib \
  && wget https://github.com/samtools/htslib/releases/download/${HTS}/${HTSDIR}.tar.bz2 \
  && tar -xf ${HTSDIR}.tar.bz2 \
  && cd /app/htslib/${HTSDIR} \
  &&  autoreconf -i &&  ./configure --prefix=`pwd` && make && make install 

ENV SAM=1.18
ENV SAMDIR=samtools-${SAM}

RUN mkdir /app/samtools &&  cd /app/samtools \
  && wget https://github.com/samtools/samtools/releases/download/${SAM}/${SAMDIR}.tar.bz2 \
  && tar -xf ${SAMDIR}.tar.bz2 \
  && cd /app/samtools/${SAMDIR} \
  &&  ./configure --prefix=`pwd` && make && make install 

# Build dnmtools
FROM alpine:latest as build_dnmtools
RUN apk update 
WORKDIR /app

ENV HTS=1.18
ENV HTSDIR=htslib-${HTS}

ENV SAM=1.18
ENV SAMDIR=samtools-${SAM}

RUN apk add build-base 
RUN apk add bash
RUN apk add autoconf 
RUN apk add automake 
RUN apk add gsl-dev zlib-dev bzip2-dev xz-dev wget gzip

COPY --from=build_sam_hts /app/htslib/${HTSDIR}/lib/. /usr/lib
COPY --from=build_sam_hts /app/htslib/${HTSDIR}/include/. /usr/include
COPY --from=build_sam_hts /app/samtools/${SAMDIR}/bin/. /usr/bin
RUN mkdir /app/dnmtools
COPY . /app/dnmtools
RUN cd /app/dnmtools && ./autogen.sh && mkdir build && cd build \ 
 && ../configure --prefix=`pwd` && make && make check && make distcheck && make install


# Build a light-weight image just with binaries
FROM alpine:latest

ENV HTS=1.18
ENV HTSDIR=htslib-${HTS}

RUN apk update 
WORKDIR /app
RUN apk add gsl-dev zlib-dev bzip2-dev xz-dev libstdc++ libgomp

COPY --from=build_sam_hts /app/htslib/${HTSDIR}/lib/. /usr/lib
COPY --from=build_sam_hts /app/htslib/${HTSDIR}/include/. /usr/include

COPY --from=build_dnmtools /app/dnmtools/build/dnmtools /usr/bin

ENTRYPOINT ["dnmtools"]
