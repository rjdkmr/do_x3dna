FROM ubuntu:24.04

RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive TZ=Europe/London apt-get -y install tzdata
RUN apt-get install -y cmake gcc g++ curl libssl-dev build-essential pkg-config

RUN mkdir /build
RUN mkdir /build/external
WORKDIR /build/external

RUN curl -L -O https://ftp.gromacs.org/gromacs/gromacs-2025.0.tar.gz && tar -zxf gromacs-2025.0.tar.gz

RUN mkdir gmx_installed
WORKDIR /build/external/gromacs-2025.0
RUN mkdir build
WORKDIR /build/external/gromacs-2025.0/build
RUN cmake .. -DGMX_SIMD=SSE2 -DGMX_GPU=off -DGMXAPI=OFF -DGMX_INSTALL_LEGACY_API=on -DGMX_FFT_LIBRARY=fftpack \
    -DCMAKE_INSTALL_PREFIX=/build/external/gmx_installed && \
    make && \
    make install

ENV GMX_PATH=/build/external/gmx_installed
ENV GMX_SRC=/build/external/gromacs-2025.0

CMD ["echo", "testing"]
