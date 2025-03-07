FROM ubuntu:24.04

RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive TZ=Europe/London apt-get -y install tzdata
RUN apt-get install -y cmake gcc g++ curl libssl-dev build-essential pkg-config

COPY scripts /scripts
COPY cmake /cmake
COPY src /src
COPY CMakeLists.txt /CMakeLists.txt

RUN mkdir /build
RUN mkdir /build/external
WORKDIR /build/external
RUN curl -L -O https://ftp.gromacs.org/gromacs/gromacs-2025.0.tar.gz && tar -zxf gromacs-2025.0.tar.gz

WORKDIR /build/external/gromacs-2025.0
ENV GMX_SRC=/build/external/gromacs-2025.0
ENV GMX_PATH=/build/external/gmx_installed
RUN mkdir /build/external/gmx_installed
WORKDIR /build/external/gromacs-2025.0
RUN mkdir build
WORKDIR /build/external/gromacs-2025.0/build
RUN cmake .. -DGMX_SIMD=SSE2 -DGMX_GPU=off -DGMXAPI=OFF -DGMX_INSTALL_LEGACY_API=on -DGMX_FFT_LIBRARY=fftpack -DCMAKE_INSTALL_PREFIX=${GMX_PATH}
RUN make -j4
RUN make install

WORKDIR /build
RUN cmake -DGMX_PATH=${GMX_PATH} -DGMX_SRC=${GMX_SRC} .. && make && make install

CMD ["do_x3dna", "-h"]
