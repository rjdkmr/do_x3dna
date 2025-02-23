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

RUN mkdir gmx_installed
WORKDIR /build/external/gromacs-2025.0


ENV GMX_SRC=/build/external/gromacs-2025.0
ENV GMX_PATH=/build/external/gmx_installed
RUN bash /scripts/build_gromacs.sh

WORKDIR /build
RUN cmake -DGMX_PATH=${GMX_PATH} -DGMX_SRC=${GMX_SRC} .. && make && sudo make install

CMD ["do_x3dna", "-h"]
