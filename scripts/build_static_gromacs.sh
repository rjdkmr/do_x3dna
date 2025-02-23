cd gromacs-2025.0

if [ -d build ]; then
    rm -rf build
fi

mkdir build && cd build

cmake .. -DGMX_SIMD=SSE2 -DGMX_GPU=off -DGMXAPI=OFF -DGMX_INSTALL_LEGACY_API=on -DGMX_FFT_LIBRARY=fftpack \
        -DCMAKE_INSTALL_PREFIX=/home/raj/software/gromacs/gmx_static \
        -DGMX_PREFER_STATIC_LIBS=ON -DGMX_OPENMP=OFF -DBUILD_SHARED_LIBS=OFF

make
make install

cd ..