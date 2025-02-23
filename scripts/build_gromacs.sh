cd $GMX_SRC

if [ -d build ]; then
    rm -rf build
fi

mkdir build && cd build

cmake .. -DGMX_SIMD=SSE2 -DGMX_GPU=off -DGMXAPI=OFF -DGMX_INSTALL_LEGACY_API=on -DGMX_FFT_LIBRARY=fftpack \
        -DCMAKE_INSTALL_PREFIX=/home/raj/software/gromacs/gmx_dynamic

make
make install

cd ..