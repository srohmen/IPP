#!/bin/bash

mkdir build

# use the fixup branch
git clone -b fixup https://github.com/srohmen/phreeqrm build/phreeqcrm
cmake -S build/phreeqcrm -B build/phreeqcrm-build -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DPHREEQCRM_BUILD_MPI=ON -DPHREEQCRM_USE_ZLIB=ON || die
cmake --build build/phreeqcrm-build || die
cmake --install build/phreeqcrm-build --prefix build/phreeqcrm-install || die

git clone https://github.com/srohmen/palabos build/palabos
cmake -S build/palabos -B build/palabos-build -G Ninja -DCMAKE_BUILD_TYPE=Release -DENABLE_MPI=ON -DENABLE_SMP_PARALLEL=ON || die
cmake --build build/palabos-build || die
cmake --install build/palabos-build --prefix build/palabos-install || die

cmake -S . -B build/ipp-build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPALABOS_ROOT_DIR=build/palabos-install -DPHREEQCRM_ROOT_DIR=build/phreeqcrm-install || die
cmake --build build/ipp-build || die

