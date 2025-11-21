# IPP
Interface Palabos-PhreeqC is a coupling code between the Lattice-Boltzmann transport solver Palabos and the chemical speciation solver PhreeqC.

There are system dependencies required to build IPP:
- OpenMPI (or any other standard conforming MPI implementation)
- VTK
- Boost
- GSL
- JsonCpp
- ZLIB

Install these using your package manger of your Linux distribution. Moreover, IPP requires modified/patched versions of both, [Palabos](https://github.com/srohmen/palabos) (based on version 2.0), and [PhreeqcRM](https://github.com/srohmen/phreeqrm/tree/fixup) (based on version 3.4).

To download the Palabos and PhreeqcRM dependencies and build the project, there is a build script located in this repository:
```
./build.sh
```

To run the code, there is a test application called `IPPTest`. There are a couple of examples as input for IPPTest available, too. After the build script has built the library code and the `IPPTest` example, you can run them using mpirun:
```
# switch to directory
cd build/ipp-build/examples/IPPTest
# run the example using 16 CPU cores
mpirun -np 16 ./IPPTest -i input/calcite_nucleation_high_SI.json
```

To cite this code, and for further information, please use this PhD thesis:

Rohmen, Stephan (2024). "Pore-scale reactive transport modeling in cementitious materials: development and application of a high-performance computing code based on the Lattice-Boltzmann method", RWTH Aachen University. DOI: [10.18154/RWTH-2025-03516](https://doi.org/10.18154/RWTH-2025-03516)
