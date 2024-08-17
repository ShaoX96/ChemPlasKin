# N2 plasma kinetics simulation
All input parameters are hard coded. `plasmaReactor.h` has been
simplified. Using external profiles of electron density and electric field. This example shows how to use ChemPlasKin purely for
plasma kinetics simulations.

**Example Source**: http://www.zdplaskin.laplace.univ-tlse.fr/external-profiles-of-electron-density-and-electric-field/index.html

To run it:
```shell
mkdir build
cd build
cmake ..
make
./ChemPlasKin
```
You don't need to specify `-case` since it doesn't read user parameter input.