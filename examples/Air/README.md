# Nanosecond spark discharges in air
All input parameters are hard coded. Original `plasmaReactor.h` from `ChemPlasKin/src/` is
used. External profiles of electron density and electric field are imposed.

**Example Source**:

[1] https://iopscience.iop.org/article/10.1088/0022-3727/46/46/464010

[2] https://doi.org/10.1016/j.combustflame.2022.111990

To run the case:

```shell
mkdir build
cd build
cmake ..
make
./ChemPlasKin
```
You don't need to specify `-case` since it doesn't read user parameter input.