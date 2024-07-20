# ChemPlasKin

ChemPlasKin is a free code optimized for zero-dimensional (0D) simulations of neutral gas chemical kinetics coupled with non-equilibrium plasma.

## Overview

ChemPlasKin integrates an electron Boltzmann equation solver, [CppBOLOS](https://github.com/ShaoX96/CppBOLOS), with the open-source combustion library [Cantera](https://cantera.org) at the source code level.
ChemPlasKin enables the calculation of species concentration and gas temperature over time within a unified gas-plasma framework.
This approach enables accurate modeling of both chemical thermal effects and plasma-induced heating, including fast gas heating and slower vibrational-translational relaxation processes.

Check our [paper](https://doi.org/10.1016/j.jaecs.2024.100280) and cite:

Shao, X., Lacoste, D. A., & Im, H. G. (2024). ChemPlasKin: A general-purpose program for unified gas and plasma kinetics simulations. Applications in Energy and Combustion Science, 100280. https://doi.org/10.1016/j.jaecs.2024.100280

## Key Features

- **Unified ODE system**: neutral gas and plasma kinetics are solved together in a single ODE system
- **Versatility**: Suitable for plasma assisted combustion (PAC) and plasma assisted fuel reforming.
- **Compatibility**: Maintains compatibility with Bolsig+ cross-section input format and ZDPlasKin input mechanism.
- **High Performance**: Written in pure C++, at least 3x faster than ZDPlasKin + Cantera/CHEMKIN method.
- **Heat Loss model**: Designed for nanosecond pulsed discharges in pin-pin electrode configurations.

## Getting Started

This section provides details on how to modify the Cantera source code and compile it for the usage of ChemPlasKin. 

1. **Prepare working directories**: (eg.`~/Destop/`)
   
   ```sh
   cd ~/Desktop
   mkdir ChemPlasProject
   cd ChemPlasKinProject
   ```

2. **Clone the Repositories**
   
   ```sh
   git clone --recursive https://github.com/Cantera/cantera.git
   git clone https://github.com/ShaoX96/ChemPlasKin.git
   ```
   
   Now under `ChemPlasKinProject/` you should have `cantera/` and `ChemPlasKin/`.

3. **Test [Compiling Cantera from Source](https://cantera.org/install/compiling-install.html#sec-compiling)**:  
   You should be familiar with [Compiling Cantera from Source](https://cantera.org/install/compiling-install.html#sec-compiling).
   A Conda environment is recommended for [Compilation Requirements](https://cantera.org/compiling/compilation-reqs.html#sec-conda).
   You should be able to compile the original Cantera source before making any modifications to it:

```sh
cd cantera
git checkout 3.0
scons build
```

4. **Obtain external libraries for ChemPlasKin**:
- Create a new branch for Cantera (recommended)

```shell
git checkout -b for_chemplaskin
git branch
```

- Fetch the source code of [CppBOLOS](https://github.com/ShaoX96/CppBOLOS) and [muParser](https://beltoforion.de/en/muparser/) and put them under `cantera/ext/bolos/` and `cantera/ext/muparser/`, respectively.
5. **Update `ext/SConscript`**: 
   
   ```shell
   cp ../ChemPlasKin/ext/Sconscript ext/Sconscript
   ```

6. **Extend Cantera kinetics module**

```shell
cp ../ChemPlaKin/include/kinetics/*.h include/cantera/kinetics/
cp ../ChemPlasKin/include/base/Solution.h include/cantera/base/
cp ../ChemPlaKin/src/kinetics/*.cpp src/kinetics/
cp ../ChemPlasKin/src/base/Solution.cpp cantera/src/base/
```

7. **Stage and Commit Changes (Recommended)**
   
   Make sure you have set your email and name in your Git configuration.
   
   ```shell
   git config --global user.email "your.email@example.com"
   git config --global user.name "Your Name"
   ```
   
   Then commit your changes, for example: 
   
   ```sh
   git add .
   git commit -m "Modify Cantera source for ChemPlasKin usage."
   ```

8. **Compile new library**
   
   ```sh
   cd cantera/
   scons build
   ```
   
   The compiled Cantera library is under `cantera/build/lib`. It should be linked to ChemPlasKin through `ChemPlasKin/CMakeLists.txt`:
   
   ```sh
   link_directories("../cantera/build/lib")
   target_link_libraries(ChemPlasKin cantera_shared ${ACCELERATE_FRAMEWORK} Threads::Threads)
   ```

9. **Build ChemPlasKin**
   
   ```sh
   cd ChemPlasKin
   mkdir build
   cmake ..
   make
   ```

10. **Run an example**:
    Use `-case` to specify the running case path where `controlDict` and `chemPlasProperties` files are located. Use optional `-log` flag to control log level: `NONE`, `WARNING`, `INFO`(default) , or `DEBUG`. 
    
    ```sh
    ./ChemPlasKin -case ../examples/H2O2He -log DEBUG
    ```

11. **Check results**:
    
    ```sh
    cd ../examples/H2O2He
    python plot.py
    ```

## Data input

Two input data files, cross section and reaction mechanism, are needed, as specified in `chemPlasProperties`:

```sh
csDataFile       "<case>/../../data/LXCat/bolsigdb_H2O2HE.dat";
mechFile         "<case>/../../data/PAC_kinetics/ZDPlasKin_kinetics/Mao-H2O2He/H2O2HE_Mao.yaml";
```

The `data/` directory contains:

- `LXCat/`: cross section data in [LXCat](https://nl.lxcat.net/home/) format for CppBOLOS
- `PAC_kinetics/`: unified gas-plasma mechanism files in YAML format

An optional parser tool, `parsePlasKin.py`, is available to convert ZDPlasKin input mechanism files 
into the human-readable YAML format automatically.

```sh
cd data/PAC_kinetics/ZDPlasKin_kinetics/kineticsParser
python parsePlasKin.py --input "plasmaH2O2.inp" --output "parsedPlasKin.yaml"
```

## Acknowledgments

ChemPlasKin uses the Cantera chemical kinetics software, which is developed and maintained by the Cantera Developers. Cantera is an open-source suite of tools for problems involving chemical kinetics, thermodynamics, and transport processes. More information about Cantera can be found at [cantera.org](https://cantera.org).

The [CppBOLOS](https://github.com/ShaoX96/CppBOLOS) solver is build upon [BOLOS](https://github.com/aluque/bolos/tree/master).

ChemPlasKin project is funded by Computational Reacting Flow Laboratory (CRFL) led by Professor [Hong G. Im](https://www.kaust.edu.sa/en/study/faculty/hong-im) at 
King Abdullah University of Science and Technology ([KAUST](https://www.kaust.edu.sa/en/)), Thuwal, Saudi Arabia.

## License

ChemPlasKin is open-sourced under the [LGPLv2 License](https://www.gnu.org/licenses/old-licenses/lgpl-2.0.html).

### Using Cantera

ChemPlasKin integrates Cantera for handling detailed chemistry-plasma kinetics. 
Users are advised that Cantera is distributed under its own license terms.

## Links

- [ChemPlasKin paper](https://doi.org/10.1016/j.jaecs.2024.100280)
- [BOLOS GitHub](https://github.com/aluque/bolos/tree/master)
- [CppBOLOS GitHub](https://github.com/ShaoX96/CppBOLOS)
- Bolsig+ Reference Paper: [G. J. M. Hagelaar and L. C. Pitchford, "Solving the Boltzmann equation to obtain electron transport coefficients and rate coefficients for fluid models", Plasma Sources Science and Technology, 2005](https://iopscience.iop.org/article/10.1088/0963-0252/14/4/011)

## Disclaimer

ChemPlasKin is an independent project that incorporates the Cantera software library. Any issues, bugs, or vulnerabilities found in ChemPlasKin do not necessarily relate to the Cantera software itself.
