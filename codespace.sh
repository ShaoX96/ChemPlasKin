conda env create -f environment.yaml

cd ..

git clone --recursive https://github.com/Cantera/cantera.git

cd cantera
git checkout 3.0
git checkout -b for_chemplaskin

cd ext
git clone https://github.com/ShaoX96/CppBOLOS
mv CppBOLOS bolos
git clone https://github.com/beltoforion/muparser

cd ..
cp ../ChemPlasKin/ext/SConscript ext/SConscript

cp ../ChemPlasKin/include/kinetics/*.h include/cantera/kinetics/
cp ../ChemPlasKin/include/base/Solution.h include/cantera/base/
cp ../ChemPlasKin/src/kinetics/*.cpp src/kinetics/
cp ../ChemPlasKin/src/base/Solution.cpp src/base/

conda run -n ct-build scons build

cd ../ChemPlasKin
mkdir build
cd build
cmake ..
make


./ChemPlasKin -case ../examples/H2O2He -log DEBUG

pip install numpy matplotlib pandas
python plot.py