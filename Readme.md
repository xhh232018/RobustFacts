# Gauging the Robustness of Facts

## Requirement
Cmake >= 3.0 C++ Boost Library >= 1.63

## Compile
Under the root directory of the project, execute the following commands to compile the source code.

```zsh
mkdir build
cd build
cmake ..
make
```

## Dataset
The Wikidata dump and OFs used in our paper can be downloaded [here](https://drive.google.com/drive/folders/1R6rH2GBbD85PTCHXL2qFISNl74y45FwT?usp=sharing).
As this is a huge dump, it takes several minutes to load the data graph. 
Please put all files into "test" folder.

## Execute
After compiling the source code, you can find the binary file 'GaugeRobFacts' under the 'build/bin' directory. 
Execute the binary with the following command ./GaugeRobFacts -d dump_file -f facts_file -t perturbation_type(entity/data) -s sampling_or_not(T/F) -r sampling_ratio(0 < r < 1), -p factPrefix, -o output_file_name -l time_limit(s).

NOTE: sampling_or_not is a bool variable. Please input T or F.

Example (Entity perturbation for 3-hop facts of Art celebrities not using samping method with 5-min time limit). 
```zsh
cd build/bin
./ptbClean -d ../../test/graph60M.dat  -f ../../test/test.txt -t entity -s F -r 0.05 -p art3 -o ./art.txt -l 300
```
