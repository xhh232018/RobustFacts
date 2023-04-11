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
The Wikidata dump used in our paper can be downloaded [here](https://drive.google.com/drive/folders/1R6rH2GBbD85PTCHXL2qFISNl74y45FwT?usp=sharing).
As this is a huge dump, it takes several minutes to construct the data graph. Please put the files into "test" folder.

## Execute
After compiling the source code, you can find the binary file 'ptbClean' under the 'build/bin' directory. 
Execute the binary with the following command ./ptbClean -e edges_file -n nodetype_file
-f facts_file -t perturbation_type -s sampling_or_not -r sampling_ratio, -p factPrefix, -o output_file.

NOTE: sampling_or_not is a bool variable. Please input T or F.

Example (Data perturbation for 3-hop facts of Art celebrities not using samping method). 
```zsh
cd build/bin
./ptbClean -e ../../test/ -n ../../test/ -f ../../test/test.txt -t data -s F -r 0.05 -p art3 -o ./art.txt
```
