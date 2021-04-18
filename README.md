Software accompanying the paper ["Fault Attacks on CCA-secure Lattice KEMs" (PDF)](https://eprint.iacr.org/2021/064.pdf). Paper published in [TCHES 2021, Issue 2](https://tches.iacr.org/index.php/TCHES/article/view/8787).

Authors:
  - [Peter Pessl](https://pessl.cc): <peter@pessl.cc>
  - [Lukas Prokop](https://www.iaik.tugraz.at/person/lukas-prokop/): <lukas.prokop@iaik.tugraz.at>

## Content
Matlab and C sources used for simulations, includes attacks for Kyber (v2, all parameter sets) and NewHope512. The scripts are currently set up to execute a single attack on Kyber512 using 7500 faults and NewHope512 using 15000 faults, but can be easily configured for other scenarios.

### Folders
* `c` contains programs used to generate linear inequalities. The programs are based on the reference implementations of Kyber and NewHope; simulated faults are injected at the positions described in the paper. The program then generates the resulting inequalities and dumps them into files.
* `matlab` constains the scripts used for solving the linear inequalities and recovering the key. First calls the appropriate program to generate the inequalities, then parses the files, and then runs the solving approach described in the paper.


## Requirements

* Matlab with the Parallel Processing Toolbox (tested with R2018b and R2020b)
* at least 64GB RAM

## How to run
For attacks on Kyber, run the following (same for NewHope):
* open the `matlab` folder
* run the matlab script kyber_attack.m (with the `matlab` folder as current directory) 

## Important configuration parameters
The first section of the Matlab scripts contains several configuration options. Some of the most important ones are:
* `kyberk` selects the used kyber parameter set. Valid options are 2 (Kyber512), 3 (Kyber768), and 4 (Kyber1024). 
* `faultrange` selects the number of faults for the sweep. Ranges used for the reported simulations can be found in the comments. Currently set up to run tests with a fixed number of faults (7500 for Kyber and 15000 for NewHope). 
* `maxiterations` sets the maximum number of iterations for the statistical solving method.
* `repeats` determines how many experiments to run per fault quantity (for determining the success rate). Currently set to 1.
* `plots` allows to turn on analysis plots, which show the current status of key recovery (current rank of the correct key for each coefficient, current probability of the correct key coefficient).
* `store` switches on the storage of results in the file `result.mat` (needed when multiple simulations are performed, i.e., in a sweep)