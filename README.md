Software accompanying the paper ["Fault Attacks on CCA-secure Lattice KEMs" (PDF)](https://eprint.iacr.org/2021/064.pdf). To appear at CHES 2021.

Authors:
  - [Peter Pessl](https://pessl.cc): <peter@pessl.cc>
  - [Lukas Prokop](https://www.iaik.tugraz.at/person/lukas-prokop/): <lukas.prokop@iaik.tugraz.at>

## Content
Matlab and C sources used for simulations, includes attacks for Kyber (v2, all parameter sets) and NewHope512. The scripts are currently set up to execute an attack on Kyber512 using 7500 faults and NewHope512 using 15000 faults, but can be easily configured for other scenarios.

## Requirements

* Matlab with the Parallel Processing Toolbox (simulations on R2018b)
* at least 64GB RAM

## How to run
For attacks on Kyber, run the following (same for NewHope):
* go to c/kyber_fault and run make
* copy the binary files "kyberXXX_fault" to the matlab directory
* run the matlab script kyber_attack.m 

## Important configuration parameters
The first section of the Matlab scripts contains several configuration options. Some of the most important ones are:
* `kyberk` selects the used kyber parameter set. Valid options are 2, 3, and 4. 
* `faultrange` selects the number of faults for the sweep. Ranges used for the reported simulations can be found in the comments.
* `maxiterations` sets the maximum number of iterations for the statistical solving method.
* `repeats` determines how many experiments to run per fault quantity (for determining the success rate)