# spectral_analysis
Spectral fit code for supernova relic neutrinos

## Installation instructions
### Requirements

This code requires Python >= 3.7, with NumPy, SciPy, and Pybind11

### Installation
You need to compile the snrate module (C++ code interfaced using Pybind11). In the Makefile, change the INC path to your Python path. Then compile using
    make

The code is now ready to use. snrate can be used like a regular Python module.


## Running the fit

The main code is fitlike.py . To run the spectral fit for all phases of SK use

```
<python executable>  fitlike.py <model name>  <output folder>  [--sys <systematics flag>] [--thr <analysis threshold for events with != 1 neutron>] 
  [--thr1n <analysis threshold for events with 1 neutron>] [--spall <spallation flag>]
  ```
  
The SN neutrino flux files for the available models can be found in the models/ folders. The flux files have the following naming convention:
  
```
flux_cross_<model name>.dat
```
  
The output folder should be created by the user before running the code.
  
The systematics flags stand for the following options:
* -1: no systematics, fit parameter initialization (should be used only inside the code, not as a command line option)
* 0 : no background systematics
* 1: systematics used in the paper (atmospheric spectral shapes)
* 2: only systematics on the mupi background
* 3: energy scale and resolution systematics only

The analysis threshold of the current analysis is 15.5 MeV in reconstructed energy, independently of the number of tagged neutrons. For events with one neutrons we 
will probably lower this threshold in the future.

The spallation flag is 1 if we account for spallation background for events with != neutron, and 0 otherwise. We assume that spallation has been entirely eliminated for events with 1 tagged neutron.

## Making plots

The main code is also fitlike.py, but this time you have to use the following options:

```
<python executable>  fitlike.py <model name>  <input folder for model>  --drawonly 1
  ```
