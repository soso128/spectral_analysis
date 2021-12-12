# spectral_analysis
Spectral fit code for the search for the [Diffuse Supernova Neutrino Background at Super-Kamiokande I-IV](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.104.122002)

## Installation instructions
### Requirements

This code has successfully worked with the following software:
* Python >= 3.6.5
* Pybind11
* Numpy 1.14.3
* Scipy 1.1.0

Using other versions, especially older ones, might work but this is not guaranteed. The code is definitely not compatible with Python 2.

### Installation
You need to compile the snrate module (C++ code interfaced using Pybind11). In the Makefile, change the INC path to your Python path. Then compile using
    make

The code is now ready to use. snrate can be used like a regular Python module.

In order to get information about the detector effects, simulated events are used for the different phases of Super-Kamiokande (SK). These events are stored in [this folder](), which is too large to be stored on Github. So you should do:
```
wget ...
cp npy_files.tar.gz mc/
cd mc/
tar -xzf npy_files.tar.gz
```

## Running the fits

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

Figures 26-28 from [the published paper](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.104.122002) can be reproduced by running this code for the "horiuchi" model, using the fit result folder from the data release.

## Making plots

The main code is also fitlike.py, but this time you have to use the following options:

```
<python executable>  fitlike.py <model name>  <input folder for model>  --drawonly 1
  ```
## Adding a new model

You just need to add a neutrino flux file in the models/ folder, named as flux_cross_\<model name\>.dat .  The file need to have at least two columns: the first column is the neutrino energy and the second column is the flux in nu/cm2/s.

**Warning**: The current code computed the flux normalization using the scipy "quad" function over the energy range considered for the DSNB spectral analysis. Moreover, the **reconstructed energy** distribution of the DSNB rate is obtained using a histogram with 1 MeV-wide bins. So if your neutrino flux is sharply peaked around a given energy the code is not likely to work as it is. 
