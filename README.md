# TSUNER -- Tsunami with Non-linear Erosion and Runup

This repository host the public version for TSUNER (**Tsu**nami with **N**on-linear **E**rosion and **R**unup). It is a geophysical model for numerical simulation of sediment evolution, mainly designed for tsunami waves, although it could be used for other geophysical events.

## Capabilities

binary map
change sediments for parametrization

## Requirements

The model is desinged to run in a NVIDIA GPU, so you need to have CUDA installed in your machine. It also uses [NetCDF](https://www.unidata.ucar.edu/software/netcdf) for IO. You would also need a standard C++ compiler, such as GNU compiler. You will also need [GNU Make](https://www.gnu.org/software/make/).


## Installation

Make sure you have all the required libraries from the previous step installed. The program include two different Makefiles that you need to modify according to your installation paths, one located in the root folder and the other in the GPU folder. In particular, make sure to set CUDA_PATH to your CUDA installation (by default in /usr/local/cuda) and your NetCDF path in BOTH makefiles. This should be easy in Linux and Mac. In Windows, you probably want to do it under a Linux subsystem installation.

After compilation, run the program with:
```
./tsuner.exe data.dat
```
To clean the installation, just perform *make clean*.

## Usage

### Controlling the simulation
The simulation are mainly controlled by two files: *Sim_parameters.cxx* and *data.dat*. In *Sim_parameters* you can controll:
- The order of the numerical simulation (1st or 2nd)
- The boundary conditions
- The Units of the initial condition files (by default they are in cm)
- It also allows to switch the behaviour of the program to read files for input or to define a set of initial conditions from the *cini.cxx* file.
- To remove all sediment evolution behaviour, turning the code into standard Shallow-water equations.
- Other general parameters you likely do not want to change, such as gravity or volumetric mean earth radius.
Note that you need to compile the code again after any changes in this file. There is no check for inconsistent choices (like periodic boundary conditions and read external files) so make sure your selection is correct.

The *data.dat* file controlls simulation parameters such as total simulation time, input and output file names, etc. The values are described in the file itself. A few things to highlight from this file:
- It includes values like xmin, xmax, ymin and ymax for domain definition. They are used if the program is set to NOT read external files as input in *Sim_parameters*, otherwise they are ignored.
- Sediment parameters characteristics are defined here.
- You can also controll the number of discretization points in *x* and *y* direction. Note that the program will interpolate the bathymetry and boundary conditions to fit your desired discretization. It is useful to run low-resolution model at first for fast results, and then scale-up to desired resolution. It does not make much sense to use more discretization points than the original resolution of your input files, but it is certainly possible to do so.
See next section for important remarks about input files. Note that if you modify *data.dat* you do NOT need to compile the program again.


### Input files

The program is desing to be compliant with NOAA Center for Tsunami Research (NCTR) input files. Failure to fully comply with these instructions for the input files are the most common case of the code not working properly. Please, make sure to follow these instructions carefully.

The model takes as input a bathymetry file in data.dat. It is a file in spheric coordinates with coordinates named (lon,lat) and NetCDF format. See the example file for reference. The file will be transformed to cartesian coordinates by the program. IMPORTANT: To be compliant with NOAA NCTR input files, the longitude values are centered in the trasnported to avoid a sign change in the Pacific Ocean. You need to add +360 to your longitude values express in standard coordinates centered in Greenwich. Follow the same order for variable definition as the variables in the example file: from left to right and from bottom to top.

Other important input files are the tsunami sources. They are three different netcdf files for water free surface and velocties. You specify the value of the free surface file in *data.dat* and the code assumes the file name for the velocities is exactly the same except they finish in *u* and *v*. Example:
```
my_source_file_h.nc
my_source_file_u.nc
my_source_file_v.nc

```
These files must be in the same folder and also in spheric coordinates. They may use as coordinates (lat, lon), (LAT, LON) or (x,y), plus name of variables *HA*, *UA* and *VA*. They are expected to be in centimeters, as this is the format used in NCTR. However, if your input files are in meters or any other system, it is not neccesary to change the input files. Instead, you can change the parameter UNITS in *Sim_parameters.hxx* file. If you have meters, you set UNITS to 1. If you have centimeters (by default) you set UNITS to 0.01 and so on. See example files.




### Output files

The program will create a NetCDF file with the name specified in *data.dat*. It saves the bathymetry evolution, free surface, velocities, the maximum water depth accross all time steps and the minimum bathymetry value across all time steps. Any NetCDF reader can be used but I recommend [Paraview](https://www.paraview.org/) or Python.

## Contributors

This program has been mainly developed by Ernesto Guerrero Fernandez, under the University of Washington and the National Oceanic and Atmospheric Administration (NOAA) and the NOAA Center for Tsunami Research (NCTR). Other relevant contributors include Manuel Jesus Castro Diaz, Marc de la Asuncion (both from University of Malaga) and Yong Wei (from NOAA/University of Washington).

## Scientific background

The model is solves the PDE system for the Shallow-Water equations with a coupled transport equation for the bathymetry function. The numerical scheme is based on the path-conservative technique. You can see the full paper describing the model in the following link: [https://repository.library.noaa.gov/view/noaa/67438/noaa_67438_DS1.pdf](https://repository.library.noaa.gov/view/noaa/67438/noaa_67438_DS1.pdf)

## Contact

For questions or collaboration you can find me at ernestog@uw.edu. You can report bugs or request new features through this Github page.

## Cite

You can cite this code with the following reference:
Fernández, E. G., Díaz, M. C., Wei, Y., & Moore, C. (2024). Modeling sediment movement in the shallow-water framework: A morpho-hydrodynamic approach with numerical simulations and experimental validation. Ocean Modelling, 192, 102445.
Or visit [this web](https://www.sciencedirect.com/science/article/pii/S146350032400132X?via%3Dihub) to download the reference.

## License
This project is licensed under the **GNU Affero General Public License v3.0** - see the [LICENSE](LICENSE) file for details.
