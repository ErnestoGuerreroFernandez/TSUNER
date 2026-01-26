# TSUNER -- Tsunami with Non-linear Erosion and Runup

This repository hosts the public version of TSUNER (**Tsu**nami with **N**on-linear **E**rosion and **R**unup). It is a geophysical model for the numerical simulation of sediment evolution, mainly designed for tsunami waves, although it could be used for other geophysical events.

## Capabilities

The model allows you to:
- Simulate a tsunami wave without any sediment evolution, from boundary input sources.
- Same as above, but including sediment morphology evolution, with customizable sediment parameters.
- Include a binary map that defines areas where there is no sediment evolution or where it occurs. This is useful to define human-made structures for the tsunami, such as harbors or ports.
- Include Sea Level Rise in the simulations.
- Choose between first and second order in space and time for the total global of accuracy (second order by default).
- Allow for a choice between reading boundary conditions or defining a synthetic, mathematically defined initial condition.
- Related with the previous item, you can run 1D simulations by setting the number of discretization points in one dimension to 1.


Other Capabilities, not yet finished:
- Include tide values in tsunami simulations.
- Define different areas with different sediments within the computational domain.

## Requirements

The model is designed to run in a NVIDIA GPU, so you need to have CUDA installed in your machine. It also uses [NetCDF](https://www.unidata.ucar.edu/software/netcdf) for IO. You would also need a standard C++ compiler, such as the GNU compiler. You will also need [GNU Make](https://www.gnu.org/software/make/).

## Installation

Make sure you have all the required libraries from the previous step installed. The program includes two different Makefiles that you need to modify according to your installation paths, one located in the root folder and the other in the GPU folder. In particular, make sure to set CUDA_PATH to your CUDA installation dir (by default in /usr/local/cuda) and your NetCDF path in both makefiles. This should be easy in Linux and Mac. In Windows, you probably want to do it under a Linux subsystem installation.

After compilation, run the program with:
```
./tsuner.exe data.dat
```
To clean the installation, just perform *make clean*.

## Usage

### Controlling the simulation

The simulations are mainly controlled by two files: *Sim\_parameters.cxx* and *data.dat*. In *Sim\_parameters* you can control:
- The order of the numerical simulation (first or second).
- The boundary conditions.
- The units of the initial condition files (by default they are in cm).
- It also allows you to switch the program's behavior to read files for input or to define a set of initial conditions from the *cini.cxx* file.
- To remove all sediment evolution behaviour, thus turning the code into standard Shallow-water equations.
- Other general parameters you likely do not want to change, such as gravity or the volumetric mean earth radius.
Note that you need to compile the code again after any changes in this file. There is no check for inconsistent choices (like periodic boundary conditions and reading external files) so make sure your selection is correct.

The *data.dat* file controls simulation parameters such as total simulation time, input and output file names, etc. The values are described in the file itself. A few things to highlight from this file:
- It includes values like xmin, xmax, ymin and ymax for domain definition. They are used if the program is set to NOT read external files as input in *Sim\_parameters*, otherwise they are ignored.
- Sediment parameters are defined here.
- You can also control the number of discretization points in the *x* and *y* directions. Note that the program will interpolate the bathymetry and boundary conditions to fit your desired discretization. It is useful to run a low-resolution model at first for fast results, and then scale up to the desired resolution. It does not make much sense to use more discretization points than the original resolution of your input files, but it is certainly possible to do so.
See the next section for important remarks about input files. Note that if you modify *data.dat* you do NOT need to compile the program again.


### Input files

The program is designed to be compliant with NOAA Center for Tsunami Research (NCTR) input files. Failure to fully comply with these instructions for the input files is the most common cause of the code not working properly. Please, make sure to follow these instructions carefully.

The model takes as input a bathymetry file specified in *data.dat*. It is a file in spherical coordinates with coordinates named *(lon, lat)* and NetCDF format. See the example file for reference. The file will be transformed to Cartesian coordinates by the program. IMPORTANT: To be compliant with NOAA NCTR input files, the longitude values are centered in the transported system to avoid a sign change in the Pacific Ocean. You need to add +360 to your longitude values expressed in standard coordinates centered in Greenwich. Follow the same order for variable definition as the variables in the example file: from left to right and from bottom to top.

Other important input files are the tsunami sources. They are three different NetCDF files for the water free surface and velocities. You specify the value of the free surface file in *data.dat* and the code assumes the file name for the velocities is exactly the same except they finish in *_u* and *_v*. Example:

```
my_source_file_h.nc
my_source_file_u.nc
my_source_file_v.nc

```

These files must be in the same folder and also in spherical coordinates. They may use coordinates *(lat, lon)*, *(LAT, LON)* or *(x, y)*, plus the variable names *HA*, *UA* and *VA* (either capital or lowercase). They are expected to be in centimeters, as this is the format used in NCTR. However, if your input files are in meters or any other system, it is not necessary to change the input files. Instead, you can change the parameter UNITS in the *Sim\_parameters.hxx* file. If you have meters, you set UNITS to 1. If you have centimeters (by default) you set UNITS to 0.01 and so on. 


In the example folder you can find a bathymetry and source files as example. Due to Github file size restrictions, the files are low-resolution and do not include many time steps. However, they are enough to give you an idea of the input file format you need to use. You can use a netcdf reader, such as *ncdump*, to visualize variables names and formats.


### Output files

The program will create a NetCDF file with the name specified in *data.dat*. It saves the bathymetry evolution, free surface, velocities, the maximum water depth across all time steps and the minimum bathymetry value across all time steps. Any NetCDF reader can be used but I recommend [Paraview](https://www.paraview.org/) or Python.

## Contributors

This program has been mainly developed by Ernesto Guerrero Fernandez, under the University of Washington and the National Oceanic and Atmospheric Administration (NOAA) and the NOAA Center for Tsunami Research (NCTR). Other relevant contributors include Manuel Jesus Castro Diaz, Marc de la Asuncion (both from University of Malaga) and Yong Wei (from NOAA/University of Washington).

## Scientific background

The model  solves the PDE system for the Shallow-Water equations with a coupled transport equation for the bathymetry function. The numerical scheme is based on the path-conservative technique. You can see the full paper describing the model in the following link: [https://repository.library.noaa.gov/view/noaa/67438/noaa_67438_DS1.pdf](https://repository.library.noaa.gov/view/noaa/67438/noaa_67438_DS1.pdf)

## Contact

For questions or collaboration you can find me at ernestog@uw.edu. You can report bugs or request new features through this Github page.

## Cite

You can cite this code with the following reference:

Fernández, E. G., Díaz, M. C., Wei, Y., & Moore, C. (2024). Modeling sediment movement in the shallow-water framework: A morpho-hydrodynamic approach with numerical simulations and experimental validation. Ocean Modelling, 192, 102445.
Or visit [this web](https://www.sciencedirect.com/science/article/pii/S146350032400132X?via%3Dihub) to download the reference.

## License
This project is licensed under the **GNU Affero General Public License v3.0** - see the [LICENSE](LICENSE) file for details.
