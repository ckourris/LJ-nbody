# Lennard Jones N-body Simulation

Our program simulates N-body systems interacting through the Lennard-Jones pair potential. The simulation uses periodic boundary conditions and the
minimum image convention.

## Getting Started

These instructions will give you an outline of the program and allow you to run the simulation on your local machine for development and testing purposes. The program comes with an optional *accelerated* version which is recommended. More on this in the **Installing** section.

### Modules

* ```Main.py``` The main method which is called to run the simulation.
* ```Particle3D.py```, ```Box.py``` The classes that are used throughout the simulation.
* ```MDUtilities.py```,```Utilities.py``` Modules that contain the initializations, integrators and analysis functions.
* ```accelerate.cpp``` C++ Accelerated functions module.
* ```accelerate_module.pyx``` The Cython3 wrapper to translate the C++ code into Python3.
* ```accelerate.h``` C++ headers file.
* ```compilec.py``` Python module to compile the C++ library.



### Prerequisites

The code is written in  ```python3``` and is therefore required for the program to run. There is an optional accelerated version which uses some c++ modules. We have supplied the compiled acceleration library so that the code will work immediately on another linux system. If it does not it needs to be recompiled. In order to do this you need ```cython3``` to wrap the ```C++``` code with python. See next section for information.


### Compiling C++

It is recommended that you use the C++ code because it runs the simulation around 30 times faster. To ensure that you have ```cython3``` use the following command:

```
sudo apt-get install cython3
```

In case the C++ library needs to be compiled before running, this is easily done by:

```
python3 compilec.py build_ext --inplace
```

Note that the above command produces very verbose output. This completes the compilation process.

## Running the tests

To run the simulation, use e.g. the ```solid.txt``` file with the desired parameters and run:

```
 python3 Main.py parameters.txt vmdoutput.xyz
```
Note that the above Python only mode is very slow. For the optional accelerated version run (this is recommended):

```
 python3 Main.py parameters.txt vmdoutput.xyz -a
```


This will produce the following files:
* ```vmdoutput.xyz``` which contains the positions of the files for every timestep. It can be loaded directly into VMD for visualization.
* ```energyfile.txt``` which contains the timestep, kinetic energy, potential energy and total energy.
* ```MSD_output.txt``` which contains the timestep and MSD values.
* ```RDF_output.txt``` which contains the timestep and RDF values.

The program will also create plots of MSD, RDF and the energies and save them in the ```Plots``` directory.

### Observables calculations

The Main method ends by the Means Square Displacement and Radial Distribution Function and Energies calculations. Those produce a plot each depending on the range of times needed.

You can alter the parameters by modifying the following part in ```Main.py```:

```
  msd_start = 7000
  msd_end = 9999
  rdf_bins = np.arange(0,5,0.1)
  rdf_start = 9900
  rdf_end = 9999
```
The ```rdf_bins``` creates an array for the radii to be plotted on the RDF diagram.

It is possible to do further tests without running the simulation again. To load the data (positions) of any data file in the simulation use the ```get_output``` method in ```Utilities.py``` and comment out the simulation run of the box. Replace the outfile string with the desired file name and uncomment the following lines in the code:

```
outfile = ""
position_list = np.array(get_output(outfile, parameters[0]))
```


## Authors

* **Christos Kourris** - [ckourris](https://github.com/ckourris)

* **Ethan van Woerkom** - [Edekje](https://github.com/Edekje)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details. The project is hosted on the following github repository:
[LJ-nbody](https://github.com/Edekje/LJ-nbody)
