# LJ-nbody
Lennard Jones N-body simulation with periodic boundary conditions to simulate properties of molecular argon.

# Lennard Jones N-body Simulation

Our program simulates N-body systems interacting through the Lennard-Jones pair potential. The simulation uses periodic boundary conditions and the
minimum image convention.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Modules

* ```Main.py``` The main method which is called to run the simulation.
* ```Particle3D.py```, ```Box.py``` The classes that are used throughout the simulation.
* ```MDUtilities.py```,```Utilities.py``` Modules that contain the initializations, integrators and analysis functions.
* ```Cforces.cpp``` Elementary calculations module in C++
* ```Cforces.pyx``` The Cython3 wrapper to translate the C++ code into Python3.
* ```getforces.h``` The newly created library with the C++ functions.
* ```compilec.py``` Python module to compile the C++ library.



### Prerequisites

The code is written in  ```python3``` and is therefore required for the program to run. Part of the modules are written in ```C++``` which has been wrapped to Python using ```cython3```. To install ```cython3``` do:

```
sudo apt-get install cython3
```

### Installing

In order to use the C++ code you need to compile the library. This is very easily done by:

```
python3 compilec.py build_ext --inplace
```

The C++ code can be avoided by changing a command in the ```Box.py``` class. Simply change Line 13 to:

```
CPP_ENABLED=False
```


## Running the tests

To run the simulation, adjust the ```parameters.txt``` file to the desired parameters and run:

```
 python3 Main.py parameters.txt vmdoutput.xyz
```

This will produce the following files:
* ```vmdoutput.xyz``` which contains the positions of the files for every timestep. It can be loaded directly into VMD for visualization.
* ```energyfile.txt``` which contains the timestep, kinetic energy, potential energy and total energy.

### Observables calculations

The Main method ends by calling the Means Square Displacement and Radial Distribution Function calculation. Those produce a plot each depending on the range of times needed.

You can alter the parameters by modifying the following commands:

```
MSD(position_list, 1,100, Simba.boxdim)
```

```
rdf_arr = RDF(position_list, 200, 300, np.arange(0,5,0.1), Simba.boxdim)
```

You can proceed with further tests without running the simulation again. To load the data (positions) of any data file in the simulation use the ```get_output``` method in ```Utilities.py``` as in:

```
position_list = np.array(get_output(outfile, parameters[0]))
```

## Deployment

Add additional notes about how to deploy this on a live system

## Authors

* **Christos Kourris** - [ckourris](https://github.com/ckourris)

* **Ethan van Woerkom** - [Edekje](https://github.com/Edekje)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
