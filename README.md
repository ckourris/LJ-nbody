# Lennard Jones N-body Simulation

Our program simulates N-body systems interacting through the Lennard-Jones pair potential. The simulation uses periodic boundary conditions and the
minimum image convention.

## Getting Started

These instructions will give you an outline of the program and allow you to run the simulationon your local machine for development and testing purposes. The program comes with an optional *accelerated* version which is recommended. More on this in the **Installing** section.

### Modules

* ```Main.py``` The main method which is called to run the simulation.
* ```Particle3D.py```, ```Box.py``` The classes that are used throughout the simulation.
* ```MDUtilities.py```,```Utilities.py``` Modules that contain the initializations, integrators and analysis functions.
* ```Cforces.cpp``` Elementary calculations module in C++
* ```Cforces.pyx``` The Cython3 wrapper to translate the C++ code into Python3.
* ```getforces.h``` The newly created library with the C++ functions.
* ```compilec.py``` Python module to compile the C++ library.



### Prerequisites

The code is written in  ```python3``` and is therefore required for the program to run. There is an optional accelerated version which uses some c++ modules. In order to use it you need ```cython3``` to wrap the c++ code with python. See next section for information.


### Installing

It is recommended that you use the C++ code because it runs the simulation around 30 times faster. Firstly download ```cython3``` using the command:

```
sudo apt-get install cython3
```

The C++ library needs to be compiled before run. This is very easily done by:

```
python3 compilec.py build_ext --inplace
```


## Running the tests

To run the simulation, adjust the ```parameters.txt``` file to the desired parameters and run:

```
 python3 Main.py parameters.txt vmdoutput.xyz
```
For the optional accelerated version run:

```
 python3 Main.py parameters.txt vmdoutput.xyz -a
```


This will produce the following files:
* ```vmdoutput.xyz``` which contains the positions of the files for every timestep. It can be loaded directly into VMD for visualization.
* ```energyfile.txt``` which contains the timestep, kinetic energy, potential energy and total energy.

### Observables calculations

The Main method ends by calling the Means Square Displacement and Radial Distribution Function calculation. Those produce a plot each depending on the range of times needed.

You can alter the parameters by modifying the following part in ```Main.py```:

```
  msd_start = 1    
  msd_end = 200
  rdf_bins = np.arange(0,5,0.1)
  rdf_start = 100
  rdf_end = 200
```
The ```rdf_bins``` creates an array for the radii to be plotted on the RDF diagram.

You can proceed with further tests without running the simulation again. To load the data (positions) of any data file in the simulation use the ```get_output``` method in ```Utilities.py``` by commenting the creation of the creation of the box and its subsequent run. Then uncomment the following and fill in the outfile string with the desired one.

```
outfile = ""
position_list = np.array(get_output(outfile, parameters[0]))
```

## Deployment

Add additional notes about how to deploy this on a live system

## Authors

* **Christos Kourris** - [ckourris](https://github.com/ckourris)

* **Ethan van Woerkom** - [Edekje](https://github.com/Edekje)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details. The project is hosted on the following github repository:
[LJ-nbody](https://github.com/Edekje/LJ-nbody)
