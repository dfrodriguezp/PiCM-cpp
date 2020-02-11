# PiCM (C++ version)

Software for simulation of the dynamics of charged particles with periodic boundary conditions using the particle-in-cell method.

## Requirements

You need the `g++` compiler, along with the `jsoncpp` library. In Ubuntu, these can be installed by typing the following command in the console

```bash
$ sudo apt install g++ libjsoncpp-dev
```

Additionally, you need **Python3** with the following modules:

- NumPy
- Click
- Matplotlib

These can be installed via

```bash
$ pip install numpy click matplotlib
```

or

```bash
$ conda install numpy click matplotlib
```

<u>**Note: unfortunately, the current version of the program only works on Ubuntu.**</u>

## Running the program

1. Clone or download the repository.
2. Compile the program with

```bash
$ make
```

3. Run a simulation with

```bash
$ ./main path/to/jsonfile.json
```

Of course, you must know how to build that JSON file and what it contains. This is further explained in the [Wiki](https://github.com/dfrodriguezp/PiCM_cpp/wiki) of this repository.

### Example

Let's perform a simulation of a two-stream instability. For this, there is a script called `build_two_stream.py` which creates the sample file as well as the `.json` file. You can set the parameters you prefer. **Note:** for this example, let's pretend you didn't modified it.

1. Open a terminal inside the repository folder.
2. Run the sample script

```bash
$ python3 build_two_stream.py
```

This will create a folder with the name specified in the `output` variable. Inside that folder you can find both the sample and the `.json` file. The sample file is called `two_stream.dat` and the JSON file `sim_two_stream.json`.

3. Run the simulation

```bash
$ ./main electrostatic/sim_two_stream.json
```

4. Plot the energy

```bash
$ python3 plotters/plot_energy.py electrostatic/sim_two_stream.json
```

![Energy](example_imgs/energy.png)

5. Plot the phase space in the step 150

```bash
$ python3 plotters/plot_phase_space.py electrostatic/sim_two_stream.json 150
```

![Phase_space](example_imgs/step_150_x_.png)
