# PiCM
## 2D3V Electrostatic Particle-in-Cell Code

This code was developed in the PCM Computational Applications research group, at the Universidad Nacional de Colombia. The aim of the code is to be a tool for the simulation of different plasma phenomena. I write the code for it to be easy to understand, and easy to use to begginers in the computational physics area.

### Pre-requisites

The user shall provide a `JSON` file, like the one shown below, wich contains all of the information needed to start the simulation.

```
{
  "N": 10000,
  "steps": 1000,
  "grid_size": [
    64,
    64
  ],
  "ss_frequency": 10,
  "dt": 0.05,
  "dx": 1.0,
  "dy": 1.0,
  "sample": "simulation_name/sample.dat",
  "seed": 53668301,
  "output": "simulation_name",
  "results": [
    "phase_space",
    "electric_potential"
  ],
  "Bfield": [
    0.0,
    0.0,
    0.0
  ]
}
```
It should be noted that both the `JSON` file and the sample file, must be inside the same folder **simulation_name**.

Furthermore, in the `"results"` option, the user can specify the data output that they want to work with.

### Instructions
1. Clone the repository.
2. Compile it with `$ make`.
3. Run the simulation with `$ ./main simulation_name/example.json`.

When the simulation is completed, the output data will be available in the `"output"` directory specified in the `JSON` file --in this case, *simulation_name*--. The user can manipulate it as their convenience.

Feel free to send your suggestions to improve the code. Remember that this is an educational project, aim to new computational physics students.
