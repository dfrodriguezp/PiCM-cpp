# PiCM
## 2D3V Electrostatic Particle-in-Cell Code

This code was developed in the PCM Computational Applications research group, at the Universidad Nacional de Colombia. The aim of the code is to be a tool for the simulation of different plasma phenomena. I write the code for it to be easy to understand, and easy to use to begginers in the computational physics area.

### Usage

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
  "dr": 1.0,
  "sample": "example/sample.dat",
  "seed": 53668301,
  "output": "example",
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

### Instructions
1. Clone the repository.
2. Compile it with `$ make`.
3. Run the simulation with `$ ./main example/example.json`.

When the simulation is completed, the output data will be available in the `"output"` directory specified in the `JSON` file --in this case, *example*--. The user can manipulate it as their convenience.

Feel free to send your suggestions to improve the code. Remember that this is an educational project, aim to new computational physics students.
