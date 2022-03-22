# SIEVE_simulator
The simulator for SIEVE is modified based on [CellCoal](https://github.com/dapogon/cellcoal) v1.0.0. It is used in the SIEVE paper to generate simulated datasets.

## Compilation

SIEVE_simulator uses cmake to build the executable.

```bash
$ mkdir build
$ cd build
$ cmake ..
$ cmake --build .
```

The built executable is located: `build/bin/SIEVE_simulator`.

## Run

To run SIEVE_simulator, use the example configuration file `parameters`:

```bash
$ ./build/bin/SIEVE_simulator -Fparameters
```

## Configuration files

SIEVE_simulator uses a similar input as CellCoal does. Please refer to the manuals of CellCoal. Those configuration files used in the SIEVE paper can be found in the [SIEVE_benchmark_pipeline](https://github.com/szczurek-lab/SIEVE_benchmark_pipeline) repository, under directory `simulation_configs`.
