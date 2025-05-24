# Particle simulation

Compilation with `make` (uses `zsh`). Needs `gfortran`(i used the one from `gcc 14.2` complier collection).

Run with:

```
./simu environment_file.dat particle_files/big_body_test.dat test.output
```

Play results with:

```
python python_src/play_simulation.py -i test.output
```
