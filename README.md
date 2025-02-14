# Rayleigh method adapted for the study of the optical response of natural photonic structures

Code to reproduce the simulations and figures in the following paper:

**Vidal, M.S.**, Dolinko, A.E. & Skigin, D.C. Rayleigh method adapted for the study of the optical response of natural photonic structures. *Eur. Phys. J. E* 44, 118 (2021).

This study explores how the [Rayleigh method](https://royalsocietypublishing.org/doi/10.1098/rspa.1907.0051), traditionally used in wave scattering problems for periodic structures, can be adapted to analyze the optical properties of natural photonic structures, which present irregularities and have complex shapes. The key aspect of the developed formalism is its capability of introducing the interface profile within the code by means of a digitalized
image of the structure, which can be either obtained from an electron microscopy image or simply by design according to the complexity of the scattering surface.

See [paper](https://link.springer.com/article/10.1140/epje/s10189-021-00124-8) for more details. 

## Organization

 - The simulations of Figure 5, 6 and 7 of the paper are generated from figures_paper_step.ipynb.
 - The simulations of Figure 9, 10 and 11 of the paper are generated from graphs_paper.ipynb.

## Simulations in C++
The chemotaxis assay simulations run in C++. To run them:

1. Compile the files (with make)

```
mingw32-make
```

2. Run the executable that was generated in the previous step

```
./main
```

3. If you want to delete the .o and .exe files and re-compile, run:

```
mingw32-make clean
```

## Other considerations

### main.cpp
In the main.cpp file you can choose running example trajectories (#define trajectories). In this case, you will save a txt file with the column values detailed in the function PrintDetail in WormAgent.cpp for every time step of the simulation. The second option is to run several simulations (#define chem_index) and save only the last time point of the simulation to compute the chemotaxis index afterwards (with the chemotaxis_index.ipynb notebook). In this case, also a txt file is generated from the PrintDetail function.
 
The NaCl concentration during cultivation can be changed with the parameter "phenotype(25)" (C_breed). We explored the values C_breed=25,50,100.

### WormAgent.h
You can choose between a step of NaCl (#define GRAD_STEP), or the gradient use in the paper (#define GRAD_GAUS). Beware that if you use the step concentration, you need to change the constant Preexposure to a value lower than 5000 for memory reasons.
