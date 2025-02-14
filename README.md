# Rayleigh method adapted for the study of the optical response of natural photonic structures

Code to reproduce the simulations and figures in the following paper:

**Vidal, M.S.**, Dolinko, A.E. & Skigin, D.C. Rayleigh method adapted for the study of the optical response of natural photonic structures. *Eur. Phys. J. E* 44, 118 (2021).

This study explores how the [Rayleigh method](https://royalsocietypublishing.org/doi/10.1098/rspa.1907.0051), traditionally used in wave scattering problems for periodic structures, can be adapted to analyze the optical properties of natural photonic structures, which present irregularities and have complex shapes. The key aspect of the developed formalism is its capability of introducing the interface profile within the code by means of a digitalized
image of the structure, which can be either obtained from an electron microscopy image or simply by design according to the complexity of the scattering surface.

See [paper](https://link.springer.com/article/10.1140/epje/s10189-021-00124-8) for more details. 

## Project Roadmap
### 1. Extracting Profile Coordinates from an image
* The folder "tutorial_imagej" contains a step-by-step explanation of how to extract profile coordinates from an image using ImageJ.
* The tutorial is demonstrated with a TEM image of the butterfly *Dione vanillae*, which is included in the folder for testing.
### 2. Processing Coordinates with Python
* The folder "tutorial_image_python" provides an example of how to process coordinates extracted with ImageJ for an *Euglena* image. It includes:
    - A Python script (.py) for processing the coordinates.
    - A text file (.txt) with the extracted coordinates from ImageJ.
### 3. Plane Wave Illumination
The folder "plane_wave" contains all the necessary scripts for simulations with plane wave illumination.
Steps to Run the Plane Wave Simulation:
  * Run script “funciones_perfil_analitico_onda_plana.py” or “funciones_perfil_imagen_onda_plana.py”.
    - funciones_perfil_analitico_onda_plana.py → Use this if the profile is defined by an analytical function.
    - funciones_perfil_imagen_onda_plana.py → Use this if the profile is defined by coordinates extracted from an image. Import a txt file with the coordinates. A sample file, euglenido_paper_listo.txt, is included for testing.
  * Run the script to set up the linear system
    - If you want to compute reflected fields, run "R_onda_plana.py" (grid parameters can be adjusted).
    - If you want to compute transmitted fields, run "T_onda_plana.py" following the same approach.
  * Energy Balance Check
     - After computing the reflected and transmitted fields, verify energy conservation by running "chequeo_energ_onda_plana_diel_sin_perdidas.py". The "e" value in the output represents the percentage error.
### 4. Gaussian Beam Illumination
The folder "gaussian_beam" contains the same types of scripts as "plane_wave", but adapted for simulations with a Gaussian beam as the incident wave.
