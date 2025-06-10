# MRI Metal Artifact Simulation and Mitigation

This repository contains scripts and notebooks used to simulate metal artifacts in MRI images and train neural networks to remove them. The workflow is:

1. **Simulation** – The MATLAB script `try_numpy_mac.m` converts real MRI slices stored as `.npy` files into synthetic images with simulated metal artifacts using a Fourier‑based approach. It is adapted from [this MathWorks file](https://www.mathworks.com/matlabcentral/fileexchange/56680-mri-simulation-using-forecast-fourier-based-off-resonance-artifact-simulation-in-the-steady-state). Each processed slice produces a pair of arrays: the original proton density slice and the corresponding artifact‑corrupted magnitude image. Output arrays are saved under `outputs/<ID>/`.
2. **Autoencoder** – `Final_project_Autoencoder.ipynb` loads the paired arrays, visualizes sample data, and trains a fully connected autoencoder in PyTorch. The model learns to map artifact‑affected inputs back to their artifact‑free counterparts.
3. **U‑Net** – `Final_Project_Unet.ipynb` implements a custom U‑Net architecture for the same restoration task. It defines the network, prepares data loaders, and demonstrates training on the simulated dataset.

The resulting models aim to reduce metal artifacts in MRI scans across various scenarios such as head & neck or orthopedic imaging.
