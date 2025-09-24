## Project Status
This repository is archived and not under active development.

# L4tt1ce

A lightweight **2D**/**3D** Python implementation of the Lattice Boltzmann Method (LBM).

This project uses the **BGK** (Bhatnagar-Gross-Krook) operator to simulate:

- Navier-Stokes fluid dynamics
- Scalar transport (advection-diffusion)

## Requirements
Install required libraries via:

```bash
pip install numpy matplotlib pyevtk
```

## Feature Overview

| Feature                          | Status   | Description                                        |
|----------------------------------|----------|----------------------------------------------------|
| 2D Simulation                    | ✅ Done  | Current implementation supports 2D flows           |
| Dirichlet & Neumann BCs         | ✅ Done  | Velocity, pressure, concentration, flux boundaries |
| Coupled Scalar Transport         | ✅ Done  | Advection-diffusion modeling                       |
| 3D Simulation                    | ✅ Done   | Load geometry from VTK files                        |
| 3D Scalar Transport                   | 🔜 Todo    | Scalar transport in 3D is currently not implemented                        |
| MRT/TRT Models                   | 🔜 Todo  | Advanced collision operators coming soon                       |

## Showcases
### 3D Porous Media Simulation
![me](https://media4.giphy.com/media/v1.Y2lkPTc5MGI3NjExb204bXdmbnFhaWd3OW1sNmE4d3A4MXdrYnZ5ZzU3Z3FhejE5MmsxZiZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9Zw/xuj1xJyPcAeX6Q7bUE/giphy.gif)
### Flow around obstacles
![me](https://media1.giphy.com/media/v1.Y2lkPTc5MGI3NjExdmE5N2tnaDZjN2Jkd3JvaXBxZXMwbTN0YXhtZXB6dTM4bmQwb2UzZiZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9Zw/gv1LI63WFKrTyV4Rcf/giphy.gif)
### Heated obstacle in channel
![me](https://media1.giphy.com/media/v1.Y2lkPTc5MGI3NjExancxc3d4Mmc0YWJxaG5hMGs4eHY3cWZ3ZmxkOTF3dnQ4Z2Z5ZnprYyZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9Zw/8oFlu2HEWHScOQqiny/giphy.gif)




