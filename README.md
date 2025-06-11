# L4tt1ce

A lightweight **2D** Python implementation of the Lattice Boltzmann Method (LBM).

This project uses the **BGK** (Bhatnagar-Gross-Krook) operator to simulate:

- Navier-Stokes fluid dynamics
- Scalar transport (advection-diffusion)

---

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
| 3D Simulation                    | 🚧 WIP   | Planned for future release                         |
| Geometry from Images             | 🔜 Todo  | Planned feature                                    |
| MRT/TRT Models                   | 🔜 Todo  | Advanced collision operators                       |

## Showcases
### Flow around obstacles
![me](https://media1.giphy.com/media/v1.Y2lkPTc5MGI3NjExdmE5N2tnaDZjN2Jkd3JvaXBxZXMwbTN0YXhtZXB6dTM4bmQwb2UzZiZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9Zw/gv1LI63WFKrTyV4Rcf/giphy.gif)
### Heated obstacle in channel
![me](https://media1.giphy.com/media/v1.Y2lkPTc5MGI3NjExancxc3d4Mmc0YWJxaG5hMGs4eHY3cWZ3ZmxkOTF3dnQ4Z2Z5ZnprYyZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9Zw/8oFlu2HEWHScOQqiny/giphy.gif)




