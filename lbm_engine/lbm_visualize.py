import matplotlib.pyplot as plt
import os 
from pyevtk.hl import imageToVTK
import numpy as np

""" This file contains some basic visualization files that can create plots using matplotlib as well as VTK Files for Paraview
The visualization methods are kinda ugly and have to be rewritten because of scaling and layout but it works for now so..... """

def visualize(sim, t, vmax_set):
    vmag = np.sqrt(sim.u[:, :, 0]**2 + sim.u[:, :, 1]**2)
    if "obstacle" in sim.geometry:
        mask = sim.geometry["obstacle"].mask
        mask = np.transpose(mask) if mask.shape != vmag.shape else mask
        vmag[mask] = np.nan
    #plt.figure(figsize=(8, 5))
    u_mag = np.sqrt(sim.u[:, :, 0]**2 + sim.u[:, :, 1]**2)
    max_velocity = np.max(u_mag)
    #print(max_velocity)
    plt.imshow(vmag, cmap="coolwarm", vmin=0, vmax=vmax_set)
    plt.title(f"Velocity magnitude – t = {t}")
    plt.axis("off")
    plt.colorbar(label="|u|")
    plt.tight_layout()
    plt.savefig(f"frames/frame_{t:05d}.png", dpi=150)
    plt.close()

def visualizeScalar(sim, t):
        u_mag = np.sqrt(sim.phi[:, :]**2 + sim.phi[:, :]**2)
        max_velocity = np.max(u_mag)
        print(max_velocity)
        plt.figure(figsize=(8, 4))
        plt.imshow(sim.phi, cmap='inferno')
        plt.title(f"Scalar field ϕ at t={t}")
        plt.colorbar(label="ϕ")
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(f"frames/phi_{t:05d}.png", dpi=150)
        plt.close()

import numpy as np
import matplotlib.pyplot as plt
import os

def visualize_combined(sim, t, overall_title=None, scale = 10, dpi =150):
    # Compute velocity magnitude
    vmag = np.sqrt(sim.u[:, :, 0]**2 + sim.u[:, :, 1]**2)
    phi = sim.phi

    # Mask obstacles with NaN
    if "obstacle" in sim.geometry:
        mask = sim.geometry["obstacle"].mask
        if mask.shape != vmag.shape:
            mask = mask.T
        vmag = vmag.copy()
        phi = phi.copy()
        vmag[mask] = np.nan
        phi[mask] = np.nan

    # Get data shape and compute figure size to match aspect ratio
    ny, nx = vmag.shape
    fig_width = scale * (nx / max(nx, ny))  # normalize by max dimension
    fig_height = scale * (ny / max(nx, ny)) * 2  # times 2 because of two subplots

    fig, axs = plt.subplots(2, 1, figsize=(fig_width, fig_height), dpi=dpi)

    if overall_title:
        fig.suptitle(overall_title, fontsize=14, fontweight='bold', y=0.98)

    # Velocity magnitude plot (fixed color scale)
    im1 = axs[0].imshow(vmag, cmap="coolwarm", vmin=0, vmax=0.19)
    axs[0].set_title("Velocity Magnitude")
    axs[0].axis("off")
    fig.colorbar(im1, ax=axs[0], fraction=0.046, pad=0.04, label="|u|")

    # Scalar field plot (autoscale color)
    im2 = axs[1].imshow(phi, cmap="inferno")
    axs[1].set_title("Scalar Field ϕ")
    axs[1].axis("off")
    fig.colorbar(im2, ax=axs[1], fraction=0.046, pad=0.04, label="ϕ")

    os.makedirs("frames", exist_ok=True)
    plt.savefig(f"frames/combined_{t:05d}.png", dpi=dpi)
    plt.close()




""" def export_fields_vti(timestep, u, phi, output_dir="output_vti"):
    
    ny, nx = phi.shape
    spacing = (1.0, 1.0, 1.0) 
    origin = (0.0, 0.0, 0.0)

    # Prepare fields as (nx, ny, 1)
    # some disgusting shaping going on here, dont ask me how this works or why (it was painful)
    phi_3d = np.ascontiguousarray(phi.T[:, :, np.newaxis]) 
    u_x = np.ascontiguousarray(u[:, :, 0].T[:, :, np.newaxis])
    u_y = np.ascontiguousarray(u[:, :, 1].T[:, :, np.newaxis])
    u_z = np.zeros_like(u_x)

    os.makedirs(output_dir, exist_ok=True)

    imageToVTK(
        os.path.join(output_dir, f"fields_{timestep:05d}"),
        origin=origin,
        spacing=spacing,
        pointData={
            "phi": phi_3d,
            "velocity": (u_x, u_y, u_z),
        }
    ) """

def export_fields_vti2D(timestep, u, phi=None, output_dir="output_vti"):
    #this creates vti files for paraview for velocity and (scalar field)
    ny, nx = u.shape[:2]
    spacing = (1.0, 1.0, 1.0) 
    origin = (0.0, 0.0, 0.0)

    # Prepare fields as (nx, ny, 1)
    # some disgusting shaping going on here, dont ask me how this works or why (it was painful)
    u_x = np.ascontiguousarray(u[:, :, 0].T[:, :, np.newaxis])
    u_y = np.ascontiguousarray(u[:, :, 1].T[:, :, np.newaxis])
    u_z = np.zeros_like(u_x)

    point_data = {
        "velocity": (u_x, u_y, u_z)
    }
    if phi is not None:
        phi_3d = np.ascontiguousarray(phi.T[:, :, np.newaxis])
        point_data["phi"] = phi_3d

    os.makedirs(output_dir, exist_ok=True)

    imageToVTK(
        os.path.join(output_dir, f"fields_{timestep:05d}"),
        origin=origin,
        spacing=spacing,
        pointData=point_data
    )

def export_fields_vti3D(timestep, u, phi=None, output_dir="output_vti"):
    # u: (nz, ny, nx, 3) – velocity field
    # phi: (nz, ny, nx) – optional scalar field (e.g., density, temperature)
    nz, ny, nx = u.shape[:3]
    spacing = (1.0, 1.0, 1.0)
    origin = (0.0, 0.0, 0.0)

    # Transpose from (z, y, x) to (x, y, z) fdor VTK
    u_x = np.ascontiguousarray(np.transpose(u[:, :, :, 0], (2, 1, 0)))
    u_y = np.ascontiguousarray(np.transpose(u[:, :, :, 1], (2, 1, 0)))
    u_z = np.ascontiguousarray(np.transpose(u[:, :, :, 2], (2, 1, 0)))

    point_data = {
        "velocity": (u_x, u_y, u_z)
    }

    if phi is not None:
        phi_3d = np.ascontiguousarray(np.transpose(phi, (2, 1, 0)))
        point_data["phi"] = phi_3d

    os.makedirs(output_dir, exist_ok=True)

    imageToVTK(
        os.path.join(output_dir, f"fields_{timestep:05d}"),
        origin=origin,
        spacing=spacing,
        pointData=point_data
    )


