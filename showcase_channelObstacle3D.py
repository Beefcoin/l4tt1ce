import numpy as np
from lbm_engine.lbm_descriptors import D3Q19
from lbm_engine.lbm_geometry import load_mask_from_vti
from lbm_engine.lbm_collisionOperators import BGK_collisionOperator3D
from lbm_engine.lbm_simulationcore import Lattice3D
from lbm_engine.lbm_operators import BounceBack3D, PressureDirichlet3D, VelocityDirichlet3D
from lbm_engine.lbm_visualize import export_fields_vti3D

#Basic Simulation Parameters
nx, ny, nz = 100, 100, 100   #sim dimensions
timesteps = 600
u_max = 0.4
plot_interval = 1

#basic simulöation setup
descriptor = D3Q19()
collision = BGK_collisionOperator3D(tau=0.75)
sim = Lattice3D(nx=nx, ny=ny, nz=nz, descriptor=descriptor, collisionOperator=collision)

#easy, get the walls by geomtry definitons from our dimensions, kinda chill, shoutout OpenLB
X, Y, Z = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
wall_mask = (Y == 0) | (Y == ny - 1) | (Z == 0) | (Z == nz - 1)

#load geomtry from file using the load_mask_from_vti method
obstacle_mask = load_mask_from_vti("geometryTest_coherent_porous.vti", field_name="Tiff Scalars").astype(bool)
sim.addOperator("obstacles", BounceBack3D(descriptor, obstacle_mask))

# set constant inlet velocity (just for velocity dirichlet)
inlet_mask = (X == 0)
def inlet_velocity(shape):
    vel = np.zeros(shape)
    vel[..., 0] = u_max  # set flow in x-direction
    return vel

#add velocity dirichlet or pressure dirichlet to inlet
#sim.addOperator("inlet", VelocityDirichlet3D(descriptor, collision, inlet_mask, inlet_velocity))
sim.addOperator("inlet", PressureDirichlet3D(descriptor, collision, inlet_mask, rho_value=2.0))
#add pressure dirichlet to outlet
outlet_mask = (X == nx - 1)
sim.addOperator("outlet", PressureDirichlet3D(descriptor, collision, outlet_mask, rho_value=1.0))

# --- Main Simulation Loop ---
for t in range(timesteps):
    sim.step()
    if t % plot_interval == 0:
        print(f"Timestep {t} / {timesteps}")
        export_fields_vti3D(t, sim.u, phi=sim.rho, output_dir="ShowCase3D_FlowAroundObject")
