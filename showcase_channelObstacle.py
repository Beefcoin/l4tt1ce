import numpy as np
from lbm_engine.lbm_descriptors import D2Q9
from lbm_engine.lbm_collisionOperators import BGK_collisionOperator
from lbm_engine.lbm_simulationcore import Lattice2D
from lbm_engine.lbm_operators import BounceBack, PressureDirichlet, VelocityDirichlet
from lbm_engine.lbm_geometry import create_triangle_mask, create_rectangle_mask
from lbm_engine.lbm_visualize import visualize, export_fields_vti


""" This is a pretty basic showcase: It creates a Navier Stokes flow through a Channel with two rectangular obstacles.
It creates matplotlib plots as well as VTK Data. 
For the inlet a constant velocity dirichlet is applied, the outlet is modeled using a constant pressure Dirichlet.
 """

#Geometry Setup
nx, ny = 400, 150
#Simulation Setup
timesteps = 20000
u_max = 0.1
plot_interval = 100

#Set up Basic Navier Stokes Simulation
descriptorNS = D2Q9()
collisionOperatorBGK = BGK_collisionOperator(tau=0.56)
sim = Lattice2D(
    nx=nx,
    ny=ny,
    descriptor=D2Q9(),
    collisionOperator=collisionOperatorBGK,
)

X, Y = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
#Here we create a wall mask for applying bounceback operator the the simulation borders (channel modeling)
top_wall_mask = (Y == ny - 1)
bottom_wall_mask = (Y == 0)
wall_mask = top_wall_mask | bottom_wall_mask

#Here we create some obstacles and combine their mask with the wallmask
mask1 = create_rectangle_mask(X, Y, center_x=60, center_y=25, width=50, height=50)
mask2 = create_rectangle_mask(X, Y, center_x=160, center_y=ny-40, width=50, height=80)
combined_mask = mask1 | mask2 

#apply bounceback to the generalized total mask
combined_mask_total = combined_mask | wall_mask
sim.addOperator("obstacle", BounceBack(sim.descriptor, combined_mask_total))


# Constant Velocity inlet on the left (Dirichlet)
inlet_mask = X == 0
def inlet_velocity(shape):
    return np.zeros(shape) + [u_max, 0]
sim.addOperator("inlet", VelocityDirichlet(sim.descriptor, collisionOperatorBGK, inlet_mask, inlet_velocity))

# Consant Pressure Dirichlet outlet on the right
outlet_mask = X == nx - 1
sim.addOperator("outlet", PressureDirichlet(sim.descriptor, collisionOperatorBGK, outlet_mask, rho_value=1.0))

for t in range(timesteps):
    sim.step()
    if t % plot_interval == 0:
            print("Timestep " + str(t) + " / " + str(timesteps))
            export_fields_vti(t, u=sim.u, output_dir="ShowCase_FlowAroundObjects")
            visualize(sim, t, 0.23)
