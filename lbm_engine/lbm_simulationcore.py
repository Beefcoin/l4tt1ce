import numpy as np

class Lattice2D:
    def __init__(self, nx, ny, descriptor, collisionOperator):
        #define the lattice dimensions in lattice units
        self.nx, self.ny = nx, ny 
        self.descriptor = descriptor
        #set the operator used for all colission on this lattice
        self.collisionOperator = collisionOperator
        #get the velocity directions from our lattice descriptor
        self.e = descriptor.e
        #get the direction opposites from descriptor definition
        self.opp = descriptor.opp
        #get the unique velocity count from descriptor definition
        self.Q = descriptor.Q

        ##Setup of all the fields needed
        self.X, self.Y = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')
        self.geometry = {}
        #setup density at 1
        self.rho = np.ones((ny, nx))
        #setup velocity as 0 in all directions (2 for xy)
        self.u = np.zeros((ny, nx, 2))
        #setup distribution function as 0 everywhere for all Q directions
        self.f = np.zeros((ny, nx, self.Q))
        #setup of timeStep variable
        self.t = 0

        ##setup the distribution function as equilibrium using the collision operator
        print("Type of collisionOperator:", type(self.collisionOperator))
        self.f = self.collisionOperator.compute_feq(self.descriptor, self.rho, self.u)


    def addOperator(self, name, operator):
        self.geometry[name] = operator

    

    def step(self):
        ##Streaming Step
        #for all Q directions we stream the distribution function across the lattice using our discrete velocity set ei
        #looks complicated af but is actually not that hard
        #this is periodic -> if you dont want this you have to continously overwrite the boundary conditions
        for i in range(self.Q):
            self.f[:, :, i] = np.roll(
                np.roll(self.f[:, :, i], self.e[i, 0], axis=1),
                self.e[i, 1], axis=0
            )
        
        ##Here we iterate over the operators in the operator List and call the apply method for each one to apply boundary conditions
        for operator in self.geometry.values():
            #print(operator)
            operator.apply(self.f, self.u, self.rho)
        
        #calculate density from distribution function 
        self.rho = np.sum(self.f, axis=2)
        #calculate x component of velocity set from the distribution function (use direction vectors as weights (TM1 flashback))
        self.u[:, :, 0] = np.sum(self.f * self.e[:, 0], axis=2) / self.rho
        #same for y component
        self.u[:, :, 1] = np.sum(self.f * self.e[:, 1], axis=2) / self.rho


        ##Collision Step
        #Calculate delta f from collision Operator
        self.f += self.collisionOperator.compute_delta_f(self.descriptor, self.f, self.rho, self.u)

        #We are done with this Stream&Collide so increment timeStep
        self.t += 1




class ScalarLattice2D:
    #this is pretty much the same as the regular Lattice but it contains a scalar field used for diffusive transport etc
    def __init__(self, nx, ny, descriptor, collisionOperator):
        #define the lattice dimensions in lattice units
        self.nx, self.ny = nx, ny 
        self.descriptor = descriptor
        #set the operator used for all colission on this lattice
        self.collisionOperator = collisionOperator
        #get the velocity directions from our lattice descriptor
        self.e = descriptor.e
        #get the direction opposites from descriptor definition
        self.opp = descriptor.opp
        #get the unique velocity count from descriptor definition
        self.Q = descriptor.Q

        ##Setup of all the fields needed
        self.X, self.Y = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')
        self.geometry = {}
        #setup velocity as 0 in all directions (2 for xy)
        self.u = np.zeros((ny, nx, 2))
        #setup scalar field as 1 everywhere
        self.phi = np.ones((ny, nx))
        #setup distribution function as 0 everywhere for all Q directions
        self.g= np.zeros((ny, nx, self.Q))
        #setup of timeStep variable
        self.t = 0

        ##setup the distribution function as equilibrium using the collision operator
        print("Type of collisionOperator:", type(self.collisionOperator))
        self.g = self.collisionOperator.compute_feq(self.descriptor, self.phi, self.u)

    def addOperator(self, name, operator):
        self.geometry[name] = operator

    def step(self):
        for i in range(self.Q):
            self.g[:, :, i] = np.roll(np.roll(self.g[:, :, i], self.e[i, 0], axis=1), self.e[i, 1], axis=0)

        self.phi = np.sum(self.g, axis=2)

        for operator in self.geometry.values():
            operator.apply(self.g, self.u, self.phi)

        self.g += self.collisionOperator.compute_delta_f(self.descriptor, self.g, self.phi, self.u)
        #We are done with this Stream&Collide so increment timeStep
        self.t +=1

def PrintLatticeInformation(nslattice, adlattice=None):
    print("\n" + "=" * 60)
    print("LATTICE INFORMATION".center(60))
    print(f"Time step: t = {nslattice.t}".center(60))
    print("-" * 60)

    # --- Navier-Stokes Lattice ---
    print("Navier-Stokes Lattice".center(60))
    velocity_magnitude = np.sqrt(nslattice.u[..., 0]**2 + nslattice.u[..., 1]**2)
    avg_velocity = np.mean(velocity_magnitude)
    avg_rho = np.mean(nslattice.rho)

    print(f"{'Average Velocity Magnitude (u)':<40}: {avg_velocity:10.3f}")
    print(f"{'Average Density (ρ)':<40}: {avg_rho:10.3f}")
    
    # --- Advection-Diffusion Lattice ---
    if adlattice is not None:
        print("-" * 60)
        print("Advection-Diffusion Lattice".center(60))
        avg_phi = np.mean(adlattice.phi)
        min_phi = np.min(adlattice.phi)
        max_phi = np.max(adlattice.phi)

        print(f"{'Average Scalar Field (φ)':<40}: {avg_phi:10.3f}")
        print(f"{'Minimum Scalar Field (φₘᵢₙ)':<40}: {min_phi:10.3f}")
        print(f"{'Maximum Scalar Field (φₘₐₓ)':<40}: {max_phi:10.3f}")

    print("=" * 60 + "\n")
