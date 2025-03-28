import numpy as np

""" This file contains Collison Operators like BGK and Advection Diffusion BGK 
The implementation is pretty simple and you should be able to add your own 
Operators based on the main CollisionOperator """

class CollisionOperator:
    def __init__(self, tau):
        pass
    def compute_feq(self, rho=None, u=None,mask = None):
        pass

class BGK_collisionOperator(CollisionOperator):
    def __init__(self, tau):
        self.tau = tau
    def compute_feq(self, descriptor, rho, u, mask=None):
        ny, nx = rho.shape
        Q = descriptor.Q
        e = descriptor.e
        w = descriptor.w

        u2 = u[:,:,0]**2 + u[:,:,1]**2
        feq = np.zeros((ny, nx, Q))

        for i in range(Q):
            #scalar product for convective term
            cu = u[:,:,0]*e[i,0] + u[:,:,1]*e[i,1]
            #lbm equilibrium for equilibrium distribution function
            feq[:,:,i] = w[i] * rho * (1 + 3*cu + 4.5*cu**2 - 1.5*u2)

        #if mask is specified -> just apply operator to masked area
        if mask is not None:
            feq_masked = np.zeros_like(feq)
            for i in range(Q):
                feq_masked[:,:,i][mask] = feq[:,:,i][mask]
            return feq_masked
        return feq
    
    def compute_delta_f(self, descriptor, f_lattice, rho, u, mask=None):
        """ Relaxation Approach, this implementation kinda sucks atm because for 
        other collision operators (TRT, MRT) the Lattice2D Structure has to be changed
        -- but it works atm so idgaf """
        feq = self.compute_feq(descriptor, rho, u, mask)
        delta_f = -(1.0 / self.tau) * (f_lattice - feq)
        return delta_f
    

class BGK_AdvectionDiffusion_collisionOperator(CollisionOperator):
    def __init__(self, tau):
        self.tau = tau
    def compute_feq(self, descriptor, phi, u, mask = None):
        Q = descriptor.Q
        e = descriptor.e
        w = descriptor.w

        feq = np.zeros((*phi.shape, Q))
        for i in range(Q):
            #scalar product for convective term
            cu = u[:, :, 0]*e[i, 0] + u[:, :, 1]*e[i, 1]
            #equilibrium distribution function for diffusion problem
            feq[:, :, i] = w[i] * phi * (1 + 3*cu)
        return feq
    
    def compute_delta_f(self, descriptor, f_lattice, phi, u, mask = None):
        """ Relaxation Approach, this implementation kinda sucks atm because for 
        other collision operators (TRT, MRT) the Lattice2D Structure has to be changed
        -- but it works atm so idgaf """
        feq = self.compute_feq(descriptor, phi, u, mask)
        delta_f = -(1.0 / self.tau) * (f_lattice - feq)
        return delta_f