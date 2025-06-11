import numpy as np
from numba import njit, prange

@njit(parallel=True)
def bounceback(f, opp, mask):
    nz, ny, nx, Q = f.shape
    for i in prange(Q):
        for z in range(nz):
            for y in range(ny):
                for x in range(nx):
                    if mask[z, y, x]:
                        f[z, y, x, i] = f[z, y, x, opp[i]]


@njit(parallel=True)
def apply_velocity_bc(f, u, rho, mask, e, w, u0):
    nz, ny, nx = mask.shape
    Q = e.shape[0]

    for z in prange(nz):
        for y in range(ny):
            for x in range(nx):
                if mask[z, y, x]:
                    rho[z, y, x] = 1.0
                    for d in range(3):
                        u[z, y, x, d] = u0[d]

                    u2 = u0[0]**2 + u0[1]**2 + u0[2]**2
                    for i in range(Q):
                        cu = 0.0
                        for d in range(3):
                            cu += u0[d] * e[i, d]
                        feq = w[i] * rho[z, y, x] * (1 + 3*cu + 4.5*cu**2 - 1.5*u2)
                        f[z, y, x, i] = feq


@njit(parallel=True)
def apply_pressure_bc(f, u, rho, mask, e, w, rho0):
    nz, ny, nx = mask.shape
    Q = e.shape[0]

    for z in prange(nz):
        for y in range(ny):
            for x in range(nx):
                if mask[z, y, x]:
                    rho[z, y, x] = rho0
                    for d in range(3):
                        u[z, y, x, d] = 0.0

                    for i in range(Q):
                        feq = w[i] * rho0
                        f[z, y, x, i] = feq
class Operator:
    def apply(self, f, u=None, rho=None):
        pass

class BounceBack2D(Operator):
    def __init__(self, descriptor, mask):
        
        self.opp = descriptor.opp
        self.mask = mask

    def apply(self, f, u=None, rho=None):
        mask = np.transpose(self.mask) if self.mask.shape != f.shape[:2] else self.mask
        for i in range(len(self.opp)):
            f[mask, i] = f[mask, self.opp[i]]

class BounceBack3D(Operator):
    def __init__(self, descriptor, mask):
        
        self.opp = descriptor.opp
        self.mask = mask

    def apply(self, f, u=None, rho=None):
        self.mask = np.transpose(self.mask, (2, 1, 0)) if self.mask.shape != f.shape[:3] else self.mask
        bounceback(f, self.opp, self.mask)



class VelocityDirichlet2D(Operator):
    def __init__(self, descriptor, collisionOperator, mask, velocity_func):
        self.descriptor = descriptor
        self.collisionOperator = collisionOperator
        self.e = descriptor.e
        self.w = descriptor.w
        self.Q = descriptor.Q
        self.mask = mask
        self.velocity_func = velocity_func

    def apply(self, f, u, rho):
        mask = np.transpose(self.mask) if self.mask.shape != f.shape[:2] else self.mask
        u[mask] = self.velocity_func(u[mask].shape)
        rho[mask] = 1.0  # assume constant pressure
        u2 = u[:,:,0]**2 + u[:,:,1]**2
        feq = self.collisionOperator.compute_feq(self.descriptor, rho, u, mask)
        f[mask] = feq[mask]

class VelocityDirichlet3D(Operator):
    def __init__(self, descriptor, collisionOperator, mask, velocity_func):
        self.e = descriptor.e
        self.w = descriptor.w
        self.mask = mask
        self.velocity_func = velocity_func  # returns constant or spatial field

    def apply(self, f, u, rho):
        self.mask = np.transpose(self.mask, (2, 1, 0)) if self.mask.shape != f.shape[:3] else self.mask
        target_shape = (*self.mask.shape, 3)
        u0 = self.velocity_func(target_shape)[0, 0, 0]  # assumes constant for now
        apply_velocity_bc(f, u, rho, self.mask, self.e, self.w, u0)


class PressureDirichlet3D(Operator):
    def __init__(self, descriptor, collisionOperator, mask, rho_value):
        self.e = descriptor.e
        self.w = descriptor.w
        self.mask = mask
        self.rho_value = rho_value

    def apply(self, f, u, rho):
        self.mask = np.transpose(self.mask, (2, 1, 0)) if self.mask.shape != f.shape[:3] else self.mask
        apply_pressure_bc(f, u, rho, self.mask, self.e, self.w, self.rho_value)



class PressureDirichlet2D(Operator):
    def __init__(self, descriptor, collisionOperator, mask, rho_value):
        self.descriptor = descriptor
        self.collisionOperator = collisionOperator
        self.e = descriptor.e
        self.w = descriptor.w
        self.Q = descriptor.Q
        self.mask = mask
        self.rho_value = rho_value

    def apply(self, f, u, rho):
        mask = np.transpose(self.mask) if self.mask.shape != f.shape[:2] else self.mask
        rho[mask] = self.rho_value
        u[mask, :] = 0.0
        u2 = u[:,:,0]**2 + u[:,:,1]**2
        feq = self.collisionOperator.compute_feq(self.descriptor, rho, u, mask)
        f[mask] = feq[mask]


#########

class PulsedConcentrationDirichlet:
    def __init__(self, descriptor, mask, base_value, pulse_value, t_start=0, t_end=None, sharpness=10.0):
        self.e = descriptor.e
        self.w = descriptor.w
        self.Q = descriptor.Q
        self.mask = mask
        self.base_value = base_value
        self.pulse_value = pulse_value
        self.t = 0
        self.t_start = t_start
        self.t_end = t_end
        self.sharpness = sharpness

    def apply(self, g, u, phi):
        mask = np.transpose(self.mask) if self.mask.shape != phi.shape else self.mask

        # Smooth pulse using tanh
        if self.t_start <= self.t <= self.t_end:
            value = self.pulse_value
        else:
            value = self.base_value

        phi[mask] = value

        for i in range(self.Q):
            cu = u[:, :, 0]*self.e[i, 0] + u[:, :, 1]*self.e[i, 1]
            geq = self.w[i] * phi * (1 + 3*cu)
            g[mask, i] = geq[mask]
        #print(str(self.t) +" ww "+str(self.t_end)+ " ww " +str(value))
        self.t += 1

class ConstantScalarDirichlet:
    def __init__(self, descriptor, mask, value):
        self.e = descriptor.e
        self.w = descriptor.w
        self.Q = descriptor.Q
        self.mask = mask
        self.value = value

    def apply(self, g, u, phi):
        mask = self.mask
        if mask.shape != phi.shape:
            mask = mask.T

        # Set scalar field
        phi[mask] = self.value

        # Compute equilibrium and update distribution
        for i in range(self.Q):
            cu = u[:, :, 0] * self.e[i, 0] + u[:, :, 1] * self.e[i, 1]
            geq = self.w[i] * phi * (1 + 3 * cu)
            g[mask, i] = geq[mask]


class ZeroGradientOutlet:
    def __init__(self, descriptor, mask):
        self.e = descriptor.e
        self.w = descriptor.w
        self.Q = descriptor.Q
        self.mask = mask

    def apply(self, g, u, phi):
        # Neumann boundary -> Boundary cell has to be the same value as the neighboring cell (is this cheating? Idk but it seems to work)
        mask = np.transpose(self.mask) if self.mask.shape != phi.shape else self.mask

        # Buffer our scalar value to apply our neighbors phi
        shifted_phi = np.roll(phi, shift=-1, axis=1)
        phi[mask] = shifted_phi[mask]

        # Compute equilibrium and update distribution
        for i in range(self.Q):
            cu = u[:, :, 0]*self.e[i, 0] + u[:, :, 1]*self.e[i, 1]
            geq = self.w[i] * phi * (1 + 3*cu)
            g[mask, i] = geq[mask]