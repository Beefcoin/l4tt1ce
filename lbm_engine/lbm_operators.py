import numpy as np

class Operator:
    def apply(self, f, u=None, rho=None):
        pass

class BounceBack(Operator):
    def __init__(self, descriptor, mask):
        self.opp = descriptor.opp
        self.mask = mask

    def apply(self, f, u=None, rho=None):
        mask = np.transpose(self.mask) if self.mask.shape != f.shape[:2] else self.mask
        for i in range(len(self.opp)):
            f[mask, i] = f[mask, self.opp[i]]

class VelocityDirichlet(Operator):
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

class PressureDirichlet(Operator):
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