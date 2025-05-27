import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps

# Möbius strip class definition
class MobiusStrip:
    def __init__(self, R=1.0, w=0.3, n=200):
        """
        Initialize Mobius strip parameters.
        :param R: Radius (distance from center to strip centerline)
        :param w: Width of the strip
        :param n: Resolution (number of mesh points)
        """
        self.R = R
        self.w = w
        self.n = n
        self.u = np.linspace(0, 2 * np.pi, n)
        self.v = np.linspace(-w / 2, w / 2, n)
        self.U, self.V = np.meshgrid(self.u, self.v)
        self.x, self.y, self.z = self._generate_mesh()

    def _generate_mesh(self):
        """
        Generate x, y, z coordinates of the Möbius strip.
        """
        U, V = self.U, self.V
        x = (self.R + V * np.cos(U / 2)) * np.cos(U)
        y = (self.R + V * np.cos(U / 2)) * np.sin(U)
        z = V * np.sin(U / 2)
        return x, y, z

    def compute_surface_area(self):
        """
        Approximate surface area numerically using vector calculus and integration.
        """
        du = 2 * np.pi / (self.n - 1)
        dv = self.w / (self.n - 1)

        xu = np.gradient(self.x, du, axis=1)
        xv = np.gradient(self.x, dv, axis=0)
        yu = np.gradient(self.y, du, axis=1)
        yv = np.gradient(self.y, dv, axis=0)
        zu = np.gradient(self.z, du, axis=1)
        zv = np.gradient(self.z, dv, axis=0)

        cross_x = yu * zv - zu * yv
        cross_y = zu * xv - xu * zv
        cross_z = xu * yv - yu * xv

        dA = np.sqrt(cross_x**2 + cross_y**2 + cross_z**2)
        return simps(simps(dA, dx=dv), dx=du)

    def compute_edge_length(self):
        """
        Approximate edge length by summing distances along v = ±w/2.
        """
        edge1 = np.array([self.x[0], self.y[0], self.z[0]]).T
        edge2 = np.array([self.x[-1], self.y[-1], self.z[-1]]).T

        def arc_length(edge):
            diff = np.diff(edge, axis=0)
            return np.sum(np.linalg.norm(diff, axis=1))

        return arc_length(edge1) + arc_length(edge2)

    def plot(self):
        """
        Plot the Möbius strip in 3D.
        """
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(self.x, self.y, self.z, cmap='viridis', edgecolor='none')
        ax.set_title('Möbius Strip Visualization')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.tight_layout()
        plt.savefig("mobius_strip_plot.png")
        plt.show()


# Run and demonstrate
if __name__ == "__main__":
    mobius = MobiusStrip(R=1.0, w=0.3, n=200)
    area = mobius.compute_surface_area()
    length = mobius.compute_edge_length()
    print(f"Surface Area ≈ {area:.3f}")
    print(f"Edge Length ≈ {length:.3f}")
    mobius.plot()
