import numpy as np
from scipy.optimize import minimize


# Part 1
print("\nPart 1")

def add_consts(x):
    x_con = np.append(np.zeros(5), np.insert(x,1,0))
    return x_con

def remove_consts(x):
    x_con = np.delete(x, np.append(range(5),6))
    return x_con

def difference(p1, p2):
    r = np.linalg.norm(p1 - p2)
    return r**6, r**12

def potential(coords):
    p = np.reshape(add_consts(coords),((-1,3)))
    n = len(p)
    energy = 0

    for i in range(n):
        for j in range(i+1, n):
            r6, r12 = difference(p[i], p[j])
            energy += 1/r12 - 2/r6
    
    return energy

# Applying Nelder-Mead optimization to find minimum energy for n=5
initial_guesses = [
    [0,0,0,1,1,1,2,2,2], [2,2,2,2,3,2,3,3,1], [1,1,1,1,1,1,1,0,1]
    # Add more initial guesses if desired
]

min_energies = []
for guess in initial_guesses:
    result = minimize(potential, guess, method="Nelder-Mead", options={"disp": True})
    if result.fun < 0: # This weeds out one of the initial guesses which is yielding a nan value.
        min_energies.append(result.fun)

print(f"Test Run: {potential([-1,-1,0,-1,0,0,-1,0,-1])}")

print("Minimum energy for n=5:", min(min_energies))
print("Number of steps required:", result.nit)










# Part 2
print("\nPart 2")
print("See generated plot.")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plot_atoms(atom_coords):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(atom_coords[:, 0], atom_coords[:, 1], atom_coords[:, 2], color='blue', zorder=2)

    # Connect the atoms with line segments to form triangles
    for i in range(len(atom_coords)):
        for j in range(i+1, len(atom_coords)):
            ax.plot([atom_coords[i, 0], atom_coords[j, 0]], [atom_coords[i, 1], atom_coords[j, 1]], [atom_coords[i, 2], atom_coords[j, 2]], color='black', zorder=1)

    # Label the points
    for i, (x, y, z) in enumerate(atom_coords):
        ax.text(x, y, z, f"Atom {i+1}", color="green")

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f"Minimum Energy Configuration of {len(atom_coords)} Atoms")
    plt.show()

min_atoms = np.reshape(add_consts(result.x),((-1,3)))
plot_atoms(min_atoms)










# Part 3
print("\nPart 3")
def gradient(points):
    p = np.reshape(add_consts(points),((-1,3)))
    n = len(p)

    partx = np.zeros((n,n))
    party = np.zeros((n,n))
    partz = np.zeros((n,n))

    for i in range(n-1):
        for j in range(i+1,n):
            v = p[j] - p[i]
            rij = np.dot(v,v)
            partx[i,j] = 12*((rij**(-7) - rij**(-4)))*v[0]
            party[i,j] = 12*((rij**(-7) - rij**(-4)))*v[1]
            partz[i,j] = 12*((rij**(-7) - rij**(-4)))*v[2]

    du = np.zeros(n*3)
    for i in range(n):
        du[i*3] = np.sum(partx[i,i:]) - np.sum(partx[:i,i])
        du[i*3+1] = np.sum(party[i,i:]) - np.sum(party[:i,i])
        du[i*3+2] = np.sum(partz[i,i:]) - np.sum(partz[:i,i])
    # print(du)

    return remove_consts(du)

# Applying Python minimization function using the gradient for n=5
result_gradient = minimize(potential, guess, jac=gradient, method='CG')
print("Minimum energy using gradient for n=5:", result_gradient.fun)
print(f"Gradient at test coordinates: {gradient([1,1,1,1,1,1,1,0,1])}")










# Part 4
print("\nPart 4")
# Using a different method, without using the gradient
result_global_min = minimize(potential, guess, method='Powell')
print("Global minimum energy for n=5:", result_global_min.fun)










# Part 5
print("\nPart 5")
guess = [1,1,1,1,1,1,1,0,1,2,1,2]
result_NM = minimize(potential, guess, method="Nelder-Mead", options={"disp": True})
print("Minimum energy using Nelder-Mead for n=6:", result_NM.fun)
print("Number of steps required:", result_NM.nit)
result_gradient = minimize(potential, guess, jac=gradient, method='CG')
print("Minimum energy using gradient for n=6:", result_gradient.fun)
print("Number of steps required:", result_gradient.nit)
result_global_min = minimize(potential, guess, method='Powell')
print("Global minimum energy for n=5:", result_global_min.fun)

print("""It seems that the Nelder-Mead method is much slower in its convergence
than the conjugate gradient method, and the conjugate gradient method
also gives us the global minimum instead of just an approximation of it
through a local minimum. Because of this, I would trust the conjugate
gradient method more.""")










# Part 6
print("\nPart 6")
print("See generated plot.")
min_atoms = np.reshape(add_consts(result_gradient.x),((-1,3)))
plot_atoms(min_atoms)










# Part 7
print("\nPart 7")
guess = [1,1,1,1,1,1,1,0,1,2,1,2,3,1,2,2,3,2,1,3,3,3,3,3] # n = 10 atoms
result_gradient = minimize(potential, guess, jac=gradient, method='CG')
print("Minimum energy using gradient for n=10:", result_gradient.fun)
print("Number of steps required:", result_gradient.nit)

min_atoms = np.reshape(add_consts(result_gradient.x),((-1,3)))
plot_atoms(min_atoms)
