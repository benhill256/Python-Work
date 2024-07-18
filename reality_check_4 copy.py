import numpy as np
from scipy.optimize import root
import sympy as sp

# The speed of light in km/s
c = 299792.458

def equations(vars):
    x, y, z, d = vars
    eq1 = np.sqrt((x - A1)**2 + (y - B1)**2 + (z - C1)**2) - c * (t1 - d)
    eq2 = np.sqrt((x - A2)**2 + (y - B2)**2 + (z - C2)**2) - c * (t2 - d)
    eq3 = np.sqrt((x - A3)**2 + (y - B3)**2 + (z - C3)**2) - c * (t3 - d)
    eq4 = np.sqrt((x - A4)**2 + (y - B4)**2 + (z - C4)**2) - c * (t4 - d)

    return [eq1, eq2, eq3, eq4]










# Part 1
print("\nPart 1")

# Given positions and time intervals for the four satellites.
A1, B1, C1, t1 = 15600, 7540, 20140, 0.07074
A2, B2, C2, t2 = 18760, 2750, 18610, 0.07220
A3, B3, C3, t3 = 17610, 14630, 13480, 0.07690
A4, B4, C4, t4 = 19170, 610, 18390, 0.07242

# This is the given initial vector for the iterative approximation.
guess = np.array([0, 0, 6370, 0])

# This should be approximately the correct solution.
result = root(equations, guess)

print(result.x)










# Part 2
print("\nPart 2")

def solve_quadratic_formula(A, B, C):
    # Solving quadratic equation using the quadratic formula
    discriminant = B**2 - 4*A*C
    x1 = (-B + np.sqrt(discriminant)) / (2*A)
    x2 = (-B - np.sqrt(discriminant)) / (2*A)
    return x1, x2

# This matrix and w vector were derived by hand, subtracting each equation from the first in turn.
matrix = np.array([[2*(A2-A1), 2*(B2-B1), 2*(C2-C1), (2*c**2)*(t1-t2)], [2*(A3-A1), 2*(B3-B1), 2*(C3-C1), (2*c**2)*(t1-t3)], [2*(A4-A1), 2*(B4-B1), 2*(C4-C1), (2*c**2)*(t1-t4)]])
w = np.array([[A2**2-A1**2 + B2**2-B1**2 + C2**2-C1**2 + c**2*(t1**2-t2**2)], [A3**2-A1**2 + B3**2-B1**2 + C3**2-C1**2 + c**2*(t1**2-t3**2)], [A4**2-A1**2 + B4**2-B1**2 + C4**2-C1**2 + c**2*(t1**2-t4**2)]])

# Find the RREF of the augmented matrix.
aug_matrix = sp.Matrix(np.concatenate((matrix, w), axis=1))
print(aug_matrix)
rref_matrix = aug_matrix.rref()[0]
# rref_matrix = sp.Matrix(np.concatenate((matrix, w), axis=1)).rref()[0]

print(rref_matrix)

# These are the resulting equations from the RREF matrix.
# d = sp.Symbol("d")
# x = rref_matrix[4] - rref_matrix[3]*d
# y = rref_matrix[9] - rref_matrix[8]*d
# z = rref_matrix[14] - rref_matrix[13]*d

# Here is the simplified result of the first equation in (4.37), as given by Mathematica.
# quadratic = 4.91198*10**7 + 1.50772*(10**10)*d - 8.28546*(10**10)*d**2

# Show the solution obtained using the quadratic formula, given the coefficients above.
# Note that the root assigned to the variable _ is not the correct intersection, so we will not use it.
d, _ = solve_quadratic_formula(-8.28546*(10**10), 1.50772*10**10, 4.91198*10**7)

# Solving the equations from earlier, given the value of d we found.
x = rref_matrix[4] - rref_matrix[3]*d
y = rref_matrix[9] - rref_matrix[8]*d
z = rref_matrix[14] - rref_matrix[13]*d

print(f"Position: {[x, y, z, d]}")










# Parts 4 and 5 are largely based on code from ChatGPT.
# Part 4
print("\nPart 4")

# Given parameter for Part 4
rho = 26570  # km

# Values for phi and theta are arbitrarily chosen to be at the north pole and equiangular around the equator.
phi_values = [np.pi/2, 0, 0, 0]
theta_values = [0, 0, np.pi/3, np.pi*2/3]
x, y, z, d = 0, 0, 6370, 0.0001

def max_error_magnification(x, y, z, d, satellite_positions, measured_times, delta_t):
    position_errors = []
    error_magnification_factors = []

    for i in range(len(satellite_positions)):
        
        # Since 10^-8 seconds time error corresponds to 3 meters position error, this lists the position errors.
        # delta_xyz = [abs(3*delta_t[i]/10**-8), abs(3*delta_t[i]/10**-8), abs(3*delta_t[i]/10**-8)]
        delta_xyz = [3*delta_t[i]/10**-8, 3*delta_t[i]/10**-8, 3*delta_t[i]/10**-8]
        position_errors.append(max(delta_xyz))

        # Calculating error magnification factor
        # emf = np.linalg.norm(delta_xyz, ord=np.inf) / np.linalg.norm(delta_t, ord=np.inf)
        emf = max(delta_xyz) / (c*max(delta_t))
        error_magnification_factors.append(emf)
    
    return max(error_magnification_factors), max(position_errors)

# Define satellite positions from spherical coordinates
def spherical_to_cartesian(rho, phi, theta):
    x = rho * np.cos(phi) * np.cos(theta)
    y = rho * np.cos(phi) * np.sin(theta)
    z = rho * np.sin(phi)
    return x, y, z

# Calculate satellite positions, distances from the receiver (R), and travel times.
satellite_positions = np.array([spherical_to_cartesian(rho, phi, theta) for phi, theta in zip(phi_values, theta_values)])
R_values = np.linalg.norm(satellite_positions - np.array([x, y, z]), axis=1)
travel_times = d + R_values / c




A1, B1, C1, t1 = 15600, 7540, 20140, 0.07074
A2, B2, C2, t2 = 18760, 2750, 18610, 0.07220
A3, B3, C3, t3 = 17610, 14630, 13480, 0.07690
A4, B4, C4, t4 = 19170, 610, 18390, 0.07242

# This is the initial vector for the iterative approximation.
guess = np.array([0, 0, 6370, 0.0001])

# This should be approximately the correct solution.
result = root(equations, guess)




# Calculate error magnification factor for perturbed times
delta_t_1 = np.array([0.5e-8, 1e-8, -1e-8, -0.5e-8])
max_emf, max_position_error = max_error_magnification(x, y, z, d, satellite_positions, travel_times, delta_t_1)

print(f"Satellite positions: {satellite_positions}")
print(f"Satellite ranges: {R_values}")
print(f"Travel times: {travel_times}")
print(f"Maximum position error: {max_position_error}")
print(f"Maximum EMF1: {max_emf}")

delta_t_2 = np.array([0.5e-8, 0.75e-8, -1e-8, -0.8e-8])
max_emf, _ = max_error_magnification(x, y, z, d, satellite_positions, travel_times, delta_t_2)
print(f"Maximum EMF2: {max_emf}")

delta_t_3 = np.array([0.3e-8, -0.3e-8, -0.5e-8, -1e-8])
max_emf, _ = max_error_magnification(x, y, z, d, satellite_positions, travel_times, delta_t_3)
print(f"Maximum EMF3: {max_emf}")

delta_t_4 = np.array([0.95e-8, 0.95e-8, -0.95e-8, -0.95e-8])
max_emf, _ = max_error_magnification(x, y, z, d, satellite_positions, travel_times, delta_t_4)
print(f"Maximum EMF4: {max_emf}")

delta_t_5 = np.array([0.2e-8, 0.2e-8, -0.2e-8, -0.2e-8])
max_emf, _ = max_error_magnification(x, y, z, d, satellite_positions, travel_times, delta_t_5)
print(f"Maximum EMF5: {max_emf}")










# Part 5
print("\nPart 5")

# Repeat Part 4 with tightly grouped satellites
phi_values_tight = [0.01, 0.03, 0.05, 0.05]
theta_values_tight = [0.0, 0.03, 0.01, 0.05]
delta_t_0 = np.array([0, 0, 0, 0])

# Calculate satellite positions for tightly grouped satellites
satellite_positions_tight = np.array([spherical_to_cartesian(rho, phi, theta) for phi, theta in zip(phi_values_tight, theta_values_tight)])
R_values_tight = np.linalg.norm(satellite_positions - np.array([x, y, z]), axis=1)
travel_times_tight = d + R_values_tight / c

# Calculate error magnification factor for perturbed times for tightly grouped satellites
max_emf_tight, max_position_error_tight = max_error_magnification(x, y, z, d, satellite_positions_tight, travel_times_tight, delta_t_1)

print(f"Satellite positions: {satellite_positions_tight}")
print(f"Satellite ranges: {R_values_tight}")
print(f"Travel times: {travel_times_tight}")
print(f"Maximum position error: {max_position_error_tight}")
print(f"Maximum EMF1: {max_emf_tight}")

max_emf, _ = max_error_magnification(x, y, z, d, satellite_positions_tight, travel_times_tight, delta_t_2)
print(f"Maximum EMF2: {max_emf_tight}")

max_emf, _ = max_error_magnification(x, y, z, d, satellite_positions_tight, travel_times_tight, delta_t_3)
print(f"Maximum EMF3: {max_emf_tight}")

max_emf, _ = max_error_magnification(x, y, z, d, satellite_positions_tight, travel_times_tight, delta_t_4)
print(f"Maximum EMF4: {max_emf_tight}")

max_emf, _ = max_error_magnification(x, y, z, d, satellite_positions_tight, travel_times_tight, delta_t_5)
print(f"Maximum EMF5: {max_emf_tight}")