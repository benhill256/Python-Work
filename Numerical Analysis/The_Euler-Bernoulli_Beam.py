import numpy as np
from scipy import sparse
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt


# Setup.

# Dimensions of the board in meters.
w = 30*0.01
d = 3*0.01
L = 2

# These are constants defined in the setup.
E = 1.3 * 10**10
I = (w * d**3) / 12
g = 9.81

# This function constructs the sparse matrix A.
def make_matrix(a, n):
    # Initiate the matrix as an nxn zero matrix.
    matrix = np.zeros((n,n))
    for i in range(7):
        # Add matrices to the zero matrix, each with a diagonal in the proper place.
        matrix += np.diag(a[i], k = -i+3) # k=-i+3 is to convert i into the index of the proper diagonal.
    
    # Return the matrix as a sparse matrix.
    return sparse.csr_matrix(matrix)

# This function creates the diagonals used in constructing the clamped-free coefficient matrix A.
def get_a(n):
    diags = [[-1/4], [8/3], [-9], [16], [-4 for i in range(n-3)], [1 for i in range(n-4)], [0 for i in range(n-4)]]

    for i in range(n-4):
        diags[0].append(0)

    for i in range(n-3):
        diags[1].append(1)
        diags[2].append(-4)
        diags[3].append(6)
    
    diags[2].append(-28/17)
    diags[3].extend((72/17, 72/17))
    diags[4].extend((-60/17, -156/17))
    diags[5].extend((16/17, 96/17))
    diags[6].append(-12/17)

    diags_array = np.array(diags, dtype=object)
    
    return diags_array

# This is the function f(x), given in the problem setup.
def f(x):
    return -480 * w * d * g

# This defines the scalar h, the length of the subintervals.
def h(n):
    return L / n

# This creates the b vector, composed only of values of f(x). Note that f(x) is constant, so x does not affect the value.
def make_b(n, f):
    x = np.linspace(L/n, L, n)
    b = []
    for i in x:
        b.append(f(i))

    return np.array(b)









# Part 1
print("Part 1:")

# y is the solution vector in this equation, Ay = b.
y = spla.spsolve(make_matrix(get_a(10), 10), (((h(10))**4) / (E * I)) * make_b(10, f))
print(f"The solution to this system for n=10 is y = {y}")

print(input("\nPress enter to continue to Part 2"))










# Part 2
print("\nPart 2: See the generated plot for a visualization of errors in the solution y.")

# The board is divided into n = 10 subintervals in this part.
n_2 = 10
x = np.linspace(L/n_2, L, n_2)

# This is the y vector that is the "true solution", as given in the instructions.
def true_soln2(x):
    return (f(x) / (24 * E * I)) * x**2 * (x**2 - 4*L*x + 6*L**2)

# Plot the solution y from step 1 against the "true solution".
fig = plt.figure()
plt.plot(x, y, "ro", label=r"Approximate y")
plt.plot(x, true_soln2(x), "-b", label=r"True y")
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.legend(loc = 'lower left')
plt.show()

print(f"Error at x=L: {np.format_float_scientific(abs(true_soln2(L) - y[-1]), unique=True, precision=2)}")

print(input("\nPress enter to continue to Part 3"))










# Part 3
print("\nPart 3:")

# This function constructs the matrix A as a Numpy array, which we will need for computing the condition number.
def make_matrix_np(a, n):
    # Initiate the matrix as an nxn zero matrix.
    matrix = np.zeros((n,n))
    for i in range(7):
        # Add matrices to the zero matrix, each with a diagonal in the proper place.
        matrix += np.diag(a[i], k = -i+3) # k=-i+3 is to convert i into the index of the proper diagonal.
    
    # Return the matrix as a sparse matrix.
    return matrix

# We will be using integers k in [1, 11] to define values of n.
k = [i for i in range(1, 12)]
n_k = [10*2**i for i in k]

# Initiate the list of errors for the k-values.
errors = []

# Find y vector solutions to Ay = b, for each k-value.
for j in n_k:
    # Find the solution y.
    y_k = spla.spsolve(make_matrix(get_a(j), j), (((h(j))**4) / (E * I)) * make_b(j, f))

    # Difference between the last entries of this y vector and the "true solution"; this is the error at x = L.
    error = abs(true_soln2(L) - y_k[-1])
    errors.append(error)

    # Generate a table of the errors of the vectors and condition numbers of the matrices.
    print(f"For n={j}, the error at x=L is {np.format_float_scientific(error, unique=True, precision=2)} and the condition number of A is {np.format_float_scientific(np.linalg.cond(make_matrix_np(get_a(j), j), np.inf), unique=True, precision=2)}")

print("\nThe error in the approximation increases with n from the start because the condition number of A enormously increases with n.")

print(input("\nPress enter to continue to Part 4"))










# Part 4
print("\nPart 4:")

from sympy import *

x = symbols("x")
p = symbols("p")

y = (((f(x) / (24 * E * I)) * x**2 * (x**2 - 4*L*x + 6*4)) - (p*g*L/(E*I*pi)) * ((L**3/(pi)**3)*sin(pi*x/L) - (x**3)/6 + (L/2)*x**2 - (L**2/(pi**2)*x)))
yprime1 = diff(y, x)
yprime2 = diff(diff(y, x), x)
yprime3 = diff(diff(diff(y, x), x), x)
EIyprime4 = diff(diff(diff(diff(E*I*y, x), x), x), x)

if EIyprime4 := -9.81*p*sin(pi*x/L) - 42.3792:
    print("The solution y(x) given satisfies the Euler-Bernoulli equation.")

if y.subs(x, 0) == 0 and yprime1.subs(x, 0) == 0 and yprime2.subs(x, L) == 0 and yprime3.subs(x, L) == 0:
    print("It also satisfies the clamped-free boundary conditions.")

print(input("\nPress enter to continue to Part 5"))










# Part 5
print("\nPart 5:")

# For some reason, these constants were simply not loading properly into the functions here, causing persistent crashes. So here they are defined again.
E = 1.3 * 10**10
I = (w * d**3) / 12

# p is given to be 100 kg/m for this part.
p = 100

# The new function f(x) with the sinusoidal pile added onto it.
def f5(x):
    return f(x) - g*p*np.sin(np.pi*x/L)

def true_soln5(x, E, I):
    return (f(x) / (24 * E * I)) * x**2 * (x**2 - 4*L*x + 6*L**2) - (p*g*L/(E*I*np.pi)) * ((L**3/(np.pi)**3)*np.sin(np.pi*x/L) - (x**3)/6 + (L/2)*x**2 - (L**2/(np.pi**2)*x))

# We will be using values of k and n as they were defined in Part 3.

# Initiate the list of errors for the k-values.
errors = []
h_values = []

# Find y vector solutions to Ay = b, for each k-value.
for j in n_k:
    # Find the solution y.
    y_k = spla.spsolve(make_matrix(get_a(j), j), float((((h(j))**4) / (E * I))) * make_b(j, f5))

    # Difference between the last entries of this y vector and the "true solution"; this is the error at x = L.
    error = abs(true_soln5(L, E, I) - y_k[-1])
    errors.append(error)

    h_values.append(h(j))

    # Generate a table of the errors of the vectors and condition numbers of the matrices.
    print(f"For n={j}, the error at x=L is {np.format_float_scientific(error, unique=True, precision=2)} and the condition number of A is {np.format_float_scientific(np.linalg.cond(make_matrix_np(get_a(j), j), np.inf), unique=True, precision=2)}")

print("\nThe error in the approximation increases with n after n=640 because the condition number of A enormously increases with n.")
print("The error at x=L is proportional to h^2 for about n<=640. For such n, the slope on the log-log plot seems to be precisely 2.")
print("For n>640, the error seems to jump around erratically, but still generally increases")
print("See the generated log-log plot of the errors and h-values for this relationship.")
print("Note that the associated values of n get larger from right to left in the plot, because h gets smaller as n grows.")

y_values = [np.log10(i) for i in errors]

# Plot the errors of the y solutions at x=L against h.
fig = plt.figure()
# plt.plot(np.log10(h_values), np.log10(errors), "-b")
plt.plot(np.log10(h_values), y_values, "-b")
plt.xlabel(r'$Log(h)$')
plt.ylabel(r'$Log(errors)$')
plt.show()

print(input("\nPress enter to continue to Part 6"))










# Part 6
print("\nPart 6:")

def f6(x):
    return f(x) - g*70/0.2

# This is the optimal n found in the previous part.
n_6 = 640
x = np.linspace(L/n_6, L, n_6)

# The other function to make the b vector was not working with the f6 function, because f6 has to be defined piecewise. I have moved the "pieces" into here.
def make_b6(n):
    x = np.linspace(L/n, L, n)
    b = []
    for i in x:
        if i > 1.8 and i < L:
            value = (((h(n_6))**4) / (E * I)) * f6(i)
        else:
            value = (((h(n_6))**4) / (E * I)) * f(i)
        b.append(value)

    return np.array(b)

y = spla.spsolve(make_matrix(get_a(n_6), n_6), make_b6(n_6))
print(f"Under the weight of a diver, the deflection of the board at x=L is {np.format_float_scientific(y[-1], unique=True, precision=5)}")

# Plot the x inputs and resulting y values.
fig = plt.figure()
plt.plot(x, y, "-b")
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.show()

print(input("\nPress enter to continue to Part 7"))










# Part 7
print("\nPart 7:")

# This function creates the diagonals used in constructing the clamped-clamped coefficient matrix A.
def get_a_clamped(n):
    diags = [[-1/4], [8/3], [-9], [16], [-4 for i in range(n-3)], [1 for i in range(n-4)], [0 for i in range(n-4)]]

    for i in range(n-4):
        diags[0].append(0)

    for i in range(n-3):
        diags[1].append(1)
        diags[2].append(-4)
        diags[3].append(6)
    
    diags[2].append(-4)
    diags[3].extend((6, 16))
    diags[4].extend((-4, -9))
    diags[5].extend((1, 8/3))
    diags[6].append(-1/4)

    diags_array = np.array(diags, dtype=object)
    
    return diags_array

# The number of entries in the b vector is now n+1 instead of n, starting at 0, so this is the amended function.
def make_b7(n, f):
    x = np.linspace(0, L, n)
    b = []
    for i in x:
        b.append(f(i))

    return np.array(b)

# Note that f5 is the function with the sinusoidal pile, so I see no need to remake it here.
def true_soln7(x):
    # return (f5(x) / (24 * E * I)) * x**2 * (L - x)**2 - (p*g*(L**2)/(E*I*(np.pi)**4)) * (L**2*np.sin(np.pi*x/L) + np.pi*x*(x - L))
    return (f(x) / (24 * E * I)) * x**2 * (L - x)**2 - (p*g*(L**2)/(E*I*(np.pi)**4)) * (L**2*np.sin(np.pi*x/L) + np.pi*x*(x - L))

# We will be using values of k as they were defined in Part 3; these are the new n values, in agreement with there being n+1 x values.
n_k = [10*2**i + 1 for i in k]

# Initiate the list of errors for the k-values.
errors = []
h_values = []

# Find y vector solutions to Ay = b, for each k-value.
for j in n_k:
    # Find the solution y.
    y_k = spla.spsolve(make_matrix(get_a_clamped(j), j), (((h(j))**4) / (E * I)) * make_b7(j, f5))

    # Difference between the last entries of this y vector and the "true solution"; this is the error at x = L.
    error = abs(true_soln7(L/2) - y_k[int((j+1)/2)]) # The (j+1)/2 is the index of the x value in the center of the beam, x=L/2
    errors.append(error)

    h_values.append(h(j))

    # Generate a table of the errors of the vectors and condition numbers of the matrices.
    print(f"For n={j}, the error at x=L/2 is {np.format_float_scientific(error, unique=True, precision=2)} and the condition number of A is {np.format_float_scientific(np.linalg.cond(make_matrix_np(get_a_clamped(j), j), np.inf), unique=True, precision=2)}")

print("\nThe error in the approximation at x=L/2 decreases with n, even though the condition number of A increases with n.")
print("It seems to be proportional to h, not h^2.")
print("See the generated log-log plot of the errors and h-values for this relationship.")

# Plot the errors of the y solutions at x=L against h.
fig = plt.figure()
plt.plot(np.log10(h_values), np.log10(errors), "-b")
plt.xlabel(r'$Log(h)$')
plt.ylabel(r'$Log(errors)$')
plt.show()
