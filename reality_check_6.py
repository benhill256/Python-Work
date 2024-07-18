# tacoma.py
# implements the code for RC6
# Inputs: inter = time interval, ic = [y, y', theta, theta'],
# n = number of steps, p = steps per point plotted
# Calls a one-step method such as trapstep or RK4step

import numpy as np
import matplotlib.pylab as plt


d = 0.01
omega = 2*np.pi*38/60

def tacoma(inter, ic, n, p, windspeed, method):
    W = windspeed

    # trapezpoid method step
    def trapstep(t, w, h, f):
        return w + (h/2)*(f(t, w, W) + f(t + h, w + h*f(t, w, W), W))
        
    # Runge Kutta step
    def RK4step(t, w, h):
        s1 = ydot(t, w, W)
        s2 = ydot(t + h/2, w + (h/2)*s1, W)
        s3 = ydot(t + h/2, w + (h/2)*s2, W)
        s4 = ydot(t + h, w + h*s3, W)
        return w + (h/6)*(s1 + 2*s2 + 2*s3 + s4)

    # this is f(t, y) in y' = f(t, y)
    def ydot(t, y, W):
        K = 1000
        m = 2500
        L = 6
        a = 0.2
        c1 = K/(m*a)
        a1 = np.exp(a*(y[0] - L*np.sin(y[2])))
        a2 = np.exp(a*(y[0] + L*np.sin(y[2])))
        ydot = np.zeros(4)
        ydot[0] = y[1]
        ydot[1] = -d*y[1] - c1*(a1+a2-2) + 0.2*W*np.sin(omega*t)
        ydot[2] = y[3]
        ydot[3] = -d*y[3] + c1*(3/L)*np.cos(y[2])*(a1-a2)
        return ydot


    # use n points
    h = (inter[1]-inter[0])/float(n)

    # build t vectors
    t = np.zeros((p+1))
    t[0] = inter[0]
    rt = np.zeros(int(n/p)+1)
    rt[0] = inter[0]

    # build y vector
    # y[0] is vertical displacement
    # y[1] is the derivative of y[0]
    # y[2] is the rotation of the roadbed
    # y[3] is the derivative of y[2]
    y = np.zeros((p+1,4))
    y[0,:] = ic
    ry = np.zeros((int(n/p)+1,4))
    ry[0] = ic

    # compute solution
    for k in np.arange(1, int(n/p)+1):
        for i in np.arange(p):
            t[i+1] = t[i] + h
            if method == "Trap":
                y[i+1,:] = trapstep(t[i], y[i,:], h, ydot)
            elif method == "RK4":
                y[i+1,:] = RK4step(t[i], y[i,:], h)
        rt[k] = t[i+1]
        ry[k,:] = y[i+1,:]
        t[0] = rt[k]
        y[0:k] = ry[k,:]

    # return time and solution
    return rt, ry

a = 0
b = 1000
n = 50000
p = 4
theta0 = 10.0**(-3)










# Part 1
print("\nPart 1\n")

ot1, oy1 = tacoma([a, b], [0, 0, theta0, 0], n, p, 80, "Trap") # W = 80 km/h

# oy[0] = y, oy[2] = theta
fig = plt.figure(1, figsize = (8, 6))
plt.plot(ot1, oy1[:,0], 'b-')
plt.xlabel('$t$', fontsize = 18, color = 'Blue')
plt.ylabel('$y(t)$', fontsize = 18, color = 'Blue')
plt.title('Vertical Displacement, Trap')

fig = plt.figure(2, figsize = (8, 6))
plt.plot(ot1, oy1[:,2], 'b-')
plt.xlabel('$t$', fontsize = 18, color = 'Blue')
plt.ylabel('$\\theta(t)$', fontsize = 18, color = 'Blue')
plt.title('Torsional Displacement, Trap')

plt.show()
print(input("Press enter to proceed to the next part."))










# Part 2
print("\nPart 2\n")

ot2, oy2 = tacoma([a, b], [0, 0, theta0, 0], n, p, 80, "RK4") # W = 80 km/h

# oy[0] = y, oy[2] = theta
fig = plt.figure(1, figsize = (8, 6))
plt.plot(ot2, oy2[:,0], 'b-')
plt.xlabel('$t$', fontsize = 18, color = 'Blue')
plt.ylabel('$y(t)$', fontsize = 18, color = 'Blue')
plt.title('Vertical Displacement, RK4')

fig = plt.figure(2, figsize = (8, 6))
plt.plot(ot2, oy2[:,2], 'b-')
plt.xlabel('$t$', fontsize = 18, color = 'Blue')
plt.ylabel('$\\theta(t)$', fontsize = 18, color = 'Blue')
plt.title('Torsional Displacement, RK4')

plt.show()
print(input("Press enter to proceed to the next part."))










# Part 3
print("\nPart 3\n")

theta0 = [10.0**(-3), 10.0**(-4), 10.0**(-5), 10.0**(-6)]

mag_factors = []

for theta in theta0:
    ot3, oy3 = tacoma([a, b], [0, 0, theta, 0], n, p, 50, "RK4") # W = 50 km/h
    mag_factors.append(max(oy3[:,2]) / theta)

print(f"The magnification factors for the respective starting angles are {mag_factors}")

print(input("\nPress enter to proceed to the next part."))










# Part 4
print("\nPart 4\n")

theta0 = 10.0**(-3)

ot4, oy4 = tacoma([a, b], [0, 0, theta0, 0], n, p, 58.993, "RK4") # W = 58.993 km/h, found by trial and error.
print(f"Magnification Factor for W=58.993: {max(oy4[:,2]) / theta0}")

# Checking whether the magnification factor is consistent.
theta0 = [10.0**(-3), 10.0**(-4), 10.0**(-5), 10.0**(-6)]

mag_factors = []

for theta in theta0:
    ot4, oy4 = tacoma([a, b], [0, 0, theta, 0], n, p, 58.993, "RK4") # W = 58.993 km/h
    mag_factors.append(max(oy4[:,2]) / theta)

print(f"The magnification factors for the respective starting angles are {mag_factors}")

print(input("\nPress enter to proceed to the next part."))










# Part 5
print("\nPart 5\n")

theta0 = 10.0**(-3)

tolerance = 0.5e-3

def magnification_factor(W):
    _, oy4 = tacoma([a, b], [0, 0, theta0, 0], n, p, W, "RK4")
    return max(oy4[:,2]) / theta0

# Wind speed lower bound and wind speed upper bound, for a given globally-defined initial theta and tolerance.
def find_minimum_wind_speed(lb, ub):
    # Bisection method loop
    while ub - lb > tolerance:
        mid = (lb + ub) / 2
        factor = magnification_factor(mid)
        if factor > 100:
            ub = mid
        else:
            lb = mid

    return (lb + ub) / 2

min_wind_speed = find_minimum_wind_speed(57, 60)
print(f"Minimum wind speed is approximately {min_wind_speed} km/h, by bisection.")

print(input("\nPress enter to proceed to the next part."))










# Part 6
print("\nPart 6\n")

# Modify the damping coefficient and angular frequency
d = 0.02
omega = 3.0

min_wind_speed = find_minimum_wind_speed(25, 30)
print(f"Minimum wind speed within tolerance is approximately {min_wind_speed} km/h, by bisection.")

# Plot the displacements as before.
ot6, oy6 = tacoma([a, b], [0, 0, theta0, 0], n, p, min_wind_speed, "RK4")

fig = plt.figure(1, figsize = (8, 6))
plt.plot(ot6, oy6[:,0], 'b-')
plt.xlabel('$t$', fontsize = 18, color = 'Blue')
plt.ylabel('$y(t)$', fontsize = 18, color = 'Blue')
plt.title('Vertical Displacement, RK4')

fig = plt.figure(2, figsize = (8, 6))
plt.plot(ot6, oy6[:,2], 'b-')
plt.xlabel('$t$', fontsize = 18, color = 'Blue')
plt.ylabel('$\\theta(t)$', fontsize = 18, color = 'Blue')
plt.title('Torsional Displacement, RK4')

plt.show()