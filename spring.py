import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

#DEFINE INITIAL CONSTANTS (all in SI)
##################################
g = -9.8 #Assumed constant gravitational field strength
k = 1.0 #Spring constant
totalMass = 1.0
M = 100.0 #Extra mass attached to bottom
N = 5 #Number of particles
simulationTime = 1
timeStep = 0.01
##################################
#INITIAL CONSTANTS DEFINED

#GENERATE CALCULATED CONSTANTS
##################################
m = totalMass / (N**2) #Particle mass
numInstants = int(1 + (simulationTime / timeStep))
beta = m * g / k
##################################
#CALCULATED CONSTANTS GENERATED

#Create a 2D array of the particle positions, where each entry is an array, giving the positions of the particles at an instant in time
blankArray = [0.0] * N
ps = [list(blankArray) for i in range(numInstants)] #Be sure to create actual copies of the array, not copies of the pointer to it
vs = [list(blankArray) for i in range(numInstants)]

#Fill the initial positions
for n in range(N):
  ps[0][n] = (g/k) * ( ((n+1)*N - (n*(n+1)/2))*m + (n+1)*M ) #Zero-indexing, so it'll look a little different than the formula in class


def updatePositions():
  for n in range(N):
    ps[instant][n] = ps[instant-1][n] + vs[instant - 1][n] * timeStep


def updateVelocities():
  fDownOnTop = k*(ps[instant][1] - ps[instant][0]) + m*g
  vs[instant][0] = vs[instant-1][0] + (fDownOnTop / m)*timeStep

  fUpOnBottom = k * (ps[instant][N-1] - ps[instant][N-2])
  fDownOnBottom = (m+M)*g
  accOnBottom = (fDownOnBottom - fUpOnBottom) / m
  vs[instant][N-1] = vs[instant-1][N-1] + accOnBottom*timeStep

  for n in range(1, N-1):
    fUp = k * (ps[instant][n] - ps[instant][n-1])
    fDown = k*(ps[instant][n + 1] - ps[instant][n]) + m*g
    acc = (fDown - fUp) / m
    vs[instant][n] = vs[instant-1][n] + acc*timeStep


#Then we release the top (massless) spring, which instantly snaps into the mass below it, and the system starts to fall.
#Now calculate subsequent ps after each time step.
for instant in range(1, numInstants):

  # Update the positions based on their current (stored) velocity
  updatePositions()

  # Now update the velocities based on their current (calculated) acceleration
  updateVelocities()


# Get the times and points in a form ready to graph
times = np.linspace(0, simulationTime, numInstants)
points = np.transpose(ps)

# Plot the particles' trajectories
for n in range(N):
  plt.plot(times, points[n], label= ("Particle " + str(n)))
  pass

# Get the polynomial fit starting after some fraction of the trajectory of the first particle
fraction = 1.0/3
startIndex = int(fraction * len(points[0]))

# DIFFERENT VERSION USING NP.POLYFIT INSTEAD OF STATS.LINREGRESS
# poly = np.polyfit(times[startIndex:], points[0][startIndex:], 1, None, True)
# coeff = poly[0]
# print ("Coefficients: " + str(poly[0]))
# print ("Residual Sum: " + str(poly[1]))

slope, intercept, r_value, p_value, std_err = stats.linregress(times[startIndex:], points[0][startIndex:])
rValueDigits = 3
rValue = int((-r_value)*(10**rValueDigits))/float(10**rValueDigits) #Truncate

# Plot the linear approximation to the first particle's trajectory,
# and a vertical line indicating where the approximation starts
plt.axvline(x = times[startIndex], label = "Lin Approx of P0 starts")
plt.plot(times, [(slope*t + intercept) for t in times], label="Lin Approx")

# Add labels and show the graph
plt.xlabel('time / s')
plt.ylabel('position / m')
plt.title("Graph of the particles' positions over time" + "\n" + "R-Value: " + str(rValue) + "; Hanging Mass: " + str(M) + "kg")
plt.legend()

plt.show(block = False)
x = raw_input() #Don't exit immediately -- wait for some stdin input
