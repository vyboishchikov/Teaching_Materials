# Non-linear regression for chemical kinetics

import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt

# Define the model function for y = A * exp(-k * t**2 - k1*t)
def FirstOrderKinetics(t, Y_inf, Y0, k):
	return Y_inf + (Y0-Y_inf)*np.exp(-k*t)
	
def SecondOrderKinetics(t, Y_inf, Y0, k):
	return Y_inf + (Y0-Y_inf)/(1+A0*k*t)

def ZerothOrderKinetics(t, Y_inf, Y0, k):
	return Y0 - (Y0-Y_inf)*(k/A0)*t

A0 = 0.01 # initial concentration (non needed for first.order kinetics)
Function = SecondOrderKinetics # choose the order you need

tData = [ # time data points
  29
 ,30
 ,31
 ,32
 ,33
 ,34
 ,35
]

yData = [ # observable (e.g. conductivity) data points
  1223
 ,1212
 ,1196
 ,1182
 ,1175
 ,1164
 ,1157
]


k_init = 1e-12 # initial k value

# Perform the non-linear regression using curve_fit
# Initial guess for parameters Y_inf, Y0, k
InitialGuess = [yData[-1], yData[0], k_init] # initial Y_inf, Y0, k values
bounds = ([-np.inf for _ in InitialGuess],[np.inf for _ in InitialGuess])
bounds[0][-1] = 0

# Fit the model to the data
OptimizedParameters = scipy.optimize.curve_fit(Function, tData, yData, p0=InitialGuess,method='dogbox',bounds=bounds)[0]

# Generate fitted curve
tFit = np.linspace(min(tData), max(tData), 1000)
yFit = Function(tFit, *OptimizedParameters)

# Print optimized parameters
print("Optimized parameters (Y_inf, Y0, k):", *OptimizedParameters)

# Plot data and fitted curve
plt.scatter(tData, yData, label='Data', color='red')
plt.plot(tFit, yFit, label='Fitted Curve', color='blue')
plt.xlabel('t')
plt.ylabel('y')
plt.legend()
plt.show()
