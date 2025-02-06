# Non-linear regression for chemical kinetics
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt

# Integrated kinetic rate law
def FirstOrderKinetics(t, Y_inf, Y0, k):
	return Y_inf + (Y0-Y_inf)*np.exp(-k*t)
	
def SecondOrderKinetics(t, Y_inf, Y0, k):
	return Y_inf + (Y0-Y_inf)/(1+A0*k*t)

def ZerothOrderKinetics(t, Y_inf, Y0, k):
	return Y0 - (Y0-Y_inf)*(k/A0)*t

Function = SecondOrderKinetics

A0 = 0.01 # define your initial concentration here; not needed for 1st-order kinetics

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


if len(tData)>len(yData):
	print(f"There are {len(tData)} time points, but only {len(yData)} data points.\nCheck your data")
	quit()
elif len(tData)<len(yData):
	print(f"There are {len(yData)} data points, but only {len(tData)} time points.\nCheck your data")
	quit()

# Initial guess for parameters Y_inf, Y0, k
k_init = 1e-2 # initial k value
InitialGuess = [yData[-1], yData[0], k_init] # initial Y_inf, Y0, k values
bounds = ([-np.inf for _ in InitialGuess],[np.inf for _ in InitialGuess])
bounds[0][-1] = 0

# Non-linear regression
try:
	OptimizedParameters, pcov = scipy.optimize.curve_fit(Function, tData, yData, p0=InitialGuess,bounds=bounds,method='dogbox')
except:
	print("Fitting failed. Try to change the initial k value")
	quit()

# Generate fitted curve
tFit = np.linspace(min(tData), max(tData), 1000)
yFit = Function(tFit, *OptimizedParameters)

print('Optimized parameters (Y_inf, Y0, k): '+ ', '.join(['%.3g +- %.1g'%p for p in zip(OptimizedParameters,np.sqrt(np.diag(pcov)))]))

# Plot
plt.scatter(tData, yData, label='Data', color='red',s=6)
plt.plot(tFit, yFit, label='Fitted Curve', color='blue')
plt.xlabel('t')
plt.ylabel('y')
plt.legend()
plt.show()
