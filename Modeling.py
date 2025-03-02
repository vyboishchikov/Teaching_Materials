# model reaction: ozone decomposition
import scipy.integrate, numpy

k1 = 4.7e-7
k_1 = 4.2e-11
k2 = 9.3e-7

def Kinetics(t, Concentrations):
	O3, O2, O = Concentrations
	v1 = k1*O3      # use **2 if you need a concentration squared
	v_1 = k_1*O2*O
	v2 = k2*O*O3
	dO3 = -v1 + v_1 - v2
	dO2 = v1 - v_1 + 2*v2
	dO = v1 - v_1 - v2
	return dO3, dO2, dO

SimulationTime = 63115200.0 # two years in seconds; replace this by a value suitable for your problem
InitialConcentrations = 7.4e12, 1.5e17, 0.0
Substances = "O3, O2, O"

Separator = "," # adjust the field separator to your system; may be "," or ";", or "\t"

Solution = scipy.integrate.solve_ivp(Kinetics,(0, SimulationTime), InitialConcentrations, 
	method='LSODA', t_eval=numpy.linspace(0,SimulationTime,1000),max_step=SimulationTime/10000)

if not Solution.success: print("Error message:",Solution.message)
print("Time"+Separator+Separator.join(Substances.replace(" ","").split(",")))
for t,Y in zip(Solution.t,Solution.y.T):
	print("%12.7g"%t+Separator+Separator.join(["%12.7g"%y for y in Y]))
