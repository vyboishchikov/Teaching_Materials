import sys, numpy as np, matplotlib.pyplot as plt, scipy.integrate

def GETSECTION(Header,a):
	Result = []
	for i,p in enumerate(a):
		if p in Header: break
	for p in a[i+1:]:
		if p != '':
			Result.append(p)
		else:
			break
	return Result

def CALCULATE_MATRIX(a):
	NReactions = len(a)
	Constants = np.array([float(r.split()[0]) for r in a])
	Noms = []

	for r in a:
		r0 = ''
		for p in r.split()[1:]:
			r0 += p
		for q in r0.replace('>','+').split('+'):
			if not q.split('*')[-1] in Noms:
				Noms.append(q.split('*')[-1])

	NSubst = len(Noms)
	ReactantsDict = dict(zip(Noms, range(NSubst)))
	
	Matriu = np.zeros((NSubst,NReactions))
	LlistaLlistes = []
	for i,r in enumerate(a):
		LlistaLlistes.append([])
		r0 = ''
		for p in r.split()[1:]:
			r0 += p
		# Left-hand part
		for q in r0.split('>')[0].split('+'):
			qsplit = q.split('*')
			if len(qsplit)==1:
				CoefSubst = 1
			else:
				CoefSubst = int(qsplit[0])
			IndSubst = ReactantsDict[qsplit[-1]]
			for j in range(CoefSubst):
				LlistaLlistes[-1].append(IndSubst)
			Matriu[IndSubst,i] += -CoefSubst
		# Right-hand part
		for q in r0.split('>')[1].split('+'):
			qsplit = q.split('*')
			if len(qsplit)==1:
				CoefSubst = 1
			else:
				CoefSubst = int(qsplit[0])
			IndSubst = ReactantsDict[qsplit[-1]]
			Matriu[IndSubst,i] += CoefSubst
	return Matriu, LlistaLlistes, Noms, Constants, ReactantsDict

def OBTAIN_INITIAL_CONCENTRATIONS(a,ReactantsDict):
	InitialConcentrations = np.zeros(len(ReactantsDict))
	ListFrozen = []
	for r in a:
		InitialConcentrations[ReactantsDict[r.split()[0]]] = float(r.split()[1])
		if len(r.split())>2 and r.split()[2].lower() in ['frozen','freeze']:
			ListFrozen.append(ReactantsDict[r.split()[0]])
	return InitialConcentrations, ListFrozen

def OBTAIN_SIMULATION_TIME(a):
	if len(a)==0:
		return 400.0
	else:
		return float(a[0].split()[0])

def OBTAIN_TIME_STEP(a):
	if len(a)==0:
		return np.inf
	else:
		return float(a[0].split()[0])

def OBTAIN_INTEGRATION_METHOD(a):
	AllMethodsList = ['RK45','RK23','DOP853','Radau','BDF','LSODA']
	if len(a)==0:
		return 'Radau' # default
	else:
		Method = a[0].split()[0]
		if not Method in AllMethodsList:
			print('Wrong method ' + Method + '\nMethods available: ' +', '.join(AllMethodsList),file=sys.stderr)
			quit() 
		return Method

def OBTAIN_SUBSTANCES_TO_PRINT(a,Noms):
	if len(a)==0:
		return range(len(Noms))
	else:
		Result = []
		for r in a:
			for q in r.split():
				Result.append(ReactantsDict[q])
		return Result 

def OBTAIN_SUBSTANCES_TO_DRAW(a):
	Result = []
	for r in a:
		for q in r.split():
			Result.append(ReactantsDict[q])
	return Result 

def OBTAIN_UNIT(a):
	if len(a)==0:
		return ''
	else:
		return a[0].strip()

#def Deviations(ActualConstants):
#	Constants = ActualConstants.copy()
#	Solution = solve_ivp(Kinetics, (0.0, SimulationTime), InitialConcentrations, method=IntegrationMethod, max_step=TimeStep)
#	Solution.y
#	ExpConcMatrix
	
def Kinetics(t, Concentration):
	ReactionBlock = Constants.copy()
	for i,r in enumerate(LlistaLlistes):
		for j in r:
			ReactionBlock[i] *= Concentration[j]
	Deriv = np.dot(Matriu,ReactionBlock)
	Deriv[ListFrozen] = 0.0
	return Deriv

def PRINT_RESULTS(t,y,Noms,Substances):
	print('   '+'\t'.join(['  Time'] + ['%18s'%Noms[r]+'=const' if r in ListFrozen else '%18s'%Noms[r] for r in Substances]))
	for i,Time in enumerate(t):
		print('\t'.join(['%12.5e'%Time] + ['%18.11e'%y[j][i] for j in Substances]))

def DRAW_RESULTS(t, y, Noms, ToDraw, TimeUnit, ConcentrationUnit):
	if len(ToDraw) > 0:
		for i in ToDraw:
			plt.plot(t,y[i],label=Noms[i])
		#plt.title('Plot ',fontsize=15)
		Xlabel = 'Time'
		Ylabel = 'Concentration'
		if TimeUnit != '': Xlabel += ' / '+TimeUnit
		if ConcentrationUnit != '': Ylabel += ' / '+ConcentrationUnit
		plt.xlabel(Xlabel)    #,fontsize=13
		plt.ylabel(Ylabel)   #,fontsize=13
		plt.legend()
		#plt.yscale('log')
		plt.show()

if len(sys.argv) < 2:
	CoreName = sys.argv[0].split('\\')[-1]
	print('Usage:\n'+CoreName+' file1\n    or\n'+CoreName+' file1 file 2',file=sys.stderr)
	quit()
with open(sys.argv[1],'r') as Fitxer:
	a = [r.split('#')[0].strip() for r in Fitxer.readlines() if r[0:1]!='#']

ReactionsText = GETSECTION(['REACTIONS:'],a)
InitialConcentrationsText = GETSECTION(['INITIAL CONCENTRATIONS:','INITIAL_CONCENTRATIONS:','INITIALCONCENTRATIONS:'],a)
SimultationTimeText = GETSECTION(['SIMULATION TIME:','SIMULATIONTIME:','SIMULATION_TIME:'],a)
TimeStepText = GETSECTION(['TIMESTEP:', 'TIME STEP:', 'TIME_STEP:'],a)
IntegrationMethodText = GETSECTION(['METHOD:'],a)
ToPrintText = GETSECTION(['TOPRINT:','TO PRINT:','TO_PRINT:'],a)
ToDrawText = GETSECTION(['TODRAW:','TO DRAW:','TO_DRAW:'],a)
TimeUnitText = GETSECTION(['TIMEUNIT:','TIME UNIT:','TIME_UNIT:'],a)
ConcentrationUnitText = GETSECTION(['CONCENTRATIONUNIT:','CONCENTRATION UNIT:','CONCENTRATION_UNIT:'],a)

Matriu, LlistaLlistes, Noms, Constants, ReactantsDict = CALCULATE_MATRIX(ReactionsText)
InitialConcentrations, ListFrozen = OBTAIN_INITIAL_CONCENTRATIONS(InitialConcentrationsText,ReactantsDict)
print('InitialConcentrations:',InitialConcentrations)
SimulationTime = OBTAIN_SIMULATION_TIME(SimultationTimeText)
print('SimulationTime:',SimulationTime)
TimeStep = OBTAIN_TIME_STEP(TimeStepText)
if TimeStep == np.inf:
	print('Timestep: automatic')
else:
	print('Timestep:',TimeStep )
IntegrationMethod = OBTAIN_INTEGRATION_METHOD(IntegrationMethodText)
print('IntegrationMethod:',IntegrationMethod)
SubstancesToPrint = OBTAIN_SUBSTANCES_TO_PRINT(ToPrintText,Noms)
SubstancesToDraw = OBTAIN_SUBSTANCES_TO_PRINT(ToDrawText,Noms)
ConcentrationUnit = OBTAIN_UNIT(ConcentrationUnitText)
TimeUnit = OBTAIN_UNIT(TimeUnitText)

Solution = scipy.integrate.solve_ivp(Kinetics, (0.0, SimulationTime), InitialConcentrations, method=IntegrationMethod, max_step=TimeStep)
PRINT_RESULTS(Solution.t, Solution.y, Noms, SubstancesToPrint)
DRAW_RESULTS(Solution.t, Solution.y, Noms, SubstancesToDraw, TimeUnit, ConcentrationUnit)
