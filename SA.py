from SALib.sample import saltelli
from SALib.analyze import sobol
import numpy as np
from scipy.integrate import odeint
from scipy import integrate
#from scipy.optimize import minimize
#import matplotlib.pyplot as plt
#
#=================================================
# No. of disacch.
# n = 35 and 40 chosen to implement with.
n = 35
#
#=================================================
# Solutions to model
def sol(p,initial_cond,t0,t_end,stpz):
	# Parameters
	k1, k2, k3, k5, k6 = p
	# time-grid-----
	t = np.arange(t0, t_end, stpz)
	# Model
	def funct(y,t):
		# if n = 4,
		# y[0] = E
		# y[i] = D[i], i=1,n
		# y[i] = EoD[i-n], i=n+1,2n, 5,6,7,8
		# y[i] = EvD[i-2n+1], i=2n+1,3n-1, 9,10,11
		# y[i] = ExD[i-3n+3], i=3n,4n-3, 12,13
		#n = 20
		#k = 4*n - 2
		# the model equations
			# sum of D[i], i=1,n
		#sumD = sum(y[i] for i in range(1,n+1))
			# sum of i*D[i], i=1,n
		#sumiD = sum(i*y[i] for i in range(1,n+1))
			# sum of EoD[i], i=1,n
		#sumC = sum(y[i] for i in range(n+1,2*n+1))
			# sum of EvD[i], i=2,n
		#sumT = sum(y[i] for i in range(2*n+1,3*n))
			# sum of ExD[i], i=3,n
		#sumX = sum(y[i] for i in range(3*n,4*n-2)
		#===============================
		f = []
		# Eq. of Enzyme = y[0]
		f.append(-k1*y[0]*sum(i*y[i] for i in range(1,n+1)) + k5*y[n+1]
			+ k2*(sum(y[i] for i in range(2*n+1,3*n))
			+ sum(y[i] for i in range(3*n,4*n-2))))
		# Eq. of D[1] = y[1]
		f.append(-k1*y[0]*y[1] + k5*y[n+1]
			+ k3*sum(y[i] for i in range(2*n+1,3*n)))
		# Eq. of D[2] = y[2]
		f.append(-k1*2*y[0]*y[2] + k2*y[2*n+1]
			+ k3*sum(y[i]/(i-3*n+1) for i in range(3*n,4*n-2)))
		# Eqs. of D[i] = y[i], i=3,n-1
		for i in range(3,n):
			f.append(-k1*i*y[0]*y[i] + k2*y[2*n+i-1] + k2*y[3*n+i-3]
				+ k3*sum(y[k]/(k-3*n+1) for k in range(3*n+1+i-3,4*n-2)))
		# Eq. of D[n] = y[n]
		f.append(-k1*n*y[0]*y[n] + k2*y[3*n-1] + k2*y[4*n-3])
		# Eq. of EoD[1] = y[n+1]
		f.append(-k5*y[n+1] + k1*y[0]*y[1] + k3*y[2*n+1]
			+ k3*sum(y[i]/(i-3*n+1) for i in range(3*n,4*n-2)))
		# Eqs. of EoD[2] = y[n+2], 2<=i<=n-2
		for i in range(n+2,2*n-1):
			f.append(-k5*y[i] + k1*y[0]*y[i-n] + k3*y[i+n] + k6*y[i+n-1]
				+ k3*sum(y[k]/(k-3*n+1) for k in range(i+2*n-1,4*n-2)))
		# Eq. of EoD[n-1] = y[2*n-1]
		f.append(-k5*y[2*n-1] + k1*y[0]*y[n-1] + k3*y[3*n-1] + k6*y[3*n-2])
		# Eq. of EoD[n] = y[2*n]
		f.append(-k5*y[2*n] + k1*y[0]*y[n] + k6*y[3*n-1])
		# Eqs. of EvD[i] = y[2*n+i-1], i=2,n
		for i in range(2*n+1,3*n):
			f.append(-(k2 + k3 + k6)*y[i] + k5*y[i-n+1] + k1*y[0]*y[i-2*n+1])
		# Eqs. of ExD[i] = y[3n+i-3], i=3,n
		for i in range(3*n,4*n-2):
			f.append(-(k2 + k3)*y[i] +k1*(i-3*n+1)*y[0]*y[i-3*n+3])
		return(f)
	#===============================
	# integrate the system---
	ds = integrate.odeint(funct,initial_cond,t)
	return(ds)
#
#===============================
# Model for SA of parameters in initial conditions
# It's impossible because of lack of data
# for each sample of initial conditions
#if 1 == 10:
	#===============================
	# Data---
	#Td = [0.0, 0.3528, 0.9119, 1.8932, 3.9141, 6.0, 24.0, 48.0]
	#Zd = [180.1559*0.025/n, 0.6969, 1.3391, 2.7056, 4.5154, 4.6042, 4.6947, 4.7562]
	#indices = [int(x*1e+4) for x in Td]
	#
	#===============================
	# Model fo SA
	#def model(ICs):
		#E, D = ICs
		# Guess paraneters
		#p0 = [10000.0, 200.0, 2000.0, 3000.0, 100.0]
		#def score(p):
			# Initial conditions
			#y0 = []
			#y0.append(E)
			#for i in range(1,n):
				#y0.append(0.0)
			#y0.append(D)
			#for i in range(n+1,4*n-2):
				#y0.append(0.0)
			#y0
			#
			# Solutions
			#y = sol(p, y0, 0.0, 48.0 + 1e-4, 1e-4)
			# Reducing ends
			#r_ends = 180.1559*sum(y[:,i] for i in range(1,4*n-2))
			#Zm = np.take(r_ends, indices)
			#ss = sum((x - y)**2 for x, y in zip(Zm, Zd))
			#return(ss)
		#bnds = ((5000,3e+4), (0,1e+3), (500,1e+4), (500,1e+4), (0,1e+2))
		#s_fit = minimize(score,p0,method='L-BFGS-B',bounds=bnds,options={'eps': 1e-4, 'disp': False, 'ftol':1e-06, 'maxiter':1000000},tol=1e-6)
		#s_fit =	minimize(score,p0,method='SLSQP',bounds=bnds,options={'eps': 1e-04, 'disp': False, 'ftol':1e-06, 'maxiter':1000000})
		#new_p = s_fit.x
		#new_fit = score(new_p)
		#return(new_p)
	#
	#===============================
	# Define problem of SA
	#problem = {
		#'num_vars': 2,
		#'names': ['E', 'D'],
		#'bounds':[[0.9*1.21e-5, 1.1*1.21e-5], [0.9*2.5e-2/n, 1.1*2.5e-2/n]]
	#}
	#
	#===============================
	# Generate samples
	#param_values = saltelli.sample(problem, 2000)
	#
	#===============================
	# Number of model points
	#for i in range(5):
		#Y = np.zeros([param_values.shape[0]])
		#for j, X in enumerate(param_values):
			#K = model(X)[i]
			#Y[j] = K
		# Perform analysis
		#Si = sobol.analyze(problem, Y, print_to_console=True)
		#print('for i = ' + str(i) + '\n')
		#print(Si['S1'])
		#print(Si['ST'])
		#
#====================================================
# Model for Sensitivity Analysis
def model(p):
	#===============================
	# initial conditions
	y0 = []
	E0 = 1.21e-5 # 0.00001204 mmol/ml
	D0 = 2.64e-2/n # 0.02491957867/n
	y0.append(E0)
	for i in range(1,n):
		y0.append(0.0)
	y0.append(D0)
	for i in range(n+1,4*n-2):
		y0.append(0.0)
	y0
	#==============================
	# Solutions
	y = sol(p,y0,0.0,6.0 + 0.02, 0.01)
	# Disaccharide in mg/ml
	F = 401.3*y[:,1]
	return(F)
	#
#=================================================
# Define problem of sensitivity analysis
problem = {
	'num_vars': 5,
	'names': ['k1', 'k2', 'k3', 'k5', 'k6'],
	'bounds': [[0.9*9998.23,1.1*9998.23 ], [0.9*109.36,1.1*109.36], [0.9*2695.43,1.1*2695.43], [0.9*2096.69,1.1*2096.69], [0.9*4.7,1.1*4.7]]
}
#
# Generate samples
param_values = saltelli.sample(problem, 1000)
#
#=================================================
# Number of model points at which implementing SA
print('For n = ' + str(n) + '\n')
L = [101, 201, 301, 451, 601]
for i in L:
	# Run model
	Y = np.zeros([param_values.shape[0]])
	for j, X in enumerate(param_values):
		K = model(X)
		Y[j] = K[i]
	# Perform analysis
	Si = sobol.analyze(problem, Y, print_to_console=True)
	#
	# Print the first-order sensitivity indices
	print(' and for i = ' + str(i) + '\n')
	print(Si['S1'])
	# Print the total-order sensitivity indices
	print(Si['ST'])
	#print(param_values)
	#print(Y)
#=========================================================
