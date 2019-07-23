# Appendix A. Supplementary material
#
#
# The program used to estimate parameters and
# calculate sensitivity indices. It can be executed
# with the following libraries in Python 2.7!
#
# Libraries for parameter estimation
import numpy as np
from scipy.integrate import odeint
from scipy import integrate
from scipy.optimize import minimize
#
# Library for plotting the results
import matplotlib.pyplot as plt
#
# Libraries for doing sensitivity anealysis
from SALib.sample import saltelli
from SALib.analyze import sobol
#
#=================================================
# n: number of disaccharides of a HA chain.
# Take care when choosing n for sensitivity analysis
n = 125
#
#=================================================
# Create solution to the model that is a vector-value function
# of p, initial_cond, t0, t_end and stpz, where p, initial_cond,
# t0, t_end and stpz are a parameter vector, a vector of initial
# concentrations, starting time, ending time, and step size of time
# interval, respectively.
#
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
		# Eq. for Enzyme = y[0]
		f.append(-k1*y[0]*sum(i*y[i] for i in range(1,n+1)) + k5*y[n+1]
			+ k2*(sum(y[i] for i in range(2*n+1,3*n))
			+ sum(y[i] for i in range(3*n,4*n-2))))
		# Eq. for D[1] = y[1]
		f.append(-k1*y[0]*y[1] + k5*y[n+1]
			+ k3*sum(y[i] for i in range(2*n+1,3*n)))
		# Eq. for D[2] = y[2]
		f.append(-k1*2*y[0]*y[2] + k2*y[2*n+1]
			+ k3*sum(y[i]/(i-3*n+1) for i in range(3*n,4*n-2)))
		# Eqs. for D[i] = y[i], i=3,n-1
		for i in range(3,n):
			f.append(-k1*i*y[0]*y[i] + k2*y[2*n+i-1] + k2*y[3*n+i-3]
				+ k3*sum(y[k]/(k-3*n+1) for k in range(3*n+1+i-3,4*n-2)))
		# Eq. for D[n] = y[n]
		f.append(-k1*n*y[0]*y[n] + k2*y[3*n-1] + k2*y[4*n-3])
		# Eq. of EoD[1] = y[n+1]
		f.append(-k5*y[n+1] + k1*y[0]*y[1] + k3*y[2*n+1]
			+ k3*sum(y[i]/(i-3*n+1) for i in range(3*n,4*n-2)))
		# Eqs. for EoD[2] = y[n+2], 2<=i<=n-2
		for i in range(n+2,2*n-1):
			f.append(-k5*y[i] + k1*y[0]*y[i-n] + k3*y[i+n] + k6*y[i+n-1]
				+ k3*sum(y[k]/(k-3*n+1) for k in range(i+2*n-1,4*n-2)))
		# Eq. for EoD[n-1] = y[2*n-1]
		f.append(-k5*y[2*n-1] + k1*y[0]*y[n-1] + k3*y[3*n-1] + k6*y[3*n-2])
		# Eq. for EoD[n] = y[2*n]
		f.append(-k5*y[2*n] + k1*y[0]*y[n] + k6*y[3*n-1])
		# Eqs. for EvD[i] = y[2*n+i-1], i=2,n
		for i in range(2*n+1,3*n):
			f.append(-(k2 + k3 + k6)*y[i] + k5*y[i-n+1] + k1*y[0]*y[i-2*n+1])
		# Eqs. for ExD[i] = y[3n+i-3], i=3,n
		for i in range(3*n,4*n-2):
			f.append(-(k2 + k3)*y[i] +k1*(i-3*n+1)*y[0]*y[i-3*n+3])
		return(f)
	#===============================
	# integrate the system---
	ds = integrate.odeint(funct,initial_cond,t)
	return(ds)
#
#===============================
# Section for Parameter estimation
# Molecular weight of a disaccharide unit
#
b = 401.30
#
# Initial conditions
#
y0 = []
E0 = 1.21e-5
D0 = 2.64e-4/26**2		# 2.64e-2/n (real)
y0.append(E0)
for i in range(1,75):
	y0.append(0.0)
for i in range(75,101):
	y0.append((i-74)*D0)
for i in range(101,n+1):
	y0.append((n+1-i)*D0)
for i in range(n+1,4*n-2):
	y0.append(0.0)
y0
#
#===============================
# Time grid 1-----
stpz = 1e-4
t0 = 0.0
t_end = 48.0 + stpz
t1 = np.arange(t0, t_end, stpz)
#
#===============================
# Original guesses
#
p0 = [9998.23, 109.36, 2695.43, 2096.69, 4.70]
#
#===============================
# Data section
	# Data
	# Time grid
Td = [0.0, 0.3528, 0.9119, 1.8932, 3.9141, 6.0, 24.0, 48.0]
        # Pneumococcal
Zd = [180.1559*D0, 0.6969, 1.3391, 2.7056, 4.5154, 4.6042, 4.6947, 4.7562]
indices = [int(x*1e+4) for x in Td]
#
#===============================
# Section for parameter estimation
	# Score fit of the model to data
def score(p):
	# Solutions
	y = sol(p, y0, 0.0, t_end, stpz)
	# Reducing ends
	r_ends = 180.1559*sum(y[:,i] for i in range(1,4*n-2))
	Zm = np.take(r_ends, indices)
	ss = sum((x - y)**2 for x, y in zip(Zm, Zd))
	return(ss)
#
#===============================
# Minimize the score
#
# Original bounds
#
#bnds = ((5000,3e+4), (0,1e+3), (500,5e+3), (500,5e+3), (0,1e+2))
#
#s_fit = minimize(score,p0,method='L-BFGS-B',bounds=bnds,options={'eps': 1e-4, 'disp': False, 'ftol':1e-06,'maxiter':1000000},tol=1e-6)
#s_fit =	minimize(score,p0,method='SLSQP',bounds=bnds,options={'eps': 1e-04, 'disp': True, 'ftol':1e-06, 'maxiter':1000000})
s_fit = minimize(score, p0, method='nelder-mead', options={'xtol':1e-4, 'disp': True})
new_p = s_fit.x
print('for n =' + str(n))
print(new_p)
#
#---------------------------------------------------------------------------------
# Section for plotting
# the model curve and the data
#
if 1 == 10:
    	y_sol = sol(p0, y0, t0, t_end, stpz)
 	#
    	# Plot Reducing sugars
	fg1 = plt.figure(1)
	plt.plot(t1, 180.1559*sum(y_sol[:,i] for i in range(1,4*n-2)))
	plt.plot(Td, Zd, 'ro')
	plt.legend(['Model curve', 'Data'])
	plt.xlabel('TIME (HOURS)')
	plt.ylabel('INCREASE IN REDUCING ENDS' + '\n' + 'AS GLUCOSE $mg/ml$')
	#fg1.savefig('data.eps')
	plt.show()

