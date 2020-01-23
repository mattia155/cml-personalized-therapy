#!/usr/bin/env python
# coding: utf-8

################################################################################
# Direct Multiple Shooting module: Direct multiple shooting methods for leukemia
# Description: Implementation of direct multiple shooting for continuous and 
#              discrete test problem                                                             
################################################################################



################################### Imports ####################################
from casadi import *
import matplotlib.pyplot as plt
import math
import functionalcost as fc
################################################################################



################################## Functions ###################################
# - heaviside step function
# - direct single shooting on discrete problems
# - direct single shooting on continuos problems
# - direct multiple shooting on discrete problems
# - direct multiple shooting on continuos problems
# - print results

#-------------------------------------------------------------------------------
# Heaviside step function
def theta(tref,time):
	return (1 * (time-tref >= 0))
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Choose discrete or continuos version of the problem
def solve(time, intervalBetweenDoses, numberOfControlIntervals, 
		  rungeKuttaSteps, functionalCost, clearance, volume, 
		  constant_ka, phi=0, ec50=0, constant_dosage=0):
	if(intervalBetweenDoses == 0):
		x, f1, f2 = solveContinuos(time, numberOfControlIntervals, 
						           rungeKuttaSteps, functionalCost, clearance, 
						           volume, constant_ka, phi=0, 
						           ec50=0, constant_dosage=0)
	#else:
		#does not exist a discrete case for test disease
		
	return x, f1, f2
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------	
# Direct multiple shooting to resolve the continuos case of the problem
def solveContinuos(time, numberOfControlIntervals, rungeKuttaSteps, functionalCost, clearance, volume, constant_ka, phi=0, ec50=0, constant_dosage=0):
	
	# Declare model constants
	T = time # Time horizon
	K = numberOfControlIntervals # number of control intervals
	ff=1
	ka = constant_ka
	CL = clearance
	V = volume
	functionalCost = None
	ctg = 0.75
	# Declare model variables
	xb = MX.sym('xb')
	d = MX.sym('d1')
	
	# Model equations (ODE)
	xdot = d-(CL/V)*xb + ka*ec50
	# Cost
	L, Mayer = fc.select_cost(functionalCost)
	L = eval(L)
	Mayer = eval(Mayer)
	
	# Formulate discrete time dynamics
	M = rungeKuttaSteps # RK4 steps per interval (x(t)-->x(t+1) in 4 steps)
	DT = T/K/M
	f = Function('f', [xb, d, t], [xdot, L]) #f: x,u,t --> xdot,L
	X0 = MX.sym('X0') 
	U = MX.sym('U') #u=[u1] continuos dosage
	X = X0 #X new state (starts as X0)
	Q = 0 #new control
	for j in range(M): 
		k1, k1_q = f(X, U)
		k2, k2_q = f(X + DT/2 * k1, U)
		k3, k3_q = f(X + DT/2 * k2, U)
		k4, k4_q = f(X + DT * k3, U)
		X=X+DT/6*(k1 +2*k2 +2*k3 +k4)
		Q = Q + DT/6*(k1_q + 2*k2_q + 2*k3_q + k4_q)
	# Inputs [x_0,U,t] and outputs [X,Q] computed through Runge-Kutta
	F = Function('F', [X0, U], [X, Q],['x0','p'],['xf','qf'])
	MY = Function('MY', [xb], [Mayer], ['xf'], ['mf'])

	# Initialize an empty NLP
	w=[] #vector of controls and state over time (u1,x1,u2,x2,...)
	w0 = [] #initial vector of w
	lbw = [] #lower bound for w
	ubw = [] #upper bound for w
	J = 0 #function on (x,u) to be minimized
	g=[] #function on (x,u) (evolution over time)
	lbg = [] #lower bounds on g
	ubg = [] #upper bounds on g

	# "Lift" initial conditions for multiple shooting
	# State Xk
	Xk = MX.sym('X0') #k=0
	w += [Xk]
	lbw += [0]
	ubw += [0]
	w0 += [0]

	# Formulate the NLP problem
	t = MX(0) #maybe useless
	for k in range(K):
		# New NLP variable for the control
		Uk = MX.sym('U_' + str(k))
		w += [Uk]
		w0 += [0]
		lbw += [0]
		ubw += [1000] #constraint control to be under lethal dosage
		# Integrate till the end of the interval
		Fk = F(x0=Xk, p=Uk)
		Xk_end = Fk['xf']
		J=J+Fk['qf']
		# Add mayer term
		Mk = MY(xf=Xk_end)
		J=J+Mk['mf'] #mayer
		# New NLP variable for state at end of interval
		Xk = MX.sym('X_'+str(k+1))
		w += [Xk]
		w0 += [0]
		lbw += [0]
		ubw += [inf]
		# Add equality constraint
		g   += [Xk_end-Xk] #only of the interval k
		lbg += [0]
		ubg += [0]

	# Define a NLP solver and solve the problem.
	prob = {'f': J, 'x': vertcat(*w), 'g': vertcat(*g)}
	solver = nlpsol('solver', 'ipopt', prob, {'ipopt':{'print_frequency_time':60.0}});

	# Solve the NLP
	sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)
	w_opt = sol['x'].full().flatten()
	# State Xb values
	x_opt = w_opt[0::2]/V #start from 0 and select a value every numberOfDoses+1
	# Controls d values
	d_opt = [w_opt[1::2]]
	fig1, fig2 = print_results(x_opt, d_opt, T, K, 1, functionalCost, ctg) #controlla questi parametri e togli quelli inutili
	return x_opt, fig1, fig2
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Print results
def print_results(x_opt, d_opt, T, K, numberOfDoses, functionalCost, ctg): #controlla sti parametri e togli gli inutili
	tgrid = [T/K*k for k in range(K+1)]
	plt.clf()
	plt.ioff()
	plt.figure(1)
	# "discretize" dosage
	d_opt2 = d_opt
	if(numberOfDoses > 0): 
		for i in range(len(d_opt2)):
			start = 0
			summ = 0
			count = 0
			for j in range(len(d_opt2[i])):
				if(d_opt2[i][j] > 1 and start == 0):
					start = 1
					summ = summ + d_opt2[i][j]
					count = count + 1 
				elif(d_opt2[i][j] > 1 and start != 0):
					summ = summ + d_opt2[i][j]
					count = count + 1
				elif(start != 0):
					#end of mean interval
					start = 0
					d_opt2[i][j-count:j] = 0
					d_opt2[i][j] = summ/count #mean
					count = 0
					summ = 0
					
	for i in range(numberOfDoses):
		plt.step(tgrid, vertcat(DM.nan(1), d_opt2[i]), color='blue')
	plt.title('Dosage')
	plt.xlabel('time')
	plt.ylabel('value')
	#legends = []optimal target
	#for i in range(numberOfDoses):
	#	legends += ['d'+str(i+1)]
	#plt.legend(legends)
	plt.grid()
	fig1 = plt.gcf()
	fig1.set_size_inches(5.1,3.75)
	# State Xb
	plt.figure(2)
	plt.plot(tgrid, x_opt, '--', color='red')
	if(functionalCost != 1):
		plt.axhline(ctg, color="gray", label='target') #target
	else:
		plt.axhline(ctg/10, color="gray", label='target') #FATTORE SCALA DA RIMUOVERE
	plt.title('State Xb')
	plt.xlabel('time')
	plt.ylabel('value')
	plt.legend(['xb','target'])
	plt.grid()
	fig2 = plt.gcf()
	
	return fig1, fig2
#-------------------------------------------------------------------------------

################################################################################
