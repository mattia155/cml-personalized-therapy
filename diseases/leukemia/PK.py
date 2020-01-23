#!/usr/bin/env python
# coding: utf-8

################################################################################
# Direct Multiple Shooting module: Direct multiple shooting methods for leukemia
# Description: Implementation of direct multiple shooting for continuous and 
#              discrete leukemia problem (Copy and re-arrange for new disease)                                                             
################################################################################



################################### Imports ####################################
from casadi import *
import matplotlib.pyplot as plt
import math
import functionalcost as fc
################################################################################



################################## Functions ###################################

#-------------------------------------------------------------------------------
# Heaviside step function
def theta(tref,time):
	return (1 * (time-tref >= 0))
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Define discrete or continuos version of the leukemia problem
def personalPK(time, intervalBetweenDoses, numberOfControlIntervals, 
               rungeKuttaSteps, functionalCost, clearance, volume,
               constant_ka, phi=0, ec50=0, constant_dosage=0):
	
	# Declare model constants
	T = time # Time horizon
	K = numberOfControlIntervals # number of control intervals
	ff=1
	ka = constant_ka
	CL = clearance
	V = volume
	if(functionalCost != 1):
		ctg=1
	else:
		ctg=0.57 # Literature
	
	# Declare model variables
	t = MX.sym('t')
	xb = MX.sym('xb')
	d = MX.sym('d1')
	if(intervalBetweenDoses != 0):
		# Discrete dosage scenario
		n_doses = math.floor(T/intervalBetweenDoses) #hours between doses 
		for i in range(1,n_doses):
			di = MX.sym('d'+str(i+1))
			d = vertcat(d, di)
	else:
		# Continuous dosage scenario
		n_doses = 0
	
	# Model equations (ODE)
	if(n_doses != 0):
		# Discrete dosage scenario
		xdot_summatory = 0
		if(constant_dosage != 0):
			for i in range(n_doses):
				xdot_summatory = xdot_summatory + (theta(i*intervalBetweenDoses,t)*ff*constant_dosage*exp((-ka)*(t-i*intervalBetweenDoses)))
		else:
			for i in range(n_doses):
				xdot_summatory = xdot_summatory + (theta(i*intervalBetweenDoses,t)*ff*d[i]*exp((-ka)*(t-i*intervalBetweenDoses)))
		xdot = ka * xdot_summatory - (CL/V) * xb
	else:
		# Continuous dosage scenario
		if(constant_dosage != 0):
			xdot = constant_dosae - (CL/V)*xb
		else:
			xdot = d - (CL/V)*xb
			
	# Objective term (Lagrange and Mayer terms)
	if(functionalCost >= 2):
		w1 = phi
		w2 = 1
		a = 0.87
		p = 0.45
		emax = 1.0
		k = 377
	L, Mayer = fc.select_cost(functionalCost)
	# Lagrange term
	L = eval(L)
	# Mayer term
	Mayer = eval(Mayer)

	# Formulate discrete time dynamics
	M = rungeKuttaSteps # RK4 steps per interval (x(t)-->x(t+1) in M steps)
	DT = T/K/M
	f = Function('f', [xb, d, t], [xdot, L]) #f: x,u,t --> xdot,L
	X0 = MX.sym('X0')
	if(n_doses != 0):
		# Discrete dosage scenario 
		U = MX.sym('U',n_doses) #u=[u1,u2,...,n_doses]
	else:
		# Continuous dosage scenario
		U = MX.sym('U')
	X = X0 #X new state (starts as X0)
	Q = 0 #new control
	for j in range(M): 
		k1, k1_q = f(X, U, t)
		k2, k2_q = f(X + DT/2 * k1, U, t)
		k3, k3_q = f(X + DT/2 * k2, U, t)
		k4, k4_q = f(X + DT * k3, U, t)
		X=X+DT/6*(k1 +2*k2 +2*k3 +k4)
		Q = Q + DT/6*(k1_q + 2*k2_q + 2*k3_q + k4_q)
	# Inputs [x_0,U,t] and outputs [X,Q] computed through Runge-Kutta
	F = Function('F', [X0, U, t], [X, Q],['x0','p', 't'],['xf','qf'])
	MY = Function('MY', [xb], [Mayer], ['xf'], ['mf'])

	# Define lower and upper bound for dosage
	lowerbound = 0
	upperbound = 1000 #more than 1000mg could be lethal dosage
	
	# Return variables
	return n_doses, F, MY, lowerbound, upperbound
	
#-------------------------------------------------------------------------------

################################################################################
