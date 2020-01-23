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
# Define test problem
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
	ctg = 0.75
	n_doses = 0
	
	# Declare model variables
	xb = MX.sym('xb')
	d = MX.sym('d1')
	t = MX.sym('t') # Time
	
	# Model equations (ODE)
	xdot = d-(CL/V)*xb + ka*ec50
	
	# Objective term (Lagrange and Mayer terms)
	if(functionalCost >= 2):
		w1 = phi
		w2 = 1
		a = 0.8
		p = 0.4
		emax = 0.9
		k = 301
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
	upperbound = 300 
	
	# Return variables
	return n_doses, F, MY, lowerbound, upperbound
	
#-------------------------------------------------------------------------------

################################################################################
