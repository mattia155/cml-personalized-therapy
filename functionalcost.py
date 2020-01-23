#!/usr/bin/env python
# coding: utf-8

################################################################################
# Functional cost module: handle functional costs
# Description: define the possible functional costs and return the selected cost
# (define here new costs)
################################################################################



################################## Functions ###################################
# Cost selector function: select which cost has to be returned
def select_cost(which):
# which could be only 0, 1 or 2 (interface settings, default 0)
	if(which == 0):
		# Functional cost 1: concentration swing around target
		Lagrange = "fabs((xb/V)-ctg)"
		Mayer = "0"
	elif(which == 1):
		# Functional cost 2: keep concentration always above target
		Lagrange = "1000*((xb/V)-ctg<0) + (fabs((xb/V)-ctg))*((xb/V)-ctg>=0)"
		Mayer = "0"
	elif(which == 2):
		# Functional cost 3: tb data
		Lagrange = "w1*((math.log10(math.exp(1))*((2*a-1)*p - k*(emax*(xb/V))/(ec50 + xb/V)))) + w2*(xb/V)"
		Mayer = "0"
	elif(which == 3):
		# Adaptive cost
		Lagrange = "fabs((math.log10(math.exp(1))*(2*a-1)*p) - k*(emax*(xb/V))/(ec50 + xb/V))"
		Mayer = "fabs(10000000 * math.exp(((2*a-1)*p - k*(emax*(xb/V))/(ec50 + xb/V))*T) - 1000)"
		#Mayer = "fabs(10000000 * math.exp(((2*a-1)*p - k*emax*xb*(V**(-1))*((ec50 + xb*(V**(-1)))**(-1)))) - 1000)"
		Mayer = "0"
	else:
		# Null cost
		Lagrange = "0"
		Mayer = "0"
	return Lagrange, Mayer
################################################################################


