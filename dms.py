#!/usr/bin/env python
# coding: utf-8

################################################################################
# Direct Multiple Shooting module: DMS methods for a disease
# Description: Implementation of direct multiple shooting for continuous and 
#              discrete optimal dosage problem                                                             
################################################################################



################################### Imports ####################################
from casadi import *
import matplotlib.pyplot as plt
import math
import functionalcost as fc
################################################################################



################################## Functions ###################################

#-------------------------------------------------------------------------------
# Direct multiple shooting to resolve the discrete case of the problem
# OBS: F encodes dynamics of the problem (ODE) and cost J
# OBS: volume is optional depending on the disease
def solve(time,
          intervalBetweenDoses, 
          n_doses,
          numberOfControlIntervals, 
          rungeKuttaSteps,
          F,
          Mayer,
          lb_w,
          ub_w,
          volume=0):
	
	# Initialize an empty NLP
	w=[] #vector of controls and state over time (u1,x1,u2,x2,...)
	w0 = [] #initial vector of w
	lbw = [] #lower bound for w
	ubw = [] #upper bound for w
	J = 0 #propagation of L
	g=[] #function on w (evolution over time)
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
	for k in range(numberOfControlIntervals):
		# Obtain time equivalent to k-interval
		t = (time/numberOfControlIntervals)*k
		#Uk = 0 #Declaration of Uk
		if(n_doses>0):
			# Discrete Dosage
			# New NLP variable for the control
			Uk = MX.sym('U_' + str(k), n_doses)
			w += [Uk]
			for i in range(n_doses):
				w0 += [0]
				lbw += [lb_w]
				if(math.ceil(t) == i*intervalBetweenDoses or math.floor(t) == i*intervalBetweenDoses):
					ubw += [ub_w] #constraint control to be under lethal dosage
				else:
					ubw += [0] #constraint control to be 0 out of dosage time
		else:
			# Continuos Dosage
			Uk = MX.sym('U_' + str(k))
			w += [Uk]
			w0 +=[0]
			lbw += [lb_w]
			ubw += [ub_w]
		# Integrate till the end of the interval
		Fk = F(x0=Xk, p=Uk, t=t)
		Xk_end = Fk['xf']
		J = J + Fk['qf'] #Add Lagrange term on evolution of cost J
		# New NLP variable for state at end of interval
		Xk = MX.sym('X_'+str(k+1))
		w += [Xk]
		w0 += [0]
		lbw += [0]
		ubw += [inf]
		# Add equality constraint
		g += [Xk_end-Xk] #only of the interval k
		lbg += [0]
		ubg += [0]
	# Add mayer term
	Mk = Mayer(xf=Xk_end)
	J = J + Mk['mf'] #Add Mayer term to the cost J

	# Define a NLP solver and solve the problem.
	prob = {'f': J, 'x': vertcat(*w), 'g': vertcat(*g)}
	solver = nlpsol('solver', 'ipopt', prob, {'ipopt':{'print_frequency_time':60.0}});

	# Solve the NLP
	sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)
	w_opt = sol['x'].full().flatten()
	if(n_doses == 0):
		# State Xb values
		if(volume != 0):
			# volume as scale factor of the state value
			x_opt = w_opt[0::2]/volume
		else:
			x_opt = w_opt[0::2]
		# Controls d values
		d_opt = [w_opt[1::2]]
	else:
		# State Xb values
		if(volume != 0):
			x_opt = w_opt[0::(n_doses+1)]/volume #start from 0 and select a value every numberOfDoses+1
		else:
			x_opt = w_opt[0::(n_doses+1)]
		# Controls d values
		d_opt = []
		for i in range(n_doses):
			di_opt = w_opt[(i+1)::(n_doses+1)] #start from i+1 and select a value every numberOfDoses+1
			d_opt += [di_opt]
			
	fig1, fig2, dosestable = export_results(x_opt, d_opt, time, numberOfControlIntervals, n_doses)
	return x_opt, fig1, fig2, dosestable
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Print results
def export_results(x_opt, d_opt, T, K, numberOfDoses): #controlla sti parametri e togli gli inutili
	tgrid = [T/K*k for k in range(K+1)]
	plt.clf()
	plt.ioff()
	plt.figure(1)
	# "discretize" dosage
	d_opt2 = d_opt
	dosestable = []
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
					d_opt2[i][j-1] = summ/count #mean
					dosestable += [[math.floor((j-1)*(T/K)), np.round(d_opt2[i][j-1])]]
					count = 0
					summ = 0
	else: 
		for j in range(len(d_opt2[0])):
			dosestable += [[np.round(j*(T/K),3), np.round(d_opt[0][j])]]
		numberOfDoses = 1 
	for i in range(numberOfDoses):
		plt.step(tgrid, vertcat(DM.nan(1), d_opt2[i]), color='blue')
	plt.title('Dosage')
	plt.xlabel('time')
	plt.ylabel('value')
	plt.grid()
	fig1 = plt.gcf()
	fig1.set_size_inches(5,4)
	plt.close(1)
	# State Xb
	plt.figure(2)
	plt.plot(tgrid, x_opt, '--', color='red')
	#plt.axhline(ctg, color="gray", label='target') #target
	plt.title('State Xb')
	plt.xlabel('time')
	plt.ylabel('value')
	plt.legend(['xb','target'])
	plt.grid()
	fig2 = plt.gcf()
	fig2.set_size_inches(5,4)
	plt.close(2)
	
	return fig1, fig2, dosestable
#-------------------------------------------------------------------------------

################################################################################
