#!/usr/bin/env python3.6
# coding: utf-8

################################# Main script ##################################
# Main script of "mytherapy" project.
# 1) Import libraries
# 2) get parameters from terminal
# 3) call every module of the project to solve the problem
# Author: Mattia Pennati
################################################################################ 



# Export results function ------------------------------------------------------
def exportresults(state, doses):
	# come cazzo si stampa?!
	plt.figure(1)
	plt.plot(state)
	plt.title("State Xb")
	plt.grid()
	plt.show()
	filename = "dosage" + str(datetime.date.today()) + ".csv"
	# ask which file
	filedirectory = "/home/mattia/Documents/Università/Tirocinio/mytherapy/"
	# write csv
	filename = filedirectory + filename
	with open(filename, "w") as f:
		writer = csv.writer(f)
		writer.writerow(['Time [h]', 'Dosage [mg]']) #header
		writer.writerows(doses)
#-------------------------------------------------------------------------------



################################## 1) Import ###################################
# Function to check if a package is installed or install it otherwise
import subprocess
def import_or_install(package):
	try:
		__import__(package)
	except ImportError:
		subprocess.check_call(["python", '-m', 'pip', 'install', package]) # install pkg

# Install packages		
import_or_install('casadi') 
import_or_install('math')
import_or_install('matplotlib')
import_or_install('pandas')
import_or_install('numpy')
import_or_install('pwlf')
import_or_install('datetime')
import_or_install('csv')

# Import packages (for this script)		
import datetime
import csv
import sys # get input parameters from console
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import dms
import diseases.leukemia.dataanalysis as da
import diseases.leukemia.PK as pk
################################################################################



################################ 2) Parameters #################################
# Default values (for input fields not filled)
sex = -1
age = -1
weight = -1
ec50 = -1
cmean = -1
tbdata = ''
optimizable = True

# Get parameters from terminal
# interval between doses (0 if continuous dosage)
intdoses = float(sys.argv[1])
# total time T to be optimized
T = float(sys.argv[2])
# which functional cost use (0, 1, 2 or 3)
whichcost = int(sys.argv[3])
# number of control intervals
K = int(sys.argv[4])
# number of runge kutta steps
M = int(sys.argv[5])
# sex of the patient (if known)
if(len(sys.argv)>6):
	sex = int(sys.argv[6])
# age of the patient (if known)
if(len(sys.argv)>7):
	age = int(sys.argv[7])
# weight of the patient (if known)
if(len(sys.argv)>8):
	weight = float(sys.argv[8])	
# tumor burden data
if(len(sys.argv)>9):
	tbdata = sys.argv[9]
# ec50 of the patient (if known)
if(len(sys.argv)>10):
	ec50 = float(sys.argv[10])
################################################################################



############################### 3) Optimization ################################
# cost 0 or cost 1 -------------------------------------------------------------
if(whichcost < 2): # cost 0 or cost 1
	# no need of tumor burden data to use optimal and lower bound target
	# A) Cost selection: variable whichcost
	# B) Data analysis
	if(sex != -1 and age != -1 and weight != -1):
		# possible to compute personalized CL, V and Ka
		cl, v, ka  = da.data_analysis_without_tb(sex, age, weight)
	else:
		# not possible to compute personalized CL, V and Ka
		cl = 14.3
		v = 347
		ka = 0.61
	# C) Optimization
	# Formulate personal PK
	n_doses, F, My, lb, ub = pk.personalPK(time = T,
										   intervalBetweenDoses = intdoses,
										   numberOfControlIntervals = K,
										   rungeKuttaSteps = M,
										   functionalCost = whichcost,
										   clearance = cl,
										   volume = v,
										   constant_ka = ka)
	# Solve the optimization problem through DMS
	# Pass also volume to the solver as scale factor for state
	state, fig1, fig2, doses = dms.solve(time = T,
								 intervalBetweenDoses = intdoses,
								 n_doses = n_doses,
								 numberOfControlIntervals = K,
								 rungeKuttaSteps = M,
								 F = F,
								 Mayer = My,
								 lb_w = lb,
								 ub_w = ub,
								 volume = v)
	exportresults(state, doses)
# cost 2 or cost 3 -------------------------------------------------------------
else: 
	# need of tumor burden data to use theese costs
	# A) Cost selection: variable whichcost
	# B) Data analysis
	if(sex != -1 and age != -1 and weight  != -1 and tbdata != ''):
		cl, v, ka, phi, ec50, cmean = da.data_analysis_with_tb(sex, age, weight, tbdata, ec50, intdoses, cmean)
		# C) Optimization
		# Formulate personal PK
		n_doses, F, My, lb, ub = pk.personalPK(time = T,
											   intervalBetweenDoses = intdoses,
											   numberOfControlIntervals = K,
											   rungeKuttaSteps = M,
											   functionalCost = whichcost,
											   clearance = cl,
											   volume = v,
											   constant_ka = ka,
											   phi = phi,
											   ec50 = ec50)
		# Solve the optimization problem through DMS
		state, fig1, fig2, doses = dms.solve(time = T,
											intervalBetweenDoses = intdoses,
											n_doses = n_doses,
											numberOfControlIntervals = K,
											rungeKuttaSteps = M,
											F = F,
											Mayer = My,
											lb_w = lb,
											ub_w = ub,
											volume = v)
		exportresults(state, doses)
	else:
		print("ERROR: impossible to use this cost without all patient data.")
################################################################################
