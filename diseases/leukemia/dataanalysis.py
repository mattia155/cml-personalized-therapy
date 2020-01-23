#!/usr/bin/env python
# coding: utf-8

################################################################################
# Data analysis module: Data analysis module for a single patient              
# Description: Functions to perform data analysis for leukemia                                               
################################################################################



#################################### Import ####################################
import pandas as pd
import numpy as np
import pwlf
from scipy import stats
import math
import matplotlib.pyplot as plt
import functionalcost as fc
import dms
import diseases.leukemia.PK as pk
################################################################################



################################## Functions ###################################

# Function for data analysis without TB data -----------------------------------
def data_analysis_without_tb(sex, age, weight):
	# Without TB data simply compute Clearance, Volume and constant Ka
	CL, V, Ka, _ = compute_CL_V_Ka_Phi(sex, age, weight)
	return CL, V, Ka
# ------------------------------------------------------------------------------

# Function for data analysis with TB data --------------------------------------
def data_analysis_with_tb(sex, age, weight, tbdata, ec50, sigma, intdoses, cmean, last_point = []):
	# With TB data: 1) compute Clearance, Volume, Phi and constant Ka
	CL, V, Ka, Phi = compute_CL_V_Ka_Phi(sex, age, weight)
	a = 0.87 #literature
	p = 0.45 #literature
	k_param = 0.377 #literature
	emax = 1 #literature
	# 2) fitting TB data points
	if(last_point == []):
		# If first time fit all the data points with a piecewise linear function
		# Import tumor burden data (tbdata is a file)
		tbdata = pd.read_csv(tbdata, header = None)
		# Pre-process tb data (remove null values, change scale, ...)
		times = tbdata.iloc[0,:].tolist()
		datas = tbdata.iloc[1,:].tolist()
		# from months to days
		for i in range(len(times)):
			times[i] = times[i]*30
		# log of data values
		for i in range(len(datas)):
			if(isinstance(datas[i],float) or isinstance(datas[i],int)):
			    datas[i] = math.log10(datas[i])
		# keep only valid values
		times2 = []
		datas2 = []
		for i in range(len(datas)):
			if((isinstance(datas[i], str)==False) and (math.isnan(datas[i])==False)):
				times2 += [times[i]]
				datas2 += [datas[i]]
		times = times2
		datas = datas2
		# Piecewise linear fit (step function)
		piecewise_fit = pwlf.PiecewiseLinFit(times, datas)
		breaks = piecewise_fit.fit(2) # Number of lines to use for the fitting
		break_time = int(breaks[1])
		slope = piecewise_fit.slopes[1] # Slopes of 2nd line
		sigma = 1 - piecewise_fit.r_squared()
		x_hat = np.linspace(np.asarray(times).min(), np.asarray(times).max(), 100)
		y_hat = piecewise_fit.predict(x_hat)
		lambda_patient = slope / np.log10(math.exp(1))
		d_patient = (2*a-1)*p - lambda_patient
		# Compute mean concentration under standard dosage (400mg)
		# OBS: first time mean concentration is unknown
		nDos, F, MY, lb, ub = pk.personalPK(time = 168.,
									        intervalBetweenDoses = intdoses,
								            numberOfControlIntervals = 5000,
									        rungeKuttaSteps = 4,
									        functionalCost = 0, #every cost is ok
									        clearance = CL,
									        volume = V,
									        constant_ka = Ka,
									        constant_dosage = 400) #ignoto nel caso di problema continuo
		xvalues, _, _, _ = dms.solve(time = 168.,
									 intervalBetweenDoses = intdoses,
									 n_doses = nDos,
									 numberOfControlIntervals = 5000,
									 rungeKuttaSteps = 4,
									 F = F,
									 Mayer = MY,
									 lb_w = lb,
									 ub_w = ub,
									 volume = V)
		cmean = np.mean(xvalues[2860:])
		# Compute ec50
		ec50 = cmean * ((k_param*emax / d_patient) - 1)
		last_point = [times[len(times)-1], datas[len(datas)-1]]
	else:
		# tbdata is a couple of values [time, value] 
		# Piecewise linear fiction for all the precedent point has already
		# be done. Fit only the last point with the new point
		times = [last_point[0], tbdata[0]] #times in month
		datas = [last_point[1], math.log10(tbdata[1])]
		# linear fit (2 point)
		slope, intercept, r_squared, p_value, std_err = stats.linregress(times, datas)
		sigma = 1 - (-r_squared) #SE FITTO SOLO 2 PUNTI R^2 SARÀ SEMPRE 1
		lambda_patient = slope / np.log10(math.exp(1))
		d_patient = (2*a-1)*p - lambda_patient
		ec50 = cmean * ((k_param*emax / d_patient) - 1)
		last_point = [times[1], datas[1]]
	return CL, V, Ka, Phi, ec50, sigma, cmean, last_point 
# ------------------------------------------------------------------------------

# Function to computer parameters from demographic data ------------------------
def compute_CL_V_Ka_Phi(sex, age, weight):
	# fixed parameters from literature
	ka = 0.437
	theta_a = 12.8
	theta_b = 258
	theta_1 = 12.7
	theta_2 = 0.8
	theta_3 = -2.1
	theta_4 = 61.0
	bw_mean = 70
	age_mean = 50
	# compute weight and age factor to compute personal CL and V
	weight_factor = (weight - bw_mean) / bw_mean
	age_factor = (age - age_mean) / age_mean
	if(sex == 1):
		# Male
		cl = theta_a + theta_1*weight_factor + theta_2 + theta_3*age_factor + 1
		v = theta_b + theta_4
		phi = 60
	else:
		#Female
		cl = theta_a + theta_1*weight_factor - theta_2 + theta_3*age_factor + 1
		v = theta_b - theta_4
		phi = 75
	return cl, v, ka, phi
# ------------------------------------------------------------------------------

################################################################################
