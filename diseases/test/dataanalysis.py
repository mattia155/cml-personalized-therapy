#!/usr/bin/env python
# coding: utf-8

################################################################################
# Data analysis module: Data analysis module for a single patient              
# Description: Functions to perform data analysis for test                                               
################################################################################



#################################### Import ####################################
import pandas as pd
import numpy as np
import pwlf
import math
import matplotlib.pyplot as plt
import functionalcost as fc
import dms
import diseases.leukemia.PK as pk
################################################################################



################################## Functions ###################################

# Function for data analysis without TB data -----------------------------------
def data_analysis_without_tb(sex, age, weight):
	CL, V, Ka, Phi = compute_CL_V_Ka_Phi(sex, age, weight)
	return CL, V, Ka
# ------------------------------------------------------------------------------

# Function for data analysis without TB data -----------------------------------
def data_analysis_with_tb(sex, age, weight, tbdata, ec50, intdoses):
	CL, V, Ka, Phi = compute_CL_V_Ka_Phi(sex, age, weight)
	# Import tumor burden data
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
	x_hat = np.linspace(np.asarray(times).min(), np.asarray(times).max(), 100)
	y_hat = piecewise_fit.predict(x_hat)
	plt.figure()
	plt.grid()
	plt.plot(times, datas, 'o', color='r')
	plt.axhline(0,color='black',linewidth=.5)
	plt.xlabel("t[days]")
	plt.ylabel("Log10[BCR-ABL%]")
	plt.plot(x_hat, y_hat, '-', color='b')
	#plt.show()
	
	a = 0.87 #literature
	p = 0.45 #literature
	lambda_patient = slope / np.log10(math.exp(1))
	d_patient = (2*a-1)*p - lambda_patient
	if(ec50 == -1): #if old ec50 is ignote compute it
		ec50 = random()
	#else: IGNOTO - bisogna aggiornare ec50 - feedback
	return CL, V, Ka, Phi, ec50
# ------------------------------------------------------------------------------

# Function to computer parameters from demographic data ------------------------
def compute_CL_V_Ka_Phi(sex, age, weight):
	# fixed parameters from literature
	ka = 0.5
	bw_mean = 75
	age_mean = 35
	# compute weight and age factor to compute personal CL and V
	weight_factor = (weight - bw_mean) / bw_mean
	age_factor = (age - age_mean) / age_mean
	if(sex == 1):
		# Male
		cl = weight_factor*10
		v = weight_factor*100
		phi = 60
	else:
		#Female
		cl = weight_factor*7
		v = weight_factor*70
		phi = 75
	return cl, v, ka, phi
# ------------------------------------------------------------------------------

################################################################################
