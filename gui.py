#!/usr/bin/env python3.6
# coding: utf-8

################################################################################
# GUI: main script of myTherapy project
# Description: interface of the project
# Author: Mattia Pennati
################################################################################



#################################### Import ####################################
# Function to check if a package is installed or install it otherwise
import subprocess
def import_or_install(package):
	try:
		__import__(package)
	except ImportError:
		subprocess.check_call(["python", '-m', 'pip', 'install', package]) # install pkg

# Install packages		
import_or_install('casadi')
import_or_install('PySimpleGUI') 
import_or_install('math')
import_or_install('matplotlib')
import_or_install('pandas')
import_or_install('numpy')
import_or_install('pwlf')
import_or_install('pymongo[tls]')
import_or_install('datetime')
import_or_install('tkinter')
import_or_install('csv')
import_or_install('scipy')

# Import packages (for this script)		
import PySimpleGUI as sg
import datetime
import csv
import numpy as np
import pymongo
from pymongo import MongoClient
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.ticker import NullFormatter  # useful for `logit` scale
import matplotlib
import matplotlib.pyplot as plt
import tkinter as tk
import dms
import bayesianupdate as bu
matplotlib.use("TkAgg")
################################################################################



################################### Functions ##################################
# Function to draw a fig in a canvas
def draw(canvas, figure, old_canvas, loc=(0,0)):
	if(old_canvas != None): # Empty old canvas (destroy)
		old_canvas.get_tk_widget().destroy()
	# Create a new canvas
	new_canvas = FigureCanvasTkAgg(figure, canvas)
	new_canvas.draw()
	new_canvas.get_tk_widget().pack(side='top', fill='both', expand=1)
	return new_canvas
################################################################################



###################################### GUI #####################################
# Change default theme
sg.change_look_and_feel('Default1')

# Create a form and put all the stuff inside it
# (layout: list of lines filled with element)

# Layout Tab 1: disease parameters area
tab1lt = [  [sg.Text('Insert the parameters of the disease'), ],

			[sg.Text('Select the disease')],
			[sg.Combo(('leukemia (continuos dosage)','leukemia (discrete dosage)','test'),
			'leukemia (discrete dosage)', 
			size=(68,1), 
			key = 'disease',
			change_submits=True)
			],
			
			[sg.Text('Select the drug')],
			[sg.Combo(('imatinib',),
			'imatinib',
			size = (68,1), 
			key = 'drug')
			], 
			
			[sg.Text('Total time (hours)')], 
			[sg.InputText(72, 
			size=(70,1),
			key = 'time')
			], 
			
			[sg.Text('Interval between doses',
			key = 'intdoseslbl')], 
			[sg.InputText(24, 
			size=(70,1),
			key = 'intdoses')
			]
		 ]
		 
# Layout Tab 2: patient parameters area
tab2lt = [  [sg.Text('Insert the parameters of the specific patient'), ],

			[sg.Text('ID (if user already exist)')], 
			[sg.InputText(size=(68,1),
			key = 'id',
			change_submits=True)
			],
			
			[sg.Text('Sex')], 
			[sg.Combo(['male','female'], 
			'male', 
			size=(68,1),
			key = 'sex')
			],
			
			[sg.Text('Age:')], 
			[sg.InputText(50, 
			size=(70,1),
			key = 'age')
			],
			
			[sg.Text('Weight (kg)')], 
			[sg.InputText(70, 
			size=(70,1),
			key = 'weight')
			],
			
			[sg.Text('Ec50 (if is known, otherwise it will be computed)', size=(72,1))],
			[sg.InputText(size=(70,1),
			key= 'ec50')
			], 
			
			[sg.Text('Time (months) from the start of the therapy:',
			key = 'tbtimelabel',
			visible = False)],
			[sg.InputText(size=(70,1),
			key = 'tbtime',
			visible = False)],
			[sg.Text('TB value:',
			key = 'tbvaluelabel',
			visible = False)],
			[sg.InputText(size=(70,1),
			key = 'tbvalue',
			visible = False)],
			[sg.Text('Tumor Burden data (csv)',
			key = 'tblabel')],
			[sg.InputText(size=(58,1),
			key = 'tb'), 
			sg.FileBrowse(file_types=(("CSV Files", "*.csv"),),
			key = 'tbfile')
			]
			
		 ]
		 
# Layout Tab 3: optimization parameters area
tab3lt = [  [sg.Text('Insert the optimization parameter'), ],

			[sg.Text('Functional cost')], 
			[sg.Combo(['optimal target','lower bound target', 'maximization of tumor burden and minimization of toxicity', 'adaptive'], 
			'optimal target', 
			size=(68,1),
			key = 'cost')
			], 
			   
			[sg.Text('Control intervals K (increase accuracy but also execution time)')], # K e M andranno tolti dall'interfaccia 
			[sg.Slider(range=(100,1000),
			default_value=100,
			size=(45,10),
			orientation='horizontal',
			font=('Helvetica', 8),
			key = 'k')
			],
			
			[sg.Text('Runge-Kutta 4th order number of steps')],
			[sg.Slider(range=(1,16),
			default_value=4,
			size=(45,10),
			orientation='horizontal',
			font=('Helvetica', 8),
			key = 'm')
			]
		 ]

# Layout Tab 4: results area	   
tab4lt = [  [sg.Combo(['Optimal dosage', 'State evolution'],
			'Optimal dosage',
			size=(68,1),
			key = 'res',
			change_submits=True)
			],
			
			[sg.Canvas(key='canvas1')],
			
			[sg.Button('Export CSV',
			key="exp",
			size = (70,1),
			visible=False)
			]
		 ]

# Full layout
layout = [  [sg.TabGroup([
			[sg.Tab('Desease parameters', tab1lt), 
			sg.Tab('Patient parameters', tab2lt),
			sg.Tab('Optimization parameters', tab3lt),
			sg.Tab('Results', tab4lt)]
			])
			],
			[sg.Button('Dosage Optimization',
			size = (73,1))
			],
			[sg.Text('Status: Waiting for database connection...',
			size = (40,1),
			key="infolbl")
			]		  
		 ]
		 
# Create the Form
form = sg.FlexForm("Personalized therapy", auto_size_text=True, auto_size_buttons=True)
form.Layout(layout)
form.Finalize()

# Define canvas elements 
canvas_elem1 = form.FindElement('canvas1')
canvas1 = canvas_elem1.TKCanvas
old_canvas = None # to update canvas image
################################################################################



############################## Mongo DB connection ############################# 
# Mongo uri from Atlas DB
MONGO_URI = "mongodb+srv://admin:admin@mytherapy-hsosd.mongodb.net/test?retryWrites=true&w=majority"
# Use MongoClient as client connecting to MONGO_URI (5000 ms before a timeout)
try:
	client = MongoClient(MONGO_URI, serverSelectionTimeoutMS=3000)
	# Connect to specific cluster (database)
	db = client['patients']
	# Check server status
	db.command("serverStatus", serverSelectionTimeoutMS=3000)
	# Connect to a collection (Collection patients in cluster patients)
	patients = db.patients
except Exception as exc:
	# If mongodb is not reachable
	client = None
	sg.Popup("Database not reachable, check your connection and re-start the application. You still can work offline.", title="ERROR")
form.FindElement('infolbl').Update("Status: Ready")
form.Refresh()
################################################################################



################################ Default values ################################
# Default values (for input fields not filled)
T = 24.
p_id = None
old_patient = None
age = -1
weight = -1
ec50 = -1
tbdata = ''
optimizable = True
resultsvisibility = False
prior_ec50_mean = 0.1234
prior_ec50_sigma = 0.2 * 0.1234
tbpoint = []
################################################################################



################################## Main Loop ###################################
# Event Loop to process "events" and get the "values" of the inputs
while True:
	
	# Read the values of the form
	event, values = form.Read()
	
	## Manage input fields -----------------------------------------------------
	
	# Disease and interval between doses
	if(values['disease'] == 'leukemia (continuos dosage)'):
		# Import modules for the selected disease
		import diseases.leukemia.dataanalysis as da
		import diseases.leukemia.PK as pk
		# Update related input fields
		form.FindElement('intdoseslbl').Update(visible=False)
		form.FindElement('intdoses').Update(visible=False)
		form.FindElement('drug').Update(values=('imatinib'))
		intdoses = 0
	elif(values['disease'] == 'leukemia (discrete dosage)'):
		# Import modules for the selected disease
		import diseases.leukemia.dataanalysis as da
		import diseases.leukemia.PK as pk
		# Update related input fields
		form.FindElement('intdoseslbl').Update(visible=True)
		form.FindElement('intdoses').Update(visible=True)
		form.FindElement('drug').Update(values=('imatinib'))
		# Change default values for interval between doses
		if(values['intdoses'] != ''): #if not empty
			try:
				intdoses = int(values['intdoses'])
			except Exception as exc:
				sg.Popup("Interval between doses has to be a number", title = "ERROR")
				optimizable = False
		else:
			form.FindElement('intdoses').Update(24) #Default
			intdoses = 24
	elif(values['disease'] == 'test'):
		# Import modules for the selected disease
		import diseases.test.dataanalysis as da
		import diseases.test.PK as pk
		# Update related input fields
		form.FindElement('intdoseslbl').Update(visible=False)
		form.FindElement('intdoses').Update(visible=False)
		form.FindElement('drug').Update(values=('test1','test2'))
		intdoses = 0
	
	# Drug
	if(values['drug'] == 'imatinib'):
		drug = 0
	elif(values['drug'] == 'test1'):
		drug = 1
	elif(values['drug'] == 'test2'):
		drug = 2
	
	# Total time
	try:
		T = float(values['time'])
	except Exception as exc:
		sg.Popup("Time has to be a number", title = "ERROR")
		optimizable = False
	
	# Functional cost
	if(values['cost'] == 'optimal target'):
		whichcost = 0
	elif(values['cost'] == 'lower bound target'):
		whichcost = 1
	elif(values['cost'] == 'maximization of tumor burden and minimization of toxicity'):
		whichcost = 2
	else:
		whichcost = 3
	
	# Number of control intervals	
	try:
		K = int(values['k'])
	except Exception as exc:
		sg.Popup("Number of control interval has to be a number", title = "ERROR")
		optimizable = False
	
	# Runge-Kutta steps
	try:
		M = int(values['m'])
	except Exception as exc:
		sg.Popup("Number of Runge-Kutta steps has to be a number", title = "ERROR")
		optimizable = False
	
	# Patient id
	if(len(values['id']) >= 5 and client != None):
		p_id = values['id']
		patient = patients.find_one({'_id': p_id})
		if(patient == None):
			sg.Popup("Patient not found, try with a different ID.", title="ERROR")
		elif(p_id != old_patient):
			firstTimePatient = False
			if(firstTimePatient == False):
				sg.Popup("Personal data from the patient imported.", title="INFO")
			form.FindElement('sex').Update(patient['sex'])
			form.FindElement('age').Update(patient['age'])
			form.FindElement('weight').Update(patient['weight'])
			if(len(patient['ec50'])>0):
				form.FindElement('ec50').Update(patient['ec50'][len(patient['ec50'])-1])
			if(len(patient['ec50'])>= 2):
				# Ask for a single tb data value and not for TB data
				form.FindElement('tblabel').Update(visible=False)
				form.FindElement('tb').Update(visible=False)
				form.FindElement('tbfile').Update(visible=False)
				form.FindElement('tbtimelabel').Update(visible=True)
				form.FindElement('tbtime').Update(visible=True)
				form.FindElement('tbvaluelabel').Update(visible=True)
				form.FindElement('tbvalue').Update(visible=True)
				newpoint_time = patient['lastpoint'][0]/30 + 3
				form.FindElement('tbtime').Update(newpoint_time)
				form.Refresh()
			else:
				# Ask for TB data
				form.FindElement('tblabel').Update(visible=True)
				form.FindElement('tb').Update(visible=True)
				form.FindElement('tbfile').Update(visible=True)
				form.FindElement('tbtimelabel').Update(visible=False)
				form.FindElement('tbtime').Update(visible=False)
				form.FindElement('tbvaluelabel').Update(visible=False)
				form.FindElement('tbvalue').Update(visible=False)
				tbpoint = []
				form.Refresh()
			old_patient = p_id
	else:
		p_id = None
		firstTimePatient = True
		# Ask for TB data
		form.FindElement('tblabel').Update(visible=True)
		form.FindElement('tb').Update(visible=True)
		form.FindElement('tbfile').Update(visible=True)
		form.FindElement('tbtimelabel').Update(visible=False)
		form.FindElement('tbtime').Update(visible=False)
		form.FindElement('tbvaluelabel').Update(visible=False)
		form.FindElement('tbvalue').Update(visible=False)
		old_patient = None
		tbpoint = []
		form.Refresh()
			
	# Sex
	if(values['sex'] == 'male'):
		sex = 1
	elif(values['sex'] == 'female'):
		sex = 0
	else:
		sex = -1
	
	# Age
	if(values['age'] != ''):
		try:
			age = int(values['age'])
		except Exception as exc:
			sg.Popup("Age has to be a number", title = "ERROR")
			optimizable = False
		
	# Weight
	if(values['weight'] != ''):
		try:
			weight = float(values['weight'])
		except Exception as exc:
			sg.Popup("Weight has to be a number", title = "ERROR")
			optimizable = False
	
	# Ec50
	if(values['ec50'] != ''):
		try:
			ec50 = float(values['ec50'])
		except Exception as exc:
			sg.Popup("Ec50 has to be a number", title = "ERROR")
			optimizable = False
	else:
		ec50 = -1
	
	# Tb data
	tbdata = values['tb']
	
	# Tb single point
	if(values['tbtime'] != '' and values['tbvalue'] != ''):
		tbpoint = [float(values['tbtime']), float(values['tbvalue'])]
	
	# Results area
	if(resultsvisibility):
		if(values['res'] == 'Optimal dosage'):
			old_canvas = draw(canvas1, fig1, old_canvas)
		else:
			old_canvas = draw(canvas1, fig2, old_canvas)
		 
	##--------------------------------------------------------------------------
	
	
	# If user closes form
	if event is None:	
		break 
 
 
	## Compute therapy: Data Analysis, Optimization ----------------------------
	if event is 'Dosage Optimization':

		# Add new patient on DB
		if(p_id == None and client != None):
			# Get number of patients
			numberofpatients = patients.count()
			# Create id as PTNT + number of new patients
			p_id = 'PTNT' + str(numberofpatients+1)
			# Create patient data
			patient_data = {
				'_id': p_id,
				'sex': values['sex'],
				'age': age,
				'weight': weight,
				'cmean': -1.0, #default
				'ec50': [prior_ec50_mean], #prior
				'sigma': [prior_ec50_sigma], #prior
				'date': datetime.datetime.utcnow() #data
			}
			# Insert patient in the collection
			patients.insert_one(patient_data)
			# Update id field
			form.FindElement('id').Update(p_id)
			sg.Popup("New patient with ID:" + p_id + ". Please store the ID to access to this patient data in the future.", title="INFO")
		#ELSE: RENDERE POSSIBILE L'UPDATE DEI CAMPI ETÀ E PESO
	
		# cost 0 or cost 1 -----------------------------------------------------
		if(whichcost < 2): # cost 0 or cost 1
			# no need of tumor burden data to use optimal and lower bound target
			
			# 1) Cost selection: variable whichcost
			
			# 2) Data analysis
			if(sex != -1 and age != -1 and weight != -1):
				# possible to compute personalized CL, V and Ka
				cl, v, ka  = da.data_analysis_without_tb(sex, age, weight)
			else:
				# not possible to compute personalized CL, V and Ka
				sg.Popup('Some deographic parameter is missing, default will be used (CL=14.3, V=347, ka=0.61)', title='WARNING')
				cl = 14.3
				v = 347
				ka = 0.61
			
			# 3) Optimization
			if(optimizable):
				form.FindElement('exp').Update(visible=False)
				form.FindElement('infolbl').Update("Status: Waiting for optimization end...")
				form.Refresh()
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
				if(values['disease'] == 'leukemia (continuos dosage)' or values['disease'] == 'leukemia (discrete dosage)'):
					# Pass also volume to the solver as scale factor for state
					xval, fig1, fig2, doses = dms.solve(time = T,
												 intervalBetweenDoses = intdoses,
												 n_doses = n_doses,
												 numberOfControlIntervals = K,
												 rungeKuttaSteps = M,
												 F = F,
												 Mayer = My,
												 lb_w = lb,
												 ub_w = ub,
												 volume = v)
				else:
		   			xval, fig1, fig2, doses = dms.solve(time = T,
												 intervalBetweenDoses = intdoses,
												 n_doses = n_doses,
												 numberOfControlIntervals = K,
												 rungeKuttaSteps = M,
												 F = F,
												 Mayer = My,
												 lb_w = lb,
												 ub_w = ub)
				old_canvas = draw(canvas1, fig1, old_canvas)
				resultsvisibility = True
				sg.Popup("Computation end, check the results section.", title="INFO")
				form.FindElement('infolbl').Update("Status: Ready")
				form.FindElement('exp').Update(visible=True)
				form.Refresh()
				if(client != None):
					# Compute mean concentration 
					# exclude first days for a more stable estimation (if possible)
					if(T >= 168.):
						cmean = np.mean(xval[(K/T)*(24*3):])
					else:
						# CASO CHE NON DOVREBBE ESISTERE
						cmean = np.mean(xval)
					patients.find_one_and_update({"_id": p_id}, {"$set": {"cmean": cmean}}) #New cmean after usage of new dosage
				
		# cost 2 or cost 3 -----------------------------------------------------
		else: 
			# need of tumor burden data to use theese costs
		
			# 1) Cost selection: variable whichcost
		
			# 2) Data analysis
			if(sex != -1 and age != -1 and weight != -1 and (tbdata != '' or tbpoint != [])):
				form.FindElement('infolbl').Update("Status: Waiting for optimization end...")
				form.Refresh()
				if(client != None):
					ec50_patient = patients.find_one({"_id": p_id})['ec50']
					sigma_patient = patients.find_one({"_id": p_id})['sigma']
					cmean = patients.find_one({"_id": p_id})['cmean']
				else:
					ec50_patient = [prior_ec50_mean] # Prior
					sigma_patient = [prior_ec50_sigma] # Prior
					cmean = -1
				if(len(ec50_patient)<=1): #Only if connection to db is not available or if is the first time
					# First time (we have only prior, compute first real ec50 and sigma values)
					# First data analysis 
					cl, v, ka, phi, ec50, sigma, cmean, last_point = da.data_analysis_with_tb(sex, age, weight, tbdata, ec50_patient, sigma_patient, intdoses, cmean)
					# Update patient record 
					form.FindElement('ec50').Update(ec50)
					if(client != None):
						patients.find_one_and_update({"_id": p_id}, {"$set": {"lastpoint": last_point}})
						patients.find_one_and_update({"_id": p_id}, {"$set": {"ec50": list(np.append(ec50_patient, ec50))}})
						patients.find_one_and_update({"_id": p_id}, {"$set": {"sigma": list(np.append(sigma_patient, sigma))}})
						patients.find_one_and_update({"_id": p_id}, {"$set": {"cmean": cmean}})
						form.Refresh()
					# Can't do bayesian update
				else: #Only if client is connected
					# At least second time (we have prior and measures)
					# Data analysis
					last_point = patients.find_one({"_id": p_id})['lastpoint']
					# Use optional parameter last_point to fit only the new data point and the last data point
					cl, v, ka, phi, ec50, sigma, cmean, last_point = da.data_analysis_with_tb(sex, age, weight, tbpoint, ec50_patient, sigma_patient, intdoses, cmean, last_point)
					if(len(ec50_patient) > 1): # prior + first measure (+ new measure still not added)
						# We need at least 2 measures and prior to do bayesian update
						# Bayesian update
						ec50, sigma = bu.update(prior_mean=ec50_patient[0], 
						                        prior_sigma=sigma_patient[0],
						                        measure_mean=np.mean(np.append(ec50_patient[1:], ec50)), 
						                        measure_sigma=np.std(np.append(ec50_patient[1:], sigma)), 
						                        n_measures=len(ec50_patient[1:]))
					optimizable2 = True
					# Update patient record 
					form.FindElement('ec50').Update(ec50)
					if(client != None):
						ec50_patient = list(np.append(ec50_patient, [ec50]))
						sigma_patient = list(np.append(sigma_patient, [sigma]))
						patients.find_one_and_update({"_id": p_id}, {"$set": {"ec50": ec50_patient}})
						patients.find_one_and_update({"_id": p_id}, {"$set": {"sigma": sigma_patient}})
						patients.find_one_and_update({"_id": p_id}, {"$set": {"lastpoint": last_point}})
					else: #teoricamente non raggiungibile come caso
						sg.Popup('Connection lost. Can\'t update user data. Restart the application.', title='ERROR')
			else:
				sg.Popup('Fill all patient data to use this cost.', title='ERROR')
				optimizable2 = False
			
			# 3) Optimization
			optimizable2 = False # RIMUOVIMI - NON VOGLIO OTTIMIZZAZIONE
			if(optimizable and optimizable2):
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
				if(values['disease'] == 'leukemia (continuos dosage)' or values['disease'] == 'leukemia (discrete dosage)'):
					# Pass also volume to the solver as scale factor for state
					xval, fig1, fig2, doses = dms.solve(time = T,
														intervalBetweenDoses = intdoses,
														n_doses = n_doses,
														numberOfControlIntervals = K,
														rungeKuttaSteps = M,
														F = F,
														Mayer = My,
														lb_w = lb,
														ub_w = ub,
														volume = v)
				else:
		   			xval, fig1, fig2, doses = dms.solve(time = T,
														intervalBetweenDoses = intdoses,
														n_doses = n_doses,
														numberOfControlIntervals = K,
														rungeKuttaSteps = M,
														F = F,
														Mayer = My,
														lb_w = lb,
														ub_w = ub)
				old_canvas = draw(canvas1, fig1, old_canvas)
				resultsvisibility = True
				form.FindElement('infolbl').Update("Status: Ready")
				sg.Popup("Computation end, check the results section.", title="INFO")
				form.Refresh()
				if(client != None):
					# Compute mean concentration 
					# exclude first days for a more stable estimation (if possible)
					if(T >= 168.):
						cmean = np.mean(xval[(K/T)*(24*3):])
					else:
						# CASO CHE NON DOVREBBE ESISTERE
						cmean = np.mean(xval)
					patients.find_one_and_update({"_id": p_id}, {"$set": {"cmean": cmean}}) #New cmean after usage of new dosage
				
	##--------------------------------------------------------------------------
	
	
	## Export results (if selected)---------------------------------------------
	if event is 'exp':
		filename = str(p_id) + "dosage" + str(datetime.date.today()) + ".csv"
		# ask which file
		filedirectory = sg.PopupGetFolder('Select a folder in which save the result (file:'+filename+')', title="Export CSV")
		# write csv
		filename = filedirectory + "/" + filename
		print(filename)
		with open(filename, "w") as f:
			writer = csv.writer(f)
			writer.writerow(['Time [h]', 'Dosage [mg]']) #header
			writer.writerows(doses)
		sg.Popup("Results exported.", title="INFO", 
				  auto_close=True, auto_close_duration=100,)
	##--------------------------------------------------------------------------
	
form.Close()
################################################################################
