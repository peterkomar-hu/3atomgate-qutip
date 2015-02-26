# This script creates summaries of outputs of "3atomgate_analyze.py" and
# "calculate_fidelity_script.m".
# It reads the information about the summaries from "findopt_schedule.txt",
# and creates output files as defined there.
#
# It requires Python 2.7 or higher


from scipy import *

################################################################################################
# Load schedule file 
################################################################################################

# loading schedule file for lists of datasets
# title = label of the list of datasets plotted together (with same gamma, gamma_g)
# numers = regnumbers of the included datasets

# load schedule file
file_schedule = open('findopt_schedule.txt','r')

# empty containers for 
# the list of output filenames
filename = []
# the list of title lines
title = []
# the list of registry numbers
numbers = []

# for every line in the file
line = file_schedule.readline()
while line != ''  :

	# save output filename
	filename.append(line.rstrip())
	# read the next line
	line = file_schedule.readline()
	
	# save title string
	title.append(line.rstrip())
	# read the next line
	line = file_schedule.readline()
	
	# split up the line into words
	words = line.split()
	
	# for every word
	i = 0
	while i < len(words):
		# convert it to a 4-digit numerical string
		words[i] = words[i].zfill(4)
		i = i+1
	# save the 4-digit strings as list of registry numbert
	numbers.append(words)
	
	# read next line
	line = file_schedule.readline()
	
# close file
file_schedule.close()


################################################################################################
# Create summaries
################################################################################################

# for each scheduled job i (from the schedule file)
i = 0
while i < len(title):

	# initialize output file
	file = open('./summary/' + filename[i],'w')
	
	# write header
	file.write(title[i] + '\n')
	file.write('beta\t' + 't_opt\t' + 'Pg_opt\t' + 'ent_opt\t' + 'fidelity_opt' + '\n')
	
	beta = 0
	# for each dataset
	for num in numbers[i]:
	
		# determine beta
		
		# open parameter file (an output of the simulation "3atomgate_evolve_ramp.py")
		paramfile = open('./results/param' + num + '.txt', 'r')
		# read the lines
		line = paramfile.readline()
		while line != '':
			line = line.split()
			
			# if line starts with 'beta', save the value that comes after it as float
			if line[0] == 'beta':
				beta = float(line[1])
				break
			line = paramfile.readline()
		# close parameter file
		paramfile.close()
		
		
		# load the numerical data from 'data####.txt' 
		# (created by "3atomgate_analyze.py")
		data = loadtxt('./results/data' + num + '.txt')
		# save the first column as the timepoints
		t = data[:,0]
		# save the second column as the probabilities of begin in |g>
		Pg = data[:,1]
		# save the third column as the entangelment
		ent = data[:,2]
		
		# load the fidelity timeseries from "fidelity_####.txt" 
		# (created by "calculate_fidelity_script.m")
		data_fidelity = loadtxt('./results/fidelity_' + num + '.txt')
		# save the second column as value of the fidelity
		fidelity = data_fidelity[:,1]
		# (the first column should be identical to the first column of the 'data####.txt'
		#  so it's not loaded here)
		
		# determin the index of the optimal point by maximizing the entanglement
		# (this is maximal at the same point, where the fidelity is maximal)
		optimum = argmax(ent)
		
		# write the record to the output file
		file.write( str(beta) + '\t' + str(t[optimum]) + '\t' + str(Pg[optimum]) + '\t' + str(ent[optimum]) + '\t' + str(fidelity[optimum]) + '\n' )
	
	# close output file
	file.close()
	
	i = i+1
	
	
	
	