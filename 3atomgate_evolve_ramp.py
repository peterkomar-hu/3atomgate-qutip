# this script numerically solves the time evolution of the exact Master
# equation of one control atom and two qubit atoms in the cavity
# it samples the time series of the density matrix, and saves it to
# an output file
#
# input parameters are read from "3atomgate_evolove.ini"
# as well as from the command line (which override values from ini file)
# 
# requires Python 2.7 or higher
# requires installation of QuTip 2.0 or higher

from qutip import *
from scipy import *
import sys


################################################################################################
# Functions 
################################################################################################

def d_eff(N, C, DE, De):
	numerator = De * 0.5j - N * C**2
	denominator = De * (DE * 0.5j - C**2) - DE * N * C**2
	return real(numerator / denominator)

def read_Delta_table(filename):
# this function loads Delta1 and Delta2 from files "Delta1(2).txt"

	# empty containers
	header_list = []	# list of header values (values of beta in this script)
	n_header = 0			# its size
	C_list = []			# list of cooperativities
	n_C = 0				# its size
	Delta_numopt = []	# lis of Delta values
	
	# open input file
	Delta_file = open(filename, "r")

	# read the header line
	line = Delta_file.readline()
	
	# split it up along '\t's into words
	word = line.split()
	
	# determine the number of header entries
	n_header = len(word)
	
	# load the numerical values of the header in the numeric container
	for header_string in word:
		header_list.append( float( header_string ) )
		
	# for every line after the header
	line = Delta_file.readline()
	while line != '' :
		# read the first word in every line as a value of C
		word = line.split()
		C_list.append( float( word[0] ) )
		n_C = n_C + 1
		
		# fill the following entries in the line into a numeric array
		numbers = []
		for num in word[1:len(word)]:
			numbers.append( float( num ) )
		
		# append the array of numbers to a list of row vectors
		Delta_numopt.append( numbers )
		
		# read the next line
		line = Delta_file.readline()

	# return the header, the first column (values of C), and the bulk of the table
	return [header_list, C_list, Delta_numopt]

def find_Delta(header, C, filename):
# this function looks up the table from file filename,
# and finds the element in the bulk
# which corresponds to the header entry and the row starting with the value C

	# read the table from file
	[header_list, C_list, D_array] = read_Delta_table(filename)
	
	# find the header entry which is closest to header
	aux_list = []
	for num in header_list:
		aux_list.append( abs( num - header) )
	header_index = argmin(aux_list)
	
	# find the index of row, which starts with the number closest to C
	aux_list = []
	for num in C_list:
		aux_list.append( abs( num - C) )
	C_index = argmin(aux_list)
	
	# retun the looked up element
	return D_array[C_index][header_index]

	
	
	
################################################################################################
# Initialize parameters
################################################################################################


# obtain register number from "reg.txt"
regfile = open("reg.txt","r")
reg = int( regfile.read() )
regfile.close()

# update register file by incrementing the number stored in "red.txt"
regfile = open("reg.txt","w")
regfile.write( str(reg + 1) + "\n" )
regfile.close()

# read parameters from ini file
inifile = open('3atomgate_evolve.ini', 'r')
Delta1_default_th = False
Delta1_default_num = False
Delta2_default_th = False
Delta2_default_num = False
g_default = False
g1_default = False
time_total_default = False
Omega_default = False
gamma_f_default = False
state_default = True
line = inifile.readline()
while line != '' :
	word = line.split()
	if word[0] == 'maxphotonnumber':
		# photon number cutoff
		maxphotonnumber = int(word[1])
	elif word[0] == 'C':
		# cooperativity 
		# (=g^2/(kappa*gamma) =  plausible values [10..1000])
		C = float(word[1])
	elif word[0] == 'kappa':
		# photon loss rate 
		# (=1, setting the time unit)
		kappa   = float(word[1])		
	elif word[0] == 'gamma':
		# qubit  atoms decay rate to state 0
		# (= plausible values {0.01, 0.001}kappa)
		gamma   = float(word[1])		
	elif word[0] == 'gamma_g':
		# decay rates of control atom to g and f states, respectively
		# (= plausible values {0.5, 0.1, 0.01}gamma)
		gamma_g = float(word[1])		
	elif word[0] == 'gamma_f':
		# decay rates of control atom to g and f states, respectively
		# (= plausible value is gamma - gamma_g)
		if word[1] == 'default':
			gamma_f_default = True
		else:
			gamma_f = float(word[1])		
	elif word[0] == 'Omega':
		# laser strength driving the control atom 
		# ( << gamma < kappa, g, weak coupling limit to satisfy adiabaticity)
		if word[1] == 'default':
			Omega_default = True
		else:
			Omega = float(word[1])		
	elif word[0] == 'Delta1':
		# g-E transition's detuning in control atom
		# (= the optimal value for balancing out the decay rates of 00 and 11, by default)
		if word[1] == 'default_th':
			Delta1_default_th = True
		elif word[1] == 'default_num':
			Delta1_default_num = True
		else:
			Delta1 = float(word[1])		
	elif word[0] == 'Delta2':
		# 0-e transition's detuning in qubit atoms
		# (= the optimal value for balancing out the decay rates of 00 and 11, by default)
		if word[1] == 'default_th':
			Delta2_default_th = True
		if word[1] == 'default_num':
			Delta2_default_num = True
		else:
			Delta2 = float(word[1])		
	elif word[0] == 'delta':
		# g-(E)-f Raman transitions detuning
		# (=0, for resonant Raman transition)
		delta = float(word[1])
	elif word[0] == 'g':
		# qubit-cavity coupling for control atom
		# (= sqrt(C * kappa * gamma), from defition of cooperativity, by default)
		if word[1] == 'default':
			g_default = True
		else:
			g = float(word[1])		
	elif word[0] ==	'g1':
		# control atom - cavity coupting for qubit atoms
		# (= g, in case of similar atomic structures, for simplicity, by default)
		if word[1] == 'default':
			g1_default = True
		else:
			g1 = float(word[1])		
	elif word[0] ==	'time_total':
		# total simulation time
		# (= theoretical gate time, by default)
		if word[1] == 'default':
			time_total_default = True
		else:
			time_total = float(word[1])		
	elif word[0] == 'N_timepoints':
		# total number of sampling timepoints
		N_timepoints = int(word[1])
	elif word[0] ==	'beta':
		# coefficient of Omega: Omega = beta * gamma * sqrt(C)
		beta = float(word[1])		
	elif word[0] ==	'tau':
		# timescale of ramping up and down of the driving field
		tau = float(word[1])		
	else :
		print('unrecognized ini file line:' + line)
	line = inifile.readline()
inifile.close()

# read command line parameters
arguments = sys.argv[1:]
i = 0
while i < len(arguments):
	if arguments[i] == 'maxphotonnumber':
		# photon number cutoff
		maxphotonnumber = int(arguments[i+1])
	elif arguments[i] == 'C':
		# cooperativity 
		# (=g^2/(kappa*gamma) =  plausible values [10..1000])
		C = float(arguments[i+1])
	elif arguments[i] == 'kappa':
		# photon loss rate 
		# (=1, setting the time unit)
		kappa   = float(arguments[i+1])		
	elif arguments[i] == 'gamma':
		# qubit  atoms decay rate to state 0
		# (= plausible values {0.01, 0.001}kappa)
		gamma   = float(arguments[i+1])		
	elif arguments[i] == 'gamma_g':
		# decay rates of control atom to g and f states, respectively
		# (= plausible values {0.5, 0.1, 0.01}gamma)
		gamma_g = float(arguments[i+1])		
	elif arguments[i] == 'gamma_f':
		# decay rates of control atom to g and f states, respectively
		# (= plausible value is gamma - gamma_g)
		if arguments[i+1] == 'default':
			gamma_f_default = True
		else:
			gamma_f = float(arguments[i+1])
			gamma_f_default = False			
	elif arguments[i] == 'Omega':
		# laser strength driving the control atom 
		# ( << gamma < kappa, g, weak coupling limit to satisfy adiabaticity)
		if arguments[i+1] == 'default':
			Omega_default = True
		else:
			Omega = float(arguments[i+1])
			Omega_default = False			
	elif arguments[i] == 'Delta1':
		# g-E transition's detuning for control atom
		# (= the optimal value for balancing out the decay rates of 00 and 11, by default)
		if arguments[i+1] == 'default_th':
			Delta1_default_th = True
		elif arguments[i+1] == 'default_num':
			Delta1_default_num = True
		else:
			Delta1 = float(arguments[i+1])
			Delta1_default_th = False
			Delta1_default_num = False
	elif arguments[i] == 'Delta2':
		# 0-e transition's detuning for qubit atom
		# (= the optimal value for balancing out the decay rates of 00 and 11, by default)
		if arguments[i+1] == 'default_th':
			Delta2_default_th = True
		elif arguments[i+1] == 'default_num':
			Delta2_default_num = True
		else:
			Delta2 = float(arguments[i+1])
			Delta2_default_th = False
			Delta2_default_num = False
	elif arguments[i] == 'delta':
		# g-E-f Raman transitions detuning
		# (=0, for resonant Raman transition)
		delta = float(arguments[i+1])
	elif arguments[i] == 'g':
		# qubit-cavity coupling
		# (= sqrt(C * kappa * gamma), from defition of cooperativity, by default)
		if arguments[i+1] == 'default':
			g_default = True
		else:
			g = float(arguments[i+1])
			g_default = False			
	elif arguments[i] ==	'g1':
		# control atom - cavity coupting
		# (= g, in case of similar atomic structures, for simplicity, by default)
		if arguments[i+1] == 'default':
			g1_default = True
		else:
			g1 = float(arguments[i+1])
			g1_default = False
	elif arguments[i] ==	'time_total':
		# total time of the simulation
		# (= theoretical gate time, by default)
		if arguments[i+1] == 'default':
			time_total_default = True
		else:
			time_total = float(arguments[i+1])
			time_total_default = False
	elif arguments[i] == 'N_timepoints':
		# number of sampling timepoints
		N_timepoints = int(arguments[i+1])
	elif arguments[i] == 'beta':
		# coefficient of Omega: Omega = beta * gamma * sqrt(C)
		beta = float(arguments[i+1])
	elif arguments[i] == 'tau':
		# timescale of ramping up and down of the driving field
		tau = float(arguments[i+1])
	elif arguments[i] == 'state':
		# index of the state of the system to start the evolution from
		state = int(arguments[i+1])
		state_default = False
	else:
		print('unrecognized command line argument: ' + arguments[i])
		sys.exit(2)
	i = i + 2
	
# set values of parameters which are set to 'default' (or 'defalut_th'/'default'_num')

# driving strength
if Omega_default:
	Omega = beta * gamma * sqrt(C)
	
# decay rate for E->f
if gamma_f_default:
	gamma_f = max([0, gamma - gamma_g])

# control atom's detuning
if Delta1_default_th:
	Delta1 = gamma/2 * sqrt( 4*C + 1 )
elif Delta1_default_num:
	Delta1 = find_Delta(beta, C, "Delta1.txt")

# qubit atoms' detuning
if Delta2_default_th:
	Delta2 = C * gamma**2 / (2 * Delta1)
elif Delta2_default_num:
	Delta2 = find_Delta(beta, C, "Delta2.txt")

# control atom's coupling to the cavity field
if g_default:
	g  = sqrt(C * kappa * gamma)
	
# qubit atoms' coupling to the cavity field
if g1_default:
	g1 = g

# the total time is determined by the expected length of the pulse * 1.1
# to cover the slow evolution after ramping down the drive to a minuscule value
if time_total_default:
	time_total = tau * 1.1
print('time_total = ')
print(time_total)


################################################################################################
# Initialize quantum objects
################################################################################################


#### Operators ####

# basis states for qubit atoms (three states are labelled: '0','1','e')
qubit_0	= basis(3,0)
qubit_1 = basis(3,1)
qubit_e = basis(3,2)

# operators for probabilities and coherences of a single qubit atom
qubit_00 = qubit_0 * qubit_0.dag()
qubit_11 = qubit_1 * qubit_1.dag()
qubit_ee = qubit_e * qubit_e.dag()
qubit_e0 = qubit_e * qubit_0.dag()
qubit_e1 = qubit_e * qubit_1.dag()

# identitiy in the Hilbert space of a single qubit atom
qubit_id = qeye(3)


# basis states for control atom (three states are labelled: 'g', 'f', 'E')
control_g = basis(3,0)
control_f = basis(3,1)
control_E = basis(3,2)

# operators from probabilities and coherences for the control atom
control_gg = control_g * control_g.dag()
control_ff = control_f * control_f.dag()
control_EE = control_E * control_E.dag()
control_Eg = control_E * control_g.dag()
control_Ef = control_E * control_f.dag()

# identity in the Hilbert space of  the control atom
control_id = qeye(3)


# vacuum state for cavity mode
cavity_0    = basis(maxphotonnumber + 1, 0)

# creation and anihilation operators for cavity
cavity_a	= destroy(maxphotonnumber + 1)
cavity_ad	= create(maxphotonnumber + 1)

# number operator for cavity
cavity_num  = num(maxphotonnumber + 1)

# identity for cavity
cavity_id   = qeye(maxphotonnumber + 1)


# full form of operators
EE   = tensor( qubit_id, qubit_id, control_EE, cavity_id )		# |E><E|
gg   = tensor( qubit_id, qubit_id, control_gg, cavity_id )		# |g><g|
ee_1 = tensor( qubit_ee, qubit_id, control_id, cavity_id )		# |e_1><e_1|	(qubit atom 1)
ee_2 = tensor( qubit_id, qubit_ee, control_id, cavity_id )		# |e_2><e_2|	(qubit atom 2)
na	 = tensor( qubit_id, qubit_id, control_id, cavity_num)		# a^dagger a

a    = tensor( qubit_id, qubit_id, control_id, cavity_a  )		# a
Ef   = tensor( qubit_id, qubit_id, control_Ef, cavity_id )		# |E><f|
Eg   = tensor( qubit_id, qubit_id, control_Eg, cavity_id )		# |E><g|
e0_1 = tensor( qubit_e0, qubit_id, control_id, cavity_id )		# |e_1><0_1|	(qubit atom 1)
e1_1 = tensor( qubit_e1, qubit_id, control_id, cavity_id )		# |e_1><1_1|	(qubit atom 1)
e0_2 = tensor( qubit_id, qubit_e0, control_id, cavity_id )		# |e_2><0_2|	(qubit atom 2)
e1_2 = tensor( qubit_id, qubit_e1, control_id, cavity_id )		# |e_2><1_2|	(qubit atom 2)


#### the Hamiltonian ####

# diagonal terms
H_qubit1   = Delta2 * ee_1
H_qubit2   = Delta2 * ee_2
H_control  = Delta1 * EE
H_cavity = delta * na
H_0 = H_qubit1 + H_qubit2 + H_control + H_cavity

# interaction terms
H_int_control = g1 * (a * Ef  +  a.dag() * Ef.dag())
H_int_qubit1   = g  * (a * e1_1  +  a.dag() * e1_1.dag())
H_int_qubit2   = g  * (a * e1_2  +  a.dag() * e1_2.dag())
H_int = H_int_control + H_int_qubit1 + H_int_qubit2

# driving term
H_drive = 0.5*Omega * ( Eg + Eg.dag() )



#### Ramping function ####

# this realizes an approximate sine curve, which does not go below 0,
# but continuously changes to a small constant value at the end
buff = 0.99;
def ramp(t,args):
	return sqrt( ( sin( pi * t / (buff*tau) ) )**2 * ( 0.5 * (sign(tau * buff - t) + 1) )  + (1 - buff)/2 )

# total Hamiltonian with both the constant and time-dependent parts	
# to be used in the numeric solver later
H = [H_0, H_int, [H_drive, ramp]]



#### Lindblad operators #####

L_photonloss = sqrt(kappa) * a					# cavity decay
L_control_to_g = sqrt(gamma_g) * Eg.dag()		# control atom's E -> g decay
L_control_to_f = sqrt(gamma_f) * Ef.dag()		# control atom's E -> f decay
L_qubit1 = sqrt(gamma) * e0_1.dag()				# qubit atom 1's decay
L_qubit2 = sqrt(gamma) * e0_2.dag()				# qubit atom 2's decay

# list of jump operators 
# to be used in the numeric solver later
jump_operators = [L_photonloss, L_control_to_g, L_control_to_f, L_qubit1, L_qubit2]


print('initialization done\n')


################################################################################################
# Solving the master equation
################################################################################################

# initial state
#  = unentangled product (0 + 1) * (0 + 1) * g * 0
rotated_qubit = (qubit_0 + qubit_1)/sqrt(2)

# default starting state
psi_def = tensor( rotated_qubit, rotated_qubit, control_g, cavity_0 )

# alternative starting states (for testing the behavior of single sectors)
psi00 = tensor( qubit_0, qubit_0, control_g, cavity_0)
psi01 = tensor( qubit_0, qubit_1, control_g, cavity_0)
psi10 = tensor( qubit_1, qubit_0, control_g, cavity_0)
psi11 = tensor( qubit_1, qubit_1, control_g, cavity_0)


# choose the starting state
# depdending on the numerical value for 'state' in the command line arguments
if state_default:
	psi0 = psi_def
else:
	if   state == 0:
		psi0 = psi00
	elif state == 1:
		psi0 = psi01
	elif state == 2:
		psi0 = psi10
	elif state == 3:
		psi0 = psi11

# generate the list of timepoints where output is obtained
tlist = linspace(0.0, time_total, N_timepoints)

# ------------------------------------------------------------------
# (possible lines of code to include, to tweak the behavior of the solver)
# set solver options
# opts = qutip.Odeoptions()
# opts.atol = 1e-8
# opts.rtol = 1e-6
# opts.nsteps=1000 
# opts.order=12 
# opts.method='adams'
# opts.first_step=0
# opts.max_step=0
# opts.min_step=0
# opts.tidy=True
# ------------------------------------------------------------------


print('reached solving part\n')

# SOLVE!
# call the solver
psi_start = psi0
result = mesolve(H, psi_start, tlist, jump_operators, [])

print('solving complete\n')

################################################################################################
# Print output to file
################################################################################################

# save the time series of density matrix
states = result.states


# sampling down to a total of ~2000 timepoints
states_sampled = []
times_sampled = []
di = max([1, int(len(states)/2000) ])
i = 0
while i < len(states):
	states_sampled.append( states[i] )
	times_sampled.append( tlist[i] )
	i = i + di

# convert the registry number to a 4-digit string
regstring = '%04d' % reg

# save the sampled timeseries to file in QuTip format
qsave(states_sampled, './results/result' + regstring)

# save timepoints into file
savetxt('./results/times' + regstring + '.txt', times_sampled)

# save parameters to file
paramfile = open('./results/param' + regstring + '.txt', 'w')
paramfile.write('maxphotonnumber ' + str(maxphotonnumber) + '\n' )
paramfile.write('C ' + str(C) + '\n' )
paramfile.write('kappa ' + str(kappa) + '\n' )
paramfile.write('gamma ' + str(gamma) + '\n' )
paramfile.write('gamma_g ' + str(gamma_g) + '\n' )
paramfile.write('gamma_f ' + str(gamma_f) + '\n' )
paramfile.write('Omega ' + str(Omega) + '\n' )
paramfile.write('Delta1 ' + str(Delta1) + '\n' )
paramfile.write('Delta2 ' + str(Delta2) + '\n' )
paramfile.write('delta ' + str(delta) + '\n' )
paramfile.write('g ' + str(g) + '\n' )
paramfile.write('g1 ' + str(g1) + '\n' )
paramfile.write('time_total ' + str(time_total) + '\n' )
paramfile.write('N_timepoints ' + str(N_timepoints) + '\n' )
paramfile.write('beta ' + str(beta) + '\n' )
paramfile.write('tau ' + str(tau) + '\n' )
paramfile.close()

print('all is well\n')
