# This script analyzes the output of "3atomgate_evolve_ramp.py"
# It calcualtes the reduced density matrix of the qubits, 
# the probability of success and the entanglement between the qubits
# as a time series, and writes them in separate files
# additionally, it also plots the success probability and entanglement
# and saves it in pdf format
#
# It requires Python 2.7 and QuTip 2.0 or higher


from qutip import *
from scipy import *
from scipy import linalg as LA
from pylab import *

import sys


################################################################################################
# Functions 
################################################################################################


# defining the proper limit of the function f(x) = x * log2(x)
def xlog2x(x):
	if x <= 0:
		return 0
	else:
		return x * log2(x)


# concurrence two qubits
# using the formulas from [Wootters, PRL vol 80, p 2245 (1998)]
def concurrence(rho4by4):
	rho = rho4by4
	flip = tensor(sigmay(), sigmay())
	rhotilde = flip * rho.conj() * flip
	rhorhotilde = (rho * rhotilde).full()
	e_vals, e_vec = LA.eig(rhorhotilde)
	lambdas = sorted(sqrt(e_vals))
	c = real(lambdas[3] - lambdas[2] - lambdas[1] - lambdas[0])
	c = max([0,c])
	return c

# entanglement two qubits
# using the formulas from [Wootters, PRL vol 80, p 2245 (1998)]
def entanglement(rho4by4):
	c = concurrence(rho4by4)
	x = ( 1 + sqrt(1 - c**2) )/2
	E = - xlog2x(x) - xlog2x(1-x)
	return E

	
################################################################################################
# Initialization 
################################################################################################
	
# read the registry number from command line arguments
arguments = sys.argv[1:]
if len(arguments) < 1 :
	print('data series number missing')
	sys.exit()
regstring = arguments[0]
# make a 4-digit string out of the registry number
regstring = regstring.zfill(4)

# Load data from the file picked by the registry number
# time points
tlist = loadtxt('./results/times' + regstring + '.txt')
# full density matrix in QObject format
rholist = qload('./results/result' + regstring)

# determine the maximum photon number in the space of the full density matrix
maxphotonnumber = rholist[0].dims[0][3] - 1
	

################################################################################################
# Quantum objects
################################################################################################

# basis states for a single qubit atoms
qubit_0	= basis(3,0)
qubit_1 = basis(3,1)
qubit_e = basis(3,2)

# operators in the Hilbert space of a single qubit atom
qubit_00 = qubit_0 * qubit_0.dag()
qubit_11 = qubit_1 * qubit_1.dag()
qubit_ee = qubit_e * qubit_e.dag()
qubit_e0 = qubit_e * qubit_0.dag()
qubit_e1 = qubit_e * qubit_1.dag()
qubit_id = qeye(3)


# basis states for the control atom
control_g = basis(3,0)
control_f = basis(3,1)
control_E = basis(3,2)

# operators in the Hilbert space of the control atom
control_gg = control_g * control_g.dag()
control_ff = control_f * control_f.dag()
control_EE = control_E * control_E.dag()
control_Eg = control_E * control_g.dag()
control_Ef = control_E * control_f.dag()
control_id = qeye(3)


# vaccum state for the cavity mode
cavity_0    = basis(maxphotonnumber + 1, 0)

# operators in the Hilbert space of the cavity mode
cavity_a	= destroy(maxphotonnumber + 1)
cavity_ad	= create(maxphotonnumber + 1)
cavity_num  = num(maxphotonnumber + 1)
cavity_id   = qeye(maxphotonnumber + 1)


# full form of operators
EE   = tensor( qubit_id, qubit_id, control_EE, cavity_id )		# |E><E|
gg   = tensor( qubit_id, qubit_id, control_gg, cavity_id )		# |g><g|
ff   = tensor( qubit_id, qubit_id, control_ff, cavity_id )		# |f><f|
ee_1 = tensor( qubit_ee, qubit_id, control_id, cavity_id )		# |e_1><e_1|	(qubit 1)
ee_2 = tensor( qubit_id, qubit_ee, control_id, cavity_id )		# |e_2><e_2|	(qubit 2)
na	 = tensor( qubit_id, qubit_id, control_id, cavity_num)		# a^dagger a

a    = tensor( qubit_id, qubit_id, control_id, cavity_a  )		# a
Ef   = tensor( qubit_id, qubit_id, control_Ef, cavity_id )		# |E><f|
Eg   = tensor( qubit_id, qubit_id, control_Eg, cavity_id )		# |E><g|
e0_1 = tensor( qubit_e0, qubit_id, control_id, cavity_id )		# |e_1><0_1|	(qubit 1)
e1_1 = tensor( qubit_e1, qubit_id, control_id, cavity_id )		# |e_1><1_1|	(qubit 1)
e0_2 = tensor( qubit_id, qubit_e0, control_id, cavity_id )		# |e_2><0_2|	(qubit 2)
e1_2 = tensor( qubit_id, qubit_e1, control_id, cavity_id )		# |e_2><1_2|	(qubit 2)



################################################################################################
# Analyzing the data
################################################################################################

# empty containers for the output

Pg = []					# probabiltiy of state |g>
Pf = []					# probability of state |f>
PE = []					# probability of state |E>
ent = []				# entanglement
photon_number = []		# number of photons in the cavity

r00 = []				# <- entries of the reduced density matrix of the qubits
r01 = []				#
r02 = []				#
r03 = []				#
r10 = []				#
r11 = []				#
r12 = []				#
r13 = []				#
r20 = []				#
r21 = []				#
r22 = []				#
r23 = []				#
r30 = []				#
r31 = []				#
r32 = []				#
r33 = []				# <-


# for every full density matrix in the time series
for rho in rholist:

	# Determine probabilities
	Pf.append( (ff * rho * ff).tr() )		# <f|rho|f> = Tr( |f><f| rho |f><f| )
	PE.append( (EE * rho * EE).tr() )		# <E|rho|E> = Tr( |E><E| rho |E><E| )
	
	# Determining number of photons
	photon_number.append( real((rho * na).tr()))	# Tr( rho * a^dagger a )
	
	# Determining success probability = P(control atom is in state g )
	rho_cond = gg * rho * gg	# conditional density matrix, conditioned on |g>
	prob_g = rho_cond.tr()		# prbability of |g>
	# save the value of Pg
	Pg.append( prob_g )
	
	
	# Determine the value entanglement
	
	# if the probability is too small, enta = 0
	if real(prob_g) < 1e-20:
		ent.append( 0 )
	else:
		# normalize the conditional density matrix
		rho_cond = rho_cond / prob_g
		
		# Partial trace of the system, 
		# keeping the Hilbert space of the qubit atoms, each with all three states (0,1,e)
		rho_qubitatoms = ptrace(rho_cond, [0,1])
		
		# Projecting onto the logical qubit basis of the qubit atoms (00, 01, 10, 11)
		# (this gets rid of the state |e> in both qubit atoms)
		
		# convert the QObject to array
		r = rho_qubitatoms.full()
		
		# select specific elements of the density matrix of the qubits
		# to get the reduced density matrix, which describes the states 00, 01, 10 and 11 only
		rho4by4_full = [[ r[0][0], r[0][1], r[0][3], r[0][4]],
						[ r[1][0], r[1][1], r[1][3], r[1][4]],
						[ r[3][0], r[3][1], r[3][3], r[3][4]],
						[ r[4][0], r[4][1], r[4][3], r[4][4]] ]
	
		# save the elements
		r00.append( r[0][0] )
		r01.append( r[0][1] )
		r02.append( r[0][3] )
		r03.append( r[0][4] )
		r10.append( r[1][0] )
		r11.append( r[1][1] )
		r12.append( r[1][3] )
		r13.append( r[1][4] )
		r20.append( r[3][0] )
		r21.append( r[3][1] )
		r22.append( r[3][3] )
		r23.append( r[3][4] )
		r30.append( r[4][0] )
		r31.append( r[4][1] )
		r32.append( r[4][3] )
		r33.append( r[4][4] )

		# convert the new 4x4 array to QObject
		rhological = Qobj( rho4by4_full, [[2,2],[2,2]] )
		
		# normalize it
		rhological = rhological / rhological.tr()

		# Calculate concurrence and entanglement 
		# from the normalized reduced density matrix of the qubits
		ent.append( entanglement(rhological) )
	

################################################################################################
# Print results to file
################################################################################################

# compile the series of Pg, ent, photon number of Pf to a single array
data = c_[real(tlist), real(Pg), real(ent), real(photon_number), real(Pf)]
# save it to a tab-delimited file "data_####.txt"
savetxt('./results/data' + regstring + '.txt', data, delimiter = '\t')

# flatten the reduced density matrix 
# with real parts of off-diagonal elements in the upper triangle
# and imag parts in the lower triangle
data_rho = c_[ real(tlist),
			   real(r00), real(r01), real(r02), real(r03), 
			   imag(r01), real(r11), real(r12), real(r13),
			   imag(r02), imag(r12), real(r22), real(r23),
			   imag(r03), imag(r13), imag(r23), real(r33),
			   ]
# save it to a tab-delimited file "data_rho_####.txt"
savetxt('./results/data_rho_' + regstring + '.txt', data_rho, delimiter = '\t')


# plot Pg, ent, Pf and PE as a function of time
plot(tlist, Pg)
plot(tlist, ent)
plot(tlist, Pf)
plot(tlist, PE)

# add lables and legend to the plot
xlabel('Time')
ylabel('')
legend(("Pg", "Entanglement", "Pf", "PE"))

# print the plot to a pdf file
savefig('./results/plot' + regstring + '.pdf', bbox_inches = 'tight')


