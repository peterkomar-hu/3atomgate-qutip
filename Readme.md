This set of codes were used in the research project "3-atom gate", in order to confirm the analytic results of the adiabatic theory by comparing it to exact numerical solutions of the Master equation.

The scripts here should be used in the following way

1. Generate/Update "Delta1.txt" and "Delta2.txt" files
(requires Matlab R2014a or higher)

	- edit C_list and beta_list variables of the "optimize_detunings_script.m" to your desired value
	
	- Run "optimize_detunings_script.m" with Matlab
	
	
2. Define the running parameters
	
	- edit "3atomgate_evolve.ini", set the parameters to the desired values
	
	- for the meaning of 'default', see the code "3atomgate_evolve_ramp.py"

	- for a list of simulations, use the command line arguments, as shown in "run.bat" to set the values of the parameters that change from run to run
	
	- make sure that the values of C and beta match up with the values for which Delta1 and Delta2 have been optimized
	

3. Run the simulation
(requires Python 2.7 and installation of QuTip 2.0 or higher)
	
	- either run "3atomgate_evolve_ramp.py" directly (w/ command line arguments)
	(or call "run.bat" to run a set of simulation defined there)
	
	- note that consecutive runs produce output files labelled by an increasing registry number, which is read from "reg.txt" at the beginning of eahc run
	

4. Analyze output
(requires Python 2.7 and installation of QuTip 2.0 or higher and Matlab R2014a or higher)
	
	- run "3atomgate_analyze.py" with a single command line argument: the registry number of the output to be analyzed. This creates files with time series of the 4x4 density matrix of the qubits, succes probability, and entanglement, as well as plots of the later two.
	(see "run.bat" for sample usage)
	
	- check the plots to confirm that the evolution look as you expected
	
	- modify the list of registy numbers in the beginning of script "calculate_fidelity_script.m", to match the registry numbers of the outputs that have just been analyzed. 
	
	- run "calculate_fidelity_script.m". This creates files containing the time series of fidelity.
	
	
5. Select and summarize a group of outputs
	
	- edit the file "findopt_schedule.txt", such that the first line contains the name of your choice for a new summary file, the second line contains a text to be printed to the first line of the summary file, and the third line contains a series of registry numbers separated by spaces. 
	
	- make sure that each registry number should refer to an already analyzed set of results.
	
	- you can schedule any number of summaries in the file "findopt_schedule.txt" by creating a similar 3 lines right after the first three. Leave no empty lines.
	
	- run "findopt_beta.py". This will create the summary files, containing the first line (as a descriptor), and a table of the maximal fidelity, and the value of entanglement and succes probability at the timepoint of maximal fidelity, for the values of beta
	
	
	
	
	
	