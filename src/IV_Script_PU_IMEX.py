# IMPORTS
import sys, os, commands
import fileinput
from numpy import *
import interface_movies

def RUN(search_bias, replace_bias):

	# search and replace bias
	searchExp = 'appliedBias \t \t \t \t \t = \t ' + str(search_bias) + '\n'
	replaceExp = 'appliedBias \t \t \t \t \t = \t ' + str(replace_bias)+ '\n'
	
	for line in fileinput.input('ddp-input.ini', inplace=1):
		if searchExp in line:
			line	= line.replace(searchExp,replaceExp)
		sys.stdout.write(line)

	# run the code
	cmd_run_sim =	'./PU_IMEX'

	failure, output = commands.getstatusoutput(cmd_run_sim)
	if failure:
		raise NameError('Failed to find PU_IMEX, check loaded modules.')

	# make the directory
	Directory_Name = 'bias_0_' + str(int(replace_bias/0.01))
	print Directory_Name
	cmd_make_dir = 'mkdir ' + Directory_Name

	failure, output = commands.getstatusoutput(cmd_make_dir)
	#if failure:
	#	raise nameError('Failed to make the simulations data directory')

	# MAKE THE IMAGES AND MOVIES
#	VARIABLES = ['u','E','J']
#	for i in range(0,3):
#		VAR = VARIABLES[i]
#		interface_movies.makeInterfaceImages(VAR)


	#move the data and input file into the directory
	cmd_move_data = 'mv *.dat ./' + Directory_Name
	failure, output = commands.getstatusoutput(cmd_move_data)
	#if failure:
	#	raise nameError('Failed to move data into directory')

	cmd_copy_input = 'cp ddp-input.ini ./' + Directory_Name
	failure, output = commands.getstatusoutput(cmd_copy_input)
	#if failure:
	#	raise nameError('Failed to move input file into directory')

###########:#########################################################
bias_list = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.25, 0.35, 0.45, 0.50, 0.55]

numBias = len(bias_list)
replace_bias = bias_list[0]


for i in range(0,numBias):
	search_bias = replace_bias
	replace_bias = bias_list[i]
	print 'Runnig bias: ' + str(replace_bias)
	RUN(search_bias, replace_bias)




