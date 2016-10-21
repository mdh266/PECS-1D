import sys, os, commands
import fileinput
from numpy import *

def RUN(search_bias, replace_bias):
	# search and replace bias
	searchExp = 'appliedBias\t=\t' + str(search_bias) + '\n'
	replaceExp = 'appliedBias\t=\t' + str(replace_bias)+ '\n'
	
	for line in fileinput.input('ddp-input.ini', inplace=1):
		if searchExp in line:
			line	= line.replace(searchExp,replaceExp)
		sys.stdout.write(line)

	# run the code
	cmd_run_sim =	'./solar_cell_app'

	failure, output = commands.getstatusoutput(cmd_run_sim)
	if failure:
		raise NameError('Failed to find CTU_IMEX, check loaded modules.')

	# make the directory
	Directory_Name = 'bias_0_' + str(int(replace_bias/0.001))
	print Directory_Name
	cmd_make_dir = 'mkdir ' + Directory_Name

	failure, output = commands.getstatusoutput(cmd_make_dir)
	if failure:
		raise nameError('Failed to make the simulations data directory')

	#move the data and input file into the directory
	cmd_move_data = 'mv *.dat ./' + Directory_Name
	failure, output = commands.getstatusoutput(cmd_move_data)
	if failure:
		raise nameError('Failed to move data into directory')

	cmd_copy_input = 'cp ddp-input.ini ./' + Directory_Name
	failure, output = commands.getstatusoutput(cmd_copy_input)
	if failure:
		raise nameError('Failed to move input file into directory')

