#!/usr/bin/env python

import sys, os, commands
import matplotlib.pyplot as plt
from numpy import *
import ConfigParser

def makeInterfaceImages(VAR):
 
	imageType =	'pdf'
	movieType	= 'mp4'

	parser = ConfigParser.ConfigParser()
	parser.read('ddp-input')
	input_dictionary = {}

	# Physical
	for name in parser.options('physical'):
		string_value = parser.get('physical',name)
		if('ivp_type' != name):
			value = string_value.split()[0]
			input_dictionary[name] = float(value)
		else:
			ivp_type = string_value.split()[0]
			print ivp_type
	# need to do the split since config parser will read 
	# in the entire line

	# computational
	for name in parser.options('computational'):
		string_value = parser.get('computational',name)
		value = string_value.split()[0]
		if('boundarylayerwidth' == name):
			input_dictionary[name] = float(value)
		elif('timestepfactor' == name):
			input_dictionary[name] = float(value)
		else:
			input_dictionary[name] = int(value)


	numTimeStamps = input_dictionary['numtimestamps']

	# Check to make sure varible is acceptable, if not throw exception
	if(VAR == 'u'):
		varLabel = '$10^{16}$ [cm$^{-3}$]'
		plotTitle = 'Density Distribution'
		index = 1
		name = 'Density'
		scale = 1.0
	elif(VAR == 'E'):
		varLabel = 'Volts/cm'
		plotTitle = 'Electric field distribution'
		index = 3
		name = 'ElectricField'
		scale = 1.0
	elif(VAR == 'p'):
		varLabel = 'Volts'
		plotTitle = 'Potential distribution'
		index = 4
		name = 'Potential'
		scale = 1.0
	elif(VAR == 'J'):
		varLabel = 'mA cm$^{-2}$'
		plotTitle = 'Initial Current Distribution'
		index = 5
		name = 'Current'
		scale = 1.0
	else:
		raise NameError('Input not valid')


	cmd2convert = 'convert ' 

###############################################################################
# Open each data file according to the time stamp and find the max/min for 
# variable being plotted
###############################################################################
	yMinValue = 0.0	
	yMaxValue = 0.0
	
#	print "Finding max and min of ouput:"
	
	for timeStamp in range(0,numTimeStamps):
#		ProgressBar(timeStamp, numTimeStamps)
		
		# Open the appropriate file according to the time stamp
		stamp = str.zfill( str(timeStamp), 4) # convert to 4 characters
		nameOfFile = 'State' + stamp + '.dat'
		f = open(nameOfFile)
	
		# Get all the lines from the file and then chose the file
		linesInFile = f.readlines()
		f.close()

		# Create numpy array
		varVals = zeros( len(linesInFile) )
		varVals2 = zeros( len(linesInFile) )
		i = 0
		
		# Get the information line by line
		for lines in linesInFile:
			varVals[i] = float((lines.split())[index]) 
			if('u' == VAR):
				varVals2[i] = float((lines.split())[index+1]) 
			i = i + 1
                 
	  # Find min and max of this file, and set to 
		# global max/min if necessary
		tempMax1 = varVals.max()
		tempMin1 = varVals.min()
		
		tempMax2 = varVals2.max()
		tempMin2 = varVals2.min()
	
		if(tempMin1 < tempMin2):
			tempMin = tempMin1
		else:
			tempMin = tempMin2

		if(tempMax1 > tempMax2):
			tempMax = tempMax1
		else:
			tempMax = tempMax2

		if(tempMin < yMinValue):
				yMinValue = tempMin
		if(tempMax > yMaxValue):
				yMaxValue = tempMax



###############################################################################
# Open each data file according to the time stamp and createa a plot based off
# of the variable the user chose
###############################################################################
#	print "\n\nMaking images:"	

	for timeStamp in range(0,numTimeStamps):

#		ProgressBar(timeStamp, numTimeStamps)
		
		# Open the appropriate file according to the time stamp
		stamp = str.zfill( str(timeStamp), 4) # convert to 4 characters
		nameOfFile = 'State' + stamp + '.dat'
		f = open(nameOfFile)
	
		# Get all the lines from the file and then chose the file
		linesInFile = f.readlines()
		f.close()

		# Create numpy array
		xVals = zeros( len(linesInFile) )	
		VarVals = zeros( len(linesInFile) )
		Var2Vals = zeros( len(linesInFile) )
		electrons = zeros( len(linesInFile) )
		holes = zeros( len(linesInFile) )
		reductants = zeros( len(linesInFile) )
		oxidants = zeros( len(linesInFile) )
		
		i = 0
	
		# Get the information line by line
		for lines in linesInFile:
			xVals[i] = float((lines.split())[0]) 
			VarVals[i] = float((lines.split())[index])
			if( VAR != 'G'): #needs to be last entry
				Var2Vals[i] = float((lines.split())[index+1])
			i = i + 1

		
		# Get the length of the arrays
		LEN = i
	
		# Plot the figure
		plt.clf() # clear window		if(VAR == 'u'):
		if(VAR == 'u'):
			for j in range(0,LEN):
				if(xVals[j] < 0.0 ):
					electrons[j] = VarVals[j]
					holes[j] = Var2Vals[j] 
					reductants[j] = 0.0
					oxidants[j] = 0.0	
				if(xVals[j] >= 0.0):
					electrons[j] = 0.0
					holes[j] = 0.0 
					reductants[j] = VarVals[j]
					oxidants[j] = Var2Vals[j]
			# end the copying over
		
			plt.plot(xVals,electrons,'b-',label='electrons', linestyle='-', linewidth=2)		
			plt.plot(xVals,holes,'r-',label='holes', linestyle='-', linewidth=2)
			plt.plot(xVals,reductants,'g-',label='reductants', linestyle='-',  linewidth=2)	
			plt.plot(xVals,oxidants,'m-',label='oxidants',linestyle='-', linewidth=2)	
				
			plt.legend(('electrons', 'holes', 'reductants',
				    'oxidants'), loc=2, fontsize=20)
			plt.ylim([yMinValue, (1.1)*yMaxValue])
		else:
			VarVals = scale * VarVals
			plt.plot(xVals, VarVals,'b-', linewidth=2)		
		#	plt.ylim([yMinValue, 50]) #(1.1)*yMaxValue])
		#	time = "t = " + tCurrent + " [ps]"

	#	plt.xlabel('[cm]')
	#	plt.ylabel()
		plt.xticks(fontsize=20)
		plt.yticks(fontsize=20)
		plt.title(plotTitle,fontsize=20)

		nameOfImage = name + stamp + '.' + imageType

		plt.savefig(nameOfImage)  # Save image

		# Add the images to be converted into 
		cmd2convert += nameOfImage + ' ' 


###############################################################################
# Convert all the images into a movie
###############################################################################

	# Add the name of the movie onto the string which will execute convert
	nameOfMovie = name + '.' + movieType
	cmd2convert += nameOfMovie

# Run the command to execute convert
	failure, ouput = commands.getstatusoutput(cmd2convert)
	if failure:
		raise NameError('Error in converting png files into mp4 movie')

# Print the name of the movie

