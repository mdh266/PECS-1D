#!/usr/bin/env python

import sys, os, commands
import matplotlib.pyplot as plt
from numpy import *
from plottingTools import ProgressBar
#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text', usetex=True)

def makeInterfaceImages(numTimeStamps, VAR, 
												ans, imageType, movieType, maxOrderMX, 
												maxOrderDG, numElements, bias): 

	# Check to make sure varible is acceptable, if not throw exception
	if(VAR == 'u'):
		varLabel = '[1/cm^3]'
		plotTitle = 'Final Density Distribution'
		index = 1
		name = 'density'
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
		varLabel = 'Amps/cm^2'
		plotTitle = 'Current Distribution'
		index = 5
		name = 'Current'
		scale = 1.0
	elif(VAR == 'G'):
		varLabel = 'Amps/cm^2'
		plotTitle = 'Generation Density distribution'
		index = 6
		name = 'Generation'	
		scale = 1.0
	else:
		raise NameError('Input not valid')


	plt.figure(figsize=(10.25,6.5), linewidth=3)

	# index corresponds to which column to use in the data files.  
	# This in turn  is the column which corresponds to the 
	# variable the user choose.

	cmd2convert = 'convert ' 

###############################################################################
# Open each data file according to the time stamp and find the max/min for 
# variable being plotted
###############################################################################
	yMinValue = 0.0	
	yMaxValue = 0.0
		
	print "Finding max and min of ouput:"
	
	for timeStamp in range(0,numTimeStamps):
		ProgressBar(timeStamp, numTimeStamps)
		
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
	print "\n\nMaking images:"	

	for timeStamp in range(0,numTimeStamps):

		ProgressBar(timeStamp, numTimeStamps)
		
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
			#plt.ylim([-1, 40])
		else:
			VarVals = scale * VarVals
			plt.plot(xVals, VarVals,'b-', linewidth=2)		
			plt.ticklabel_format(style='sci', axis='y', scilimits=(-1,1))
		#	plt.ylim([yMinValue, 50]) #(1.1)*yMaxValue])
#			plt.ylim([-3.14159, 3.14159])		
#		plt.ylim([0.158, 0.16])
		if(VAR == 'J'):
			plt.ylim([-0.1, 0.1])
#		time = "t = " + tCurrent + " [ps]"
#		polyOrder = "$P_{MX}$= " + str(maxOrderMX+1) + \
#			    ",\t\t\t\t\t $P_{DG}$ = " + str(maxOrderDG)
#		elements = "Num Elem = " + str(numElements)
#		appliedBias = "App. Bias = " + str(bias) + " [V]"

#		plt.text(xVals[0.05*xVals.size] , 1.40*yMaxValue,time)
#		plt.text(xVals[0.05*xVals.size] , 1.35*yMaxValue, polyOrder)
#		plt.text(xVals[0.05*xVals.size] , 1.30*yMaxValue, elements)
#		plt.text(xVals[0.05*xVals.size] , 1.25*yMaxValue, appliedBias)

#		plt.xlabel('[cm]')
#		plt.ylabe()
		plt.xticks(fontsize=20)
		plt.yticks(fontsize=20)
		plt.title(plotTitle,fontsize=20)

		# Save the figure
		if(ans == 'y'):
			nameOfImage =  DirName + '/'+ 'State' + stamp + '.' + imageType
		else:	
			nameOfImage = 'State' + stamp + '.' + imageType

		plt.savefig(nameOfImage)  # Save image

		# Add the images to be converted into 
		cmd2convert += nameOfImage + ' ' 


###############################################################################
# Convert all the images into a movie
###############################################################################

	print "\n\nConverting to movie.  Please wait.\n"

	# Add the name of the movie onto the string which will execute convert
	nameOfMovie = name + '.' + movieType
	cmd2convert += nameOfMovie

# Run the command to execute convert
	failure, ouput = commands.getstatusoutput(cmd2convert)
	if failure:
		raise NameError('Error in converting png files into mp4 movie')

# Print the name of the movie
	print 'Your movie is complete.  It is titled: ' + nameOfMovie

