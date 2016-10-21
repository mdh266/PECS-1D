#!/usr/bin/env python

import sys, os, commands
import matplotlib.pyplot as plt
from numpy import *
from plottingTools import ProgressBar

def makeUnipolarImages(numTimeStamps, VAR, ans, imageType, movieType):
	# Check to make sure varible is acceptable, if not throw exception
	if(VAR == 'u'):
		varLabel = '$10^{17} $ cm$^{-3}$'
		plotTitle = 'Density distribution'
		index = 1
		name = 'density'
		scale = 1.0
	elif(VAR == 'E'):
		varLabel = 'V/cm'
		plotTitle = 'Electric field distribution'
		index = 2
		name = 'ElectricField'
		scale = 1.0 #260.0
	elif(VAR == 'p'):
		varLabel = 'V'
		plotTitle = 'Potential distribution'
		index = 3
		name = 'Potential'
		scale = 1.0 #0.0026
	elif(VAR == 'J'):
		varLabel = 'mA cm$^{-2}$'
		plotTitle = 'Steady State Current distribution'
		index = 4
		scale = 1.0 #0.625
		name = 'Current'
	else:
		raise NameError('Input not valid')


	# index corresponds to which column to use in the data files.  This in turn
	# is the column which corresponds to the variable the user choose.
	

		# COMMAND TO CONVER TO MOVIE
	cmd2convert = 'convert '
		

	
	###############################################################################
	# Open each data file according to the time stamp and find the max/min for 
	# variable being plotted
	###############################################################################

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
		i = 0
		
		# Used to keep track of max and min value of y
		yMinValue = 0.0
		yMaxValue = 0.0
		
		# Get the information line by line
		for lines in linesInFile:
			varVals[i] = float((lines.split())[index]) 
			i = i + 1
      
	    # Find min and max of this file, and set to 
			# global max/min if necessary
			tempMax = scale * varVals.max()
			tempMin = scale * varVals.min()
			if(tempMax > yMaxValue):
				yMaxValue = tempMax
			if(tempMin < yMinValue):
				yMinValue = tempMin
	
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
		varVals = zeros( len(linesInFile) )
		varVals2 = zeros( len(linesInFile) )
		i = 0
		
		# Get the information line by line
		for lines in linesInFile:
			xVals[i] = float((lines.split())[0]) 
			varVals[i] = float((lines.split())[index]) 

			if(VAR=='b'):
				varVals2[i] = -1.0*float((lines.split())[index+1])
			# time = float((lines.split())[7]) # not the best way
			i = i + 1

		varVals = scale * varVals
		varVals2 = scale * varVals2

    # Plot the figure
		plt.clf() # clear window
		plt.plot(xVals,varVals, 'b-', linewidth=2)
#		plt.ylim([-5, 5])
		plt.ylim([yMinValue, (1.1)*yMaxValue])
		plt.xticks(fontsize=20)
		plt.yticks(fontsize=20)
		plt.xlabel('$10^{-4}$ cm', fontsize=20)
		plt.title(plotTitle + ' [' + varLabel + ']', fontsize=20)

		# Save the figure
		if(ans == 'y'):
			nameOfImage =  DirName + '/'+ 'State' + stamp + '.' + imageType
		else:	
			nameOfImage = 'State' + stamp + '.' + imageType
	
		plt.savefig(nameOfImage)  # Save image	
	
		# Add the images to be converted into 
		cmd2convert += nameOfImage + ' ' 

	print "\n\nConverting to movie.  Please wait.\n"
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
	print 'Your movie is complete.  It is titled: ' + nameOfMovie
	

# Erase images
#if(ans == 'n'):
#	failure, ouput = commands.getstatusoutput('rm *.' + imageType)
#if failure:
#	raise NameError('Error removing images.')

# Open the movie
#failure, ouput = commands.getstatusoutput('totem ' + nameOfMovie)
#if failure:
#	raise NameError('Error in launching movie')

###############################################################################
