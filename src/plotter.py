#!/usr/bin/env python

###############################################################################
#  Plotter.py 
#  Plots and makes movies for main.cpp
#
#  Must load python
#
###############################################################################

# IMPORTS
import sys, os, commands
import matplotlib.pyplot as plt
import ConfigParser
from numpy import *

import bipolar
import unipolar
import interface

###############################################################################
# User choices
###############################################################################

# Image type
imageType = 'pdf'

# Movie type
movieType = 'mp4'


###############################################################################
# Simulation information
##############################################################################
# Read the simulation information from the input file "ddp-input"
parser = ConfigParser.ConfigParser()
parser.read('ddp-input.ini')

# NOTE:  The input parser will change everything to lower case.

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




# Note:  Python will uncapitalize everything when scanning in
#xLeftEndPoint = input_dictionary['xleftendpoint']
#xRightEndPoint = input_dictionary['xrightendpoint']
characteristicLength = input_dictionary['characteristiclength']
#xRightEndPoint  = xRightEndPoint * characteristicLength;
#xRightEndPoint  = xRightEndPoint * characteristicLength;
#timeInitial = input_dictionary['timeinitial']
#timeFinal = input_dictionary['timefinal']
maxOrderMX = input_dictionary['maxordermx']
maxOrderDG = input_dictionary['maxorderdg']
numElements = input_dictionary['numelements']
bias = input_dictionary['appliedbias']

# Get the number of time stamps and interface location

if(ivp_type == 'Continuation'):
	numTimeStamps = 2*input_dictionary['numtimestamps']
else:
	numTimeStamps = input_dictionary['numtimestamps']

##############################################################################
# see if user wants to keep the images and movie in a file
##############################################################################
ans = 'n'
#ans = raw_input("Do you want to keep a record of this run? (y or n)\n\n")

# if yes create the directory
#if(ans == 'y'):
#	DirName = raw_input("\nWhat do you want to name the record\n\n")
#	cmd2MakeDir = 'mkdir ' + DirName
#	failure, ouput = commands.getstatusoutput(cmd2MakeDir)
#	if failure:
#		raise NameError('Failed to make the directory.')

##############################################################################
# get variables to be plotted from the user
##############################################################################

model = raw_input("\nPost Processing Visualization For DDP 1.0\n\n"\
		  "\nPlease enter device type: \n\n"\
		  "u : \t Unipolar Device \nb : \t Bipolar Device\n"\
		  "i : \t Semiconductor Electrolyte Interface\n\n")


if( ('b' == model) | ( 'u' == model) ):
  # str to get the input to get variable to plot
  str4input = "\nWhat do you want to plot? Enter: \nu :\t density \n"\
	      "J :\t current\n"\
 	      "E :\t electric field \np :\t potential \n\n" \

elif( (model == 'i') ):
  # str to get the input to get variable to plot
  str4input = "\nWhat do you want to plot? Enter: \n u :\t density \n"\
	      " J :\t current \n"\
 	      " E :\t electric field \n p :\t potential \n\n"\

else:
  raise NameError('No device of that type')


  # Get variable to plot
VAR = raw_input(str4input)


print "\n\nProcessing...\n\n"

if(model == 'b'):
    bipolar.makeBipolarImages(numTimeStamps, VAR, ans, imageType, movieType,
															maxOrderMX, maxOrderDG, numElements, bias) 
if(model == 'u'):
    unipolar.makeUnipolarImages(numTimeStamps,VAR, ans, imageType, movieType)
  
if(model == 'i'):
    interface.makeInterfaceImages(numTimeStamps, VAR, 
			 ans, imageType, movieType, maxOrderMX,
			 maxOrderDG, numElements, bias) 


###############################################################################
