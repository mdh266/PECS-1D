#!/usr/bin/env/python

from numpy import * 
import sys, os, commands
import matplotlib.pyplot as plt

def DopingImage(input, output):

	# open the file
	f = open(input)
	linesInFile = f.readlines()
	f.close()

	# make numpy arrays
	xVals = zeros(len(linesInFile))
	CVals = zeros(len(linesInFile))
	i = 0

	for lines in linesInFile:
		xVals[i] = float( (lines.split())[0])
		CVals[i] = float( (lines.splot())[1])



def CarrierImage(testNum):
	#input file name
	fileName = 'testCarrier' + str(testNum) + '.dat'

	#output file name
	output = 'testCarrier' + str(testNum) + '.png'
	
	#carrier type
	u = 'u' 

	# open the file  
	f = open(fileName)
	linesInFile = f.readlines()
	f.close()

	# make numpy arays
	xVals1 = zeros(len(linesInFile))
	TESTuVals = zeros(len(linesInFile))
	TESTqVals = zeros(len(linesInFile))
	TRUEuVals = zeros(len(linesInFile))
	TRUEqVals = zeros(len(linesInFile))
	eVals = zeros(len(linesInFile))
	pVals = zeros(len(linesInFile))

	i = 0

	for lines in linesInFile:
		xVals1[i] = float( (lines.split())[0])
		TESTuVals[i] = float( (lines.split())[1])
		TESTqVals[i] = float( (lines.split())[2])
		TRUEuVals[i] = float( (lines.split())[3])
		TRUEqVals[i] = float( (lines.split())[4])
		eVals[i] = float( (lines.split())[5])
		pVals[i] = float( (lines.split())[6])
		i = i + 1


	#plot it!	
	plt.figure()
	plt.title(r'$n_{t} = n_{xx}$')
	plt.subplot(3,2,1)  # number of row, columns, image number
	plt.plot(xVals1,TESTuVals, 'b-')
	plt.title(r'test $' + u + '_{h}(x,t)$')

	plt.subplot(3,2,2) 
	plt.plot(xVals1, TESTqVals, 'r-')
	plt.title(r'test $q_{h} = -\frac{\partial ' + u + '_{h}}{\partial x}$')
	

	plt.subplot(3,2,3)
	plt.plot(xVals1, TRUEuVals, 'b-')
	plt.title(r'true $' + u + '(x,t)$') 

	plt.subplot(3,2,4)
	plt.plot(xVals1, TRUEqVals, 'r-')
	plt.title(r'true $q = -\frac{\partial ' + u + '}{\partial x}$')

	plt.subplot(3,2,5)
	plt.plot(xVals1, eVals, 'g-')
	plt.title(r'$E_{h}(x,t)$')

	plt.subplot(3,2,6)
	plt.plot(xVals1, pVals, 'm-')
	plt.title(r'$\Phi_{h}(x,t)$')

	plt.savefig(output)



def PoissonImage(testNum):
	#input file name
	fileName = 'testPoisson' + str(testNum) + '.dat'

	#output file name
	output = 'testPoisson' + str(testNum) + '.png'


	# open the file  
	f = open(fileName)
	linesInFile = f.readlines()
	f.close()

	# make numpy arays
	xVals = zeros(len(linesInFile))
	TESTEVals = zeros(len(linesInFile))
	TESTPVals = zeros(len(linesInFile))
	TRUEEVals = zeros(len(linesInFile))
	TRUEPVals = zeros(len(linesInFile))

	i = 0

	for lines in linesInFile:
		xVals[i] = float( (lines.split())[0])
		TESTEVals[i] = float( (lines.split())[1])
		TESTPVals[i] = float( (lines.split())[2])
		TRUEEVals[i] = float( (lines.split())[3])
		TRUEPVals[i] = float( (lines.split())[4])
		i = i + 1


	#plot it!	
	plt.figure()
	plt.title(r'$$')

	plt.subplot(2,2,1)  # number of row, columns, image number
	plt.plot(xVals,TESTEVals, 'b-')
	plt.title(r'test $' + 'E(x,t)$')


	plt.subplot(2,2,2)  # number of row, columns, image number
	plt.plot(xVals,TESTPVals, 'r-')
	plt.title(r'test $' + '\Phi(x,t)$')


	plt.subplot(2,2,3)  # number of row, columns, image number
	plt.plot(xVals,TRUEEVals, 'g-')
	plt.title(r'true $' + 'E(x,t)$')


	plt.subplot(2,2,4)  # number of row, columns, image number
	plt.plot(xVals,TRUEPVals, 'm-')
	plt.title(r'true $' + '\Phi(x,t)$')


	plt.savefig(output)



def CarrierMovie(endNum, carrier_type):
  cmd2Convert = 'convert '

  for num in range(0, endNum):

    #input file name
    fileName = 'MovieCarrier' + str(num) + '.dat'
	
    #output file name
    output = 'MovieCarrier' + str(num) + '.png'
    
    #carrier type
    u = str(carrier_type) 

    # open the file  
    f = open(fileName)
    linesInFile = f.readlines()
    f.close()

    # make numpy arays
    xVals1 = zeros(len(linesInFile))
    TESTuVals = zeros(len(linesInFile))
    TESTqVals = zeros(len(linesInFile))
    eVals = zeros(len(linesInFile))
    pVals = zeros(len(linesInFile))

    i = 0

    for lines in linesInFile:
	xVals1[i] = float( (lines.split())[0])
	TESTuVals[i] = float( (lines.split())[1])
	TESTqVals[i] = float( (lines.split())[2])
	eVals[i] = float( (lines.split())[3])
	pVals[i] = float( (lines.split())[4])
	i = i + 1


    #plot it!	
    plt.figure()
    plt.subplot(2,2,1)  # number of row, columns, image number
    plt.plot(xVals1,TESTuVals, 'b-')

    plt.subplot(2,2,2) 
    plt.plot(xVals1, TESTqVals, 'r-')
	
    plt.subplot(2,2,3)
    plt.plot(xVals1, eVals, 'g-')

    plt.subplot(3,2,6)
    plt.plot(xVals1, pVals, 'm-')

    plt.savefig(output)
    cmd2Convert += output + ' '

  cmd2Convert += 'Carrier.mp4'

  failure, result = commands.getstatusoutput(cmd2Convert) 


def PoissonCoupledImage(input1, input2, output):

	# open the file  
	f = open(input1)
	linesInFile = f.readlines()
	f.close()

	# make numpy arrays
	xVals1 = zeros(len(linesInFile))
	ElecVals = zeros(len(linesInFile))
	PotVals = zeros(len(linesInFile))
	i = 0

	for lines in linesInFile:
		xVals1[i] = float( (lines.split())[0])
		ElecVals[i] = float( (lines.split())[1])
		PotVals[i] = float( (lines.split())[2])
		i = i + 1

	# open the file  
	f = open(input2)
	linesInFile = f.readlines()
	f.close()

	# make numpy arrays
	xVals2 = zeros(len(linesInFile))
	CarrierVals = zeros(len(linesInFile))
	i = 0

	for lines in linesInFile:
		xVals2[i] = float( (lines.split())[0])
		CarrierVals[i] = float( (lines.split())[1])
		i = i + 1


	#plot it!	
	plt.figure()
	plt.subplot(1,3,2)  # number of row, columns, image number
	plt.plot(xVals1,ElecVals, 'b-')
	plt.title(r'$E(x,t)$')
	plt.xlabel('x values')

	plt.subplot(1,3,3) 
	plt.plot(xVals1, PotVals, 'r-')
	plt.title(r'$\Phi(x,t)$')
	plt.xlabel('x values')
	

	plt.subplot(1,3,1)
	plt.plot(xVals2, CarrierVals, 'g-')
	plt.title(r'$C(x) - [n(x,t) - p(x,t)]$') 
	plt.xlabel('x values')
	
	plt.savefig(output)
