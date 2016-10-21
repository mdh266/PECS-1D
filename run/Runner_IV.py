# IMPORTS
from Run import *
	

# list of range of bias values to run  solar_cell_app
# NOTE: THE FIRST VOLTAGE VALUE IN THIS LIST SHOULD BE 
# THE VALUE OF "appliedBias" in "ddp-input.ini"
bias_list = [0.0, 0.1, 0.2]

numBias = len(bias_list)
replace_bias = bias_list[0]

for i in range(0,numBias):
	search_bias = replace_bias
	replace_bias = bias_list[i]
	print 'Runnig bias: ' + str(replace_bias)
	RUN(search_bias, replace_bias)






