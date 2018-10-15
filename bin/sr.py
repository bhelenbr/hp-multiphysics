#------------------------------------------------------------------------------
#           Name: sr.py
#         Author: Kevin Harris
#  Last Modified: 02/13/04
#    Description: This Python script demonstrates how to perform a
#                 search-and-replace on a file via the command line.
#------------------------------------------------------------------------------

import os
import sys
import fileinput
    
def file_sr(textToSearchFor, textToReplaceWith, searchFile, outputFile = ""):
 
	if (outputFile == searchFile or outputFile == ""):
		outputFile = open( searchFile +".tmp",'w')
		same_flag = 1
	else:
		outputFile = open(outputFile,'w')
		same_flag = 0
		

	for line in fileinput.input(searchFile):
	    outputFile.write( line.replace( textToSearchFor, textToReplaceWith ) )
	    
	outputFile.close()
	
	if (same_flag):
		os.rename(searchFile +".tmp",searchFile)
	

