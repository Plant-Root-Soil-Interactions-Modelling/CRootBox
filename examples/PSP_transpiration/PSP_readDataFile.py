#PSP_readDataFile.py   
import csv
import numpy as np
             
def scanDataFile(myFile, myDelimiter):
    myReader = csv.reader(open(myFile, "rt"), delimiter=myDelimiter)
    nrCols = 0
    nrRows = 0
    for myRow in myReader:
        if (nrRows == 0):
            nrCols = len(myRow)
        else:
            if (len(myRow) != nrCols): 
                return (nrRows, 0, False)
        nrRows += 1
    
    return (nrRows, nrCols, True)
        
                
def readDataFile(myFile, nrHeaderRows, myDelimiter, printOnScreen):
    nrRows, nrCols, isFileOK = scanDataFile(myFile, myDelimiter)
    if (printOnScreen): print ('nrRows =', nrRows, ' nrCols =', nrCols)
    
    if (isFileOK == False): return (nrRows, False)
    myReader = csv.reader(open(myFile, "rt"), delimiter=myDelimiter)
    
    A = np.zeros(((nrRows-nrHeaderRows), nrCols))
    i = 0
    for myRow in myReader:
        if (printOnScreen): print(myRow)
        
        if (i >= nrHeaderRows):
            for j in range(0, len(myRow)):
                A[(i-nrHeaderRows), j] = float(myRow[j])
        i += 1    
    
    return(A, True)
