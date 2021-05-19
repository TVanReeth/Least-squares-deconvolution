import subprocess as sp
import numpy as np
import matplotlib.pyplot as plt
import math as m
import matplotlib as mpl
import sys
from shutil import move
from os import remove
mpl.rc('font', size=15)

def ReadUserInput():
    """"Reading the content of the file 'UserInput.dat' """
    commentchar1='#'
    commentchar2='!'
    file=open('UserInput.dat', 'r')
    data = [[],[],[]]
    dates = [[],[],[]]
    i=-1
    #-- read until empty line or end of file
    specname = '/'
    varname  = '/'
    wb       = -1.
    we       = -1.
    result   = '//'
    modcount = 0
    
    while 1:
        line = file.readline() 
        i += 1
        #-- special cases: 
        if not line: break  # end of file
        if line.isspace(): 
            continue # blank line and treat next line
        if line[0] == commentchar1: 
            continue # treat next line
        elif line[0] == commentchar2: 
            continue # treat next line
        #-- real info:
        if(i == 0):
          specname = str(line.strip().split()[2])
        if(i == 1):
          varname = str(line.strip().split()[2])
        if(i == 2):
          wb = float(line.strip().split()[2])
        if(i == 3):
          we = float(line.strip().split()[2])
        if(i == 6):
          result = str(line.strip().split()[2])
        if(i == 8):
          nlsd = int(line.strip().split()[4])
        if(i == 9):
          if(nlsd > 1):
            lim = []
            for j in np.arange(nlsd-1):
              lim.append(float(line.strip().split()[3+j]))
          else:
            lim = []
        if((i >= 11) & (i%4 == 3)):
          modcount += 1
        
    file.close()
    return specname,varname,result,wb,we, modcount, nlsd, lim




    
def PlotResults(specname,result,wb,we,modcount, nlsd, limm):
    plt.figure()
    clrs = ['b','r','y','g','c','m']
    symb = ['-','--','-.']
    
    expl = []
    print(limm)
    limm = [0.] + list(limm) + [1.]
    lim = np.array(limm)[::-1]
    print(lim)
    for i in np.arange(1,modcount+1):
      lsd = np.loadtxt(result+'comp'+str(i)+'.lsd')
      if(i == 1):
          for j in np.arange(1,nlsd+1):
            fmtstr = f'{clrs[j-1]}{symb[j-1]}'
            plt.plot(lsd[:,0],lsd[:,j],fmtstr)
            expl.append('flux at line centre = '+str(int(50.*(float(lim[j-1])+float(lim[j]))))+'%')
          if(nlsd > 1):
            plt.legend(tuple(expl),loc='upper right')
      for j in np.arange(1,1+nlsd):
          fmtstr = f'{clrs[j-1]}{symb[j-1]}'
          plt.errorbar(lsd[:,0],lsd[:,j],yerr=lsd[:,j+nlsd],fmt=fmtstr)
          plt.plot(lsd[:,0],lsd[:,j],fmtstr)
    plt.ylabel('LSD')
    plt.xlabel('Doppler velocity [km/s]')
    
    plt.title(specname)





if __name__ == "__main__":
  
  if(('-h' in sys.argv) | ('--help' in sys.argv)):
      print(' ')
      print('Code for LSD-based analysis of stellar spectra')
      print('----------------------------------------------')
      print('This code computes the least-squares deconvolution profile(s) for the')
      print('spectrum and mask(s) given in the file "UserInput.dat", using the')
      print('matrix method. It can also be used to compute an LSD-based fit to the')
      print('observed spectrum. For more information on the used algorithms, see:')
      print('Donati et al. (1997), MNRAS 291, 658.')
      print('Kochukhov et al. (2010), A&A 524, A5.')
      print('Tkachenko et al. (2013), A&A, in press.')
      print(' ')
      print('Usage: python LSDProject.py [-h] [-c] [-m] [-p]')
      print(' ')
      print('-h, --help            show this help message and exit')
      print('-c, --compile         (re)compile the underlying Fortran 90/95 routines')
      print('-p, --plot            plot the results. When also computing a fit to')
      print('                      the observed spectrum, this includes a comparison')
      print('                      of the observed spectrum and the fit, and a plot')
      print('                      of the residuals.')
      print(' ')
      
      sys.exit()
      
      
      
  if(('-c' in sys.argv) | ('--compile' in sys.argv)):
      comp = True
  else:
      comp = False
  
  if(('-p' in sys.argv) | ('--plot' in sys.argv)):
      plot = True
  else:
      plot = False
  
  ### COMPILING NECESSARY?
  if(comp):
      sp.call("gfortran -O3 -std=legacy ./Source/LSD/*.f90 -o LSDCalc -lblas -llapack", shell=True)
      sys.exit()
  
  ### READING THE INPUT FILE
  specname,varname,result,wb,we,modcount, nlsd, lim = ReadUserInput()
  
  specname = './Data/Spectra'+specname
  varname  = './Data/Spectra'+varname
  
  if((wb != -1.) & (we != -1.)):
      limit = True
  else:
      limit = False
  
  if(result == '//'):
      result = './Data/Results/'
  elif(result[-1] == '/'):
      result = './Data/Results'+result
  else:
      result = './Data/Results'+result+'/'
  
  
  ### RUNNING THE FORTRAN CODE
  print("\nCOMPUTING THE LSD PROFILES...")
  sp.call("./LSDCalc")
  
  ### PLOTTING REQUESTED?
  if(plot):
      PlotResults(specname,result,wb,we,modcount, nlsd, lim)
      plt.show()
