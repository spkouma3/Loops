import numpy as np
#from math import factorial
import math as math
import matplotlib.pyplot as plt 
import scipy.special as sp
plt.rcParams.update({'font.size': 18})

# In[] Constants of Approximation

res=[]

reslog_Long=[]

#constant μ, provided in Madras et al table 1.1
# mu_2=2.638
mu_2=3
print("mu_2=",mu_2)


s_lp=10
A_lp=10
#constant μ, provided in Madras et al table 1.1
# mu_3=4.684
mu3=5
# print("mu_3=",mu_3)

#Amplitude constant, provided by res paper
#https://polymerdatabase.com/polymer%20physics/SAW.html
# A_2=1.1771
A2=1



#Amplitude constant, provided by res paper
#https://polymerdatabase.com/polymer%20physics/SAW.html
# A_3=1.205
A3=1
# print("A_3=",A_3)

# number of short-mers in a large-mer, chain length "r"
# 6200/27=229, but here we want everything in terms of mers
r_Long=6200 



# number of single-mers in a arbitrary mer, chain length "r"
r_Arb=100



# number of single-mers in a arbitrary mer, chain length "r"
r_Short=27
#
#
nsurf=142600

lamda=4

#na is the number of chains so to get a value we must divide nsurf/r and then make it a multiple
#for longs 27 is 4 more than the surface pancake max...
n=27

B3=1.1213


# In[]
def _na(_r):
    return int(nsurf/_r)

def nvol(_na):
    return 1

def qvol(_r):
    return (A3*mu3)^_r

def nlp(_r,_na):
    return (_r-Ell(_na))/(s_lp-1)

def nb(_na):
    return n-_na

def Ell(_na,_r,_slp):
    return _r-nlp(_r,_na)*(_slp-1)

def stirl(x):
    if x<6:
        prod=np.log(sp.factorial(x))
    elif x>5:
        prod=np.log(sp.factorial(5))+((x-10)*np.log(x-10)-(x-10))
    return prod

#yes lambda is spelled wrong, as lambda is a cmd in python
def logwnaplusone(_na):
    return stirl(nsurf-_na*Ell(_na))\
    -stirl(nsurf-(_na+1)*Ell(_na))\
    +np.log(nsurf)
    +np.log(lamda)\
    -2*np.log(lamda-1)\
    +Ell(_na)*\
    (\
    np.log(lamda-1)\
    -\
    np.log(nsurf)
    )
        
def ScrSlp(_r,_na):
    return stirl(Ell(_na))\
    -stirl(nlp(_r,_na))\
        

def logRvol(_na):
    return np.log(nb)\
    -\
    np.log(nvol(_na))\
    -\
    np.log(qvol(_na))
        

def logRF(_na):
    return logwnaplusone(_na)-np.log(_na)

def logRlp(_na,_r,_slp):
    return (_na+1)*\
    ScrSlp(_na+1)\
    -_na*ScrSlp(_na)\
    +(_r/(_slp-1))*\
    (\
    np.log(mu_3)\
    +np.log(B3)\
    +3.7*np.log(_slp)\
    )

def DeltaQlog(_na):
    return logRF(_na)+logRlp(_na)+logRvol(_na)


