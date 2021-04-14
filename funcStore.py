import numpy as np
#from math import factorial
import math as math
import matplotlib.pyplot as plt 
import scipy.special as sp
plt.rcParams.update({'font.size': 18})
from math import factorial, pow, sqrt, pi, e
from decimal import Decimal

# In[] Constants of Approximation

#constant μ, provided in Madras et al table 1.1
# mu_2=2.638
mu2=3

#A_lp=10
#constant μ, provided in Madras et al table 1.1
# mu_3=4.684
mu3=5

#Amplitude constant, provided by res paper
#https://polymerdatabase.com/polymer%20physics/SAW.html
# A_2=1.1771
A2=1

#Amplitude constant, provided by res paper
#https://polymerdatabase.com/polymer%20physics/SAW.html
# A_3=1.205
A3=1

nsurf=2*142600

#nsurf=100

lamda=4

B3=1.1213

# In[]
def na(_r):
    return int(nsurf/_r)

def qvol(_r):
    return (A3*mu3)^_r

def nlp(_r,_na,_slp):
    return (_r-(nsurf/_na))/(_slp-1)

def Ell(_r,_na,_slp):
    return nsurf/_na

def nb(_na,n):
    return n-_na

def stirl(x):
    math.floor(x)
    if x == 0:
        return 1
    return x*math.log(x)-(x)+(0.5)*(math.log(2*math.pi*x))+(1/(12*x))




#yes lambda is spelled wrong, as lambda is a cmd in python
def logw(_j,_na,_r,_slp):
    # print("entered logw")
    # print("(nsurf-j*ELL)",(\
    #     nsurf\
    #     -\
    #     _j\
    #     *\
    #     Ell(_r,_na,_slp)\
    #     ))
    # print("s(nsurf-j*ELL)",stirl(\
    #     nsurf\
    #     -\
    #     _j\
    #     *\
    #     Ell(_r,_na,_slp)\
    #     ))
    # print("s(nsurf-(j+1)ELL)",
    # -stirl(
    #     nsurf-(_j+1)*Ell(_r,_na,_slp)
    #     ))
    return stirl(\
        nsurf\
        -\
        _j\
        *\
        Ell(_r,_na,_slp)\
        )\
    -stirl(
        nsurf-(_j+1)*Ell(_r,_na,_slp)
        )\
    +np.log(nsurf)
    +np.log(lamda)\
    -2*np.log(lamda-1)\
    +Ell(_r,_na,_slp)*\
    (\
    np.log(lamda-1)\
    -\
    np.log(nsurf)
    )


        
def ScrSlp(_na,_r,_slp):
    print("na",_na)
    return stirl(Ell(_r,_na,_slp))\
    -\
    stirl(nlp(_r,_na,_slp))\
    -\
    stirl(\
        Ell(_r,_na,_slp) 
        -\
        nlp(_r,_na,_slp)\
        )

def logRvol(_na,_nvol,_n):
    return np.log(nb(_na,_n))\
    -\
    np.log(_nvol)\
    -\
    np.log(qvol(_na))
        
def logRF(_r,_na,_slp):
    s=0
    splus=0
    narange=np.arange(1,_na)
    #if you print youll see arange only gives one less
    naplusone=np.arange(1,_na+1)
    for j in narange:
        # print("entered rf")
        # print("na: ",_na)
        # print("j: ",j)
        # print("logw=", logw(j,_na,_r,_slp))
        s+=logw(j,_na,_r,_slp)
        # print("s",s)
    for j in naplusone:
        # print("entered naplusone")
        # print("j: ",j)
        splus+=logw(j,_na+1,_r,_slp)
    return splus\
        -\
        s\
        -\
        np.log(_na+1)

def logRlp(_na,_r,_slp):
    return (_na+1)*ScrSlp(_na+1,_r,_slp)-_na*ScrSlp(_na,_r,_slp)+(_r/(_slp-1))*(np.log(mu3)+np.log(B3)+3.7*np.log(_slp))

def DeltaQlog(_na,_r,_slp,_nvol,_n):
    return logRF(_r,_na,_slp)+logRlp(_na,_r,_slp)+logRvol(_na,_nvol,_n)