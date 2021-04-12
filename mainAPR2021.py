import numpy as np
#from math import factorial
import math as math
import matplotlib.pyplot as plt 
import scipy.special as sp
plt.rcParams.update({'font.size': 18})
from funcStore import *
#

# In[] Some new things, kind of like a main.

# number of short-mers in a large-mer, chain length "r"
# 6200/27=229, but here we want everything in terms of mers
rLong=6200 



# number of single-mers in a arbitrary mer, chain length "r"
rArb=100



# number of single-mers in a arbitrary mer, chain length "r"
rShort=27

slp=2
#na is the number of chains so to get a value we must divide nsurf/r and then make it a multiple
#for longs 27 is 4 more than the surface pancake max...
# n=27

naNew=na(_r=rLong)+2

#total chains
n=30

nvol=nsurf


dQArr=[]



for k in range(naNew,naNew+4):
    dQArr.append(DeltaQlog(k,_r=rLong,_slp=slp,_nvol=nvol,_n=n))

# In[]
# n_0=n_a(r_Long)+2
# r=np.arange(n_0,n_0+5)
# for j in r:
#         # n_b=n-n_a
#         # print(f_logratio(_n_a,n_b(n,_n_a),n_surf,r_Long,n,n_vol))
#         reslog_Long.append(f_logratio(r_Long,j,n_surf))
#         print(reslog_Long)
#
# In[] Multi plot option
#
#
# plt.figure()   
# x_Long = np.array(r)
# y_Long = reslog_Long
# x_Line = r
# y_Line = np.zeros(r.size)
# # plt.xticks(range(min(x_Long), max(x_Long)+1))
# # # Multiple plots
# plt.xlabel('${n_a}$')
# plt.ylabel('$Log(f(n_a))$')
# plt.grid(b=None, which='major', axis='both')
# # plt.title('Looping Behaviour Long') 
# plt.plot(x_Long, y_Long,linewidth=2) 
# #threshhold
# plt.plot(x_Line, y_Line,'--r',linewidth=2)
#