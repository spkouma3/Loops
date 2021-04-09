import numpy as np
#from math import factorial
import math as math
import matplotlib.pyplot as plt 
import scipy.special as sp
plt.rcParams.update({'font.size': 18})
#
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
mu_3=5
print("mu_3=",mu_3)

#Amplitude constant, provided by res paper
#https://polymerdatabase.com/polymer%20physics/SAW.html
# A_2=1.1771
A_2=1



#Amplitude constant, provided by res paper
#https://polymerdatabase.com/polymer%20physics/SAW.html
# A_3=1.205
A_3=1
print("A_3=",A_3)

# number of short-mers in a large-mer, chain length "r"
# 6200/27=229, but here we want everything in terms of mers
r_Long=6200 



# number of single-mers in a arbitrary mer, chain length "r"
r_Arb=100



# number of single-mers in a arbitrary mer, chain length "r"
r_Short=27
#
#
# In[] 
# Grabbed 143750 from other project, "a"
#this value is single mer sites (not short mer)
#
n_surf=142600
print(n_surf)
# In[] Some new things, kind of like a main.



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
plt.figure()   
x_Long = np.array(r)
y_Long = reslog_Long
x_Line = r
y_Line = np.zeros(r.size)
# plt.xticks(range(min(x_Long), max(x_Long)+1))
# # Multiple plots
plt.xlabel('${n_a}$')
plt.ylabel('$Log(f(n_a))$')
plt.grid(b=None, which='major', axis='both')
# plt.title('Looping Behaviour Long') 
plt.plot(x_Long, y_Long,linewidth=2) 
#threshhold
plt.plot(x_Line, y_Line,'--r',linewidth=2)
#