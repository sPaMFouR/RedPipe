##########################################
## EPM Distance estimation of SN 2016gfy ##
##########################################

import numpy as np
import matplotlib.pyplot as plt
import subprocess
import scipy.optimize as optimize
from scipy.optimize import curve_fit
##########################################
## Parameters to be provided by the user #
## E(B-V)                                #
##########################################
ebv = 0.21
# Boltzmann constant in cm^2.g.s^-2.K^-1
k = 1.3807e-16
# Planck's constant in erg.s 
h = 6.6261e-27
# Speed of light in cm.s^-1
c = 2.99792458e10   
# Extinction calculated from E(B-V)
A_B = 4.05 * ebv
A_V = 3.1 * ebv
A_I = 0.48 * A_V

# 10**32 is multiplied to convert Angstrom (lambda in denominator) to cm
c1 = 2 * h * (c ** 2) * 1e32 
# 10**8 multiplied to convert cm (in c) to Angstrom
c2 = (h * c * 1e8) / k  

nTrial = 100

S = raw_input("Choose the filter set - 1. {BVI} 2. {BV} 3. {VI}: ")

# This code currently works for BVI only. So choose option 1.
if int(S) == 1:
   phase,mag_B, B_err, mag_V, V_err, mag_I, I_err= np.loadtxt('BVI_2016gfy.dat', unpack=True, usecols=(0,1,2,3,4,5,6))
elif int(S) == 2:
   phase,mag_B, B_err, mag_V, V_err = np.loadtxt('BVI_2016gfy.dat', unpack=True, usecols=(0,1,2,3,4))
elif int(S) == 3:
   phase,mag_V, V_err, mag_I, I_err = np.loadtxt('BVI_2016gfy.dat', unpack=True, usecols=(0,1,2,3,4))

#these are taken from Hamuy 2001 ApJ,615,642 Table 13
c0B=-45.144
c1B=7.159
c2B=-4.301
c3B=2.639
c4B=-0.811
c5B=0.098
c0V=-44.766
c1V=6.793
c2V=-4.523
c3V=2.695
c4V=-0.809
c5V=0.096
c0I=-44.345
c1I=6.347
c2I=-4.732
c3I=2.739
c4I=-0.811
c5I=0.096

tmpT = []
tmp_THETA = []
temp = []
zt = []
for j in range(len(mag_V)):
 a = np.array([]).reshape(0,2) 
 for iTrial in range(nTrial):
    mag_BT = mag_B + np.random.normal(scale=B_err,size=np.size(mag_B))
    mag_VT = mag_V + np.random.normal(scale=V_err,size=np.size(mag_V))
    mag_IT = mag_I + np.random.normal(scale=I_err,size=np.size(mag_I))
    def f(g):
        yb = c0B + c1B*(10/(g[0])) + c2B*(10/(g[0]))**2 + c3B*(10/(g[0]))**3 + c4B*(10/(g[0]))**4 + c5B*(10/(g[0]))**5
        yv = c0V + c1V*(10/(g[0])) + c2V*(10/(g[0]))**2 + c3V*(10/(g[0]))**3 + c4V*(10/(g[0]))**4 + c5V*(10/(g[0]))**5
        yi = c0I + c1I*(10/(g[0])) + c2I*(10/(g[0]))**2 + c3I*(10/(g[0]))**3 + c4I*(10/(g[0]))**4 + c5I*(10/(g[0]))**5
        print yb, yv, yi
# temp in kK, theta in 10**-11
#Dessert 2005 439,671
        if int(S) == 1:
           return ((mag_BT[j] - A_B + 5*np.log10(g[1]) -55 - yb)**2) + ((mag_VT[j] - A_V + 5*np.log10(g[1]) - 55 - yv)**2) + ((mag_IT[j] - A_I + 5*np.log10(g[1]) - 55 - yi)**2)
        elif int(S) == 2:
           return ((mag_BT[j] - A_B + 5*np.log10(g[1]) -55 - yb)**2)/ (B_err[j]**2) + ((mag_VT[j] - A_V + 5*np.log10(g[1]) - 55 - yv)**2)/ (V_err[j]**2)
        elif int(S) == 3:
           return ((mag_VT[j] - A_V + 5*np.log10(g[1]) -55 - yv)**2)/ (V_err[j]**2) + ((mag_IT[j] - A_I + 5*np.log10(g[1]) - 55 - yi)**2)/ (I_err[j]**2)

    result  = optimize.minimize(f, [6,1], method='L-BFGS-B')
    param_values = result.x
    print param_values
    a = np.vstack([a, param_values])
 np.savetxt('temp_zeta_theta_tmp.dat',a)
 T,THETA = np.loadtxt('temp_zeta_theta_tmp.dat',unpack=True, usecols=(0,1))

 temp.append(np.median(T))
 zt.append(np.median(THETA))
 tmpT.append(np.std(T))
 tmp_THETA.append(np.std(THETA))
 
#a = np.array([]).reshape(0,2)
#for j in range(len(mag_V)):
#    def f(g):
#        yb = c0B + c1B*(10/(g[0])) + c2B*(10/(g[0]))**2 + c3B*(10/(g[0]))**3 + c4B*(10/(g[0]))**4 + c5B*(10/(g[0]))**5
#        yv = c0V + c1V*(10/(g[0])) + c2V*(10/(g[0]))**2 + c3V*(10/(g[0]))**3 + c4V*(10/(g[0]))**4 + c5V*(10/(g[0]))**5
#        yi = c0I + c1I*(10/(g[0])) + c2I*(10/(g[0]))**2 + c3I*(10/(g[0]))**3 + c4I*(10/(g[0]))**4 + c5I*(10/(g[0]))**5
## temp in kK, theta in 10**-11
##Dessert 2005 439,671
#        if int(S) == 1:
#           return ((mag_B[j] - A_B + 5*np.log10(g[1]) -55 - yb)**2)/ (B_err[j]**2) + ((mag_V[j] - A_V + 5*np.log10(g[1]) - 55 - yv)**2)/ (V_err[j]**2) + ((mag_I[j] - A_I + 5*np.log10(g[1]) - 55 - yi)**2)/ (I_err[j]**2)
#        elif int(S) == 2:
#           return ((mag_B[j] - A_B + 5*np.log10(g[1]) -55 - yb)**2)/ (B_err[j]**2) + ((mag_V[j] - A_V + 5*np.log10(g[1]) - 55 - yv)**2)/ (V_err[j]**2)
#        elif int(S) == 3:
#           return ((mag_V[j] - A_V + 5*np.log10(g[1]) -55 - yv)**2)/ (V_err[j]**2) + ((mag_I[j] - A_I + 5*np.log10(g[1]) - 55 - yi)**2)/ (I_err[j]**2)
#    result  = optimize.minimize(f, [6,1], method='L-BFGS-B')
#    param_values = result.x
#    a = np.vstack([a, param_values])

#np.savetxt('temp_zeta_theta.dat',a)
#temp, zt = np.loadtxt('temp_zeta_theta.dat', unpack=True, usecols=(0,1))
np.savetxt('temp_zeta_theta_err.dat',np.c_[temp,tmpT,zt,tmp_THETA,phase])
temp,temp_err, zt, zt_err = np.loadtxt('temp_zeta_theta_err.dat', unpack=True, usecols=(0,1,2,3))
plt.plot(phase,temp)
plt.errorbar(phase,temp,temp_err,linestyle="None")
plt.show()

# Hamuy 2001 table 14
# {BVI}
a0_BVI_H01 = 0.7336
a1_BVI_H01 = -0.6942
a2_BVI_H01 = 0.3740
# {BV}
a0_BV_H01 = 0.7557
a1_BV_H01 = -0.8997
a2_BV_H01 = 0.5199
# {VI}
a0_BV_H01 = 0.7013
a1_BV_H01 = -0.5304
a2_BV_H01 = 0.2646
# Dessart & Hillier, A&A, 439, 671 (2005)
# {BVI}
a0_BVI_D05 = 0.63241
a1_BVI_D05 = -0.38375
a2_BVI_D05 = 0.28425
# {BV}
a0_BV_D05 = 0.47188
a1_BV_D05 = -0.25399
a2_BV_D05 = 0.32630
# {VI}
a0_VI_D05 = 0.81662
a1_VI_D05 = -0.62896
a2_VI_D05 = 0.33852

import uncertainties as unc  
import uncertainties.unumpy as unumpy  
temp = unumpy.uarray(( temp, temp_err ))
zt =   unumpy.uarray(( zt, zt_err ))
if int(S) == 1:
#these are now unumpy
   zeta_BVI_H01 = a0_BVI_H01 + a1_BVI_H01*(10/temp) + a2_BVI_H01*(10/temp)**2
   zeta_BVI_D05 = a0_BVI_D05 + a1_BVI_D05*(10/temp) + a2_BVI_D05*(10/temp)**2
   theta_BVI_H01 = zt / zeta_BVI_H01
   theta_BVI_D05  = zt / zeta_BVI_D05
   theta_BVI_D05err=unumpy.std_devs(theta_BVI_D05) 
   theta_err = np.sqrt((unumpy.nominal_values(zt))**2 + (unumpy.std_devs(zeta_BVI_D05))**2)
   np.savetxt('temp_theta.dat',np.c_[(phase),unumpy.nominal_values(temp),unumpy.std_devs(temp),unumpy.nominal_values(zeta_BVI_D05),unumpy.std_devs(zeta_BVI_D05),((unumpy.nominal_values(theta_BVI_D05)*10)/3.24),((theta_BVI_D05err*10)/3.24)])
# {BV}
elif int(S) == 2:
   zeta_BV_H01 = a0_BV_H01 + a1_BV_H01*(10/temp) + a2_BV_H01*(10/temp)**2
   zeta_BV_D05 = a0_BV_D05 + a1_BV_D05*(10/temp) + a2_BV_D05*(10/temp)**2
   theta_BV_H01  = zt / zeta_BV_H01
   theta_BV_D05  = zt / zeta_BV_D05
# {VI}
elif int(S) == 3:
   zeta_VI_H01 = a0_BVI_H01 + a1_BVI_H01*(10/temp) + a2_BVI_H01*(10/temp)**2
   zeta_VI_D05 = a0_VI_D05 + a1_VI_D05*(10/temp) + a2_VI_D05*(10/temp)**2
   theta_VI_H01  = zt / zeta_VI_H01
   theta_VI_D05  = zt / zeta_VI_D05

y,z,zerr = np.loadtxt("VelErr_2016gfy.dat", unpack = True, usecols = [0,1,2])
y1,z1 = np.loadtxt("Vel_2016gfy.dat", unpack = True, usecols = [0,1])
plt.plot(y,z)
plt.errorbar(y,z,zerr,linestyle="None")
plt.scatter(y1,z1,c='r')
plt.show()

tmp_H01 = np.array([]).reshape(0)
tmp_D05 = np.array([]).reshape(0)
tmp_H01_err = np.array([]).reshape(0)
tmp_D05_err = np.array([]).reshape(0)
tmp_z = np.array([]).reshape(0)
tmp_zerr = np.array([]).reshape(0)
tmp_y = np.array([]).reshape(0)
for i in range(len(y)):
    for j in range(len(phase)):
        if abs(y[i]-phase[j]) < 0.0131 and int(S) == 1:
           tmp_H01 = np.hstack((tmp_H01,theta_BVI_H01[j]))
           tmp_D05 = np.hstack((tmp_D05,theta_BVI_D05[j]))
           #tmp_D05_err = np.hstack((tmp_D05_err,theta_BVI_D05err[j]))
           #tmp_D05_err = np.hstack((tmp_D05_err,theta_err[j]))
           tmp_y = np.hstack((tmp_y,y[i]))
           tmp_z = np.hstack((tmp_z,z[i]))
           tmp_zerr = np.hstack((tmp_zerr,zerr[i]))
        elif abs(y[i]-phase[j]) < 0.01 and int(S) == 2:
           tmp_H01 = np.hstack((tmp_H01,theta_BV_H01[j]))
           tmp_D05 = np.hstack((tmp_D05,theta_BV_D05[j]))
           tmp_y = np.hstack((tmp_y,y[i]))
           tmp_z = np.hstack((tmp_z,z[i]))
        elif abs(y[i]-phase[j]) < 0.01 and int(S) == 3:
           tmp_H01 = np.hstack((tmp_H01,theta_VI_H01[j]))
           tmp_D05 = np.hstack((tmp_D05,theta_VI_D05[j]))
           tmp_y = np.hstack((tmp_y,y[i]))
           tmp_z = np.hstack((tmp_z,z[i]))
# unumpy arrays
x_H01 = tmp_H01
x_D05 = tmp_D05
print x_D05
y = tmp_y
z = tmp_z
zerr = tmp_zerr
np.savetxt('Spline_2016gfy.dat',np.c_[y,z,zerr])
z = unumpy.uarray((z,zerr))

y1,z1 = np.loadtxt("Vel_2016gfy.dat", unpack = True, usecols = [0,1])
ynew,znew,znewerr = np.loadtxt("Spline_2016gfy.dat", unpack = True, usecols = [0,1,2])
plt.plot(ynew,znew)
plt.errorbar(ynew,znew,znewerr,linestyle="None")
plt.scatter(y1,z1,c='r')
plt.show()
def line(x, a, b):
    return a * x + b
xpar_H01 = x_H01/z
xpar_D05 = x_D05/z
#print xpar
ypar = y 
popt, pcov = curve_fit(line,unumpy.nominal_values(xpar_H01),ypar)
yfit = line(xpar_H01,*popt)
print "a =", popt[0], "+/-", pcov[0,0]**0.5
print "b =", popt[1], "+/-", pcov[1,1]**0.5
#1m = 3.24*10**-23 Mpc 
print "distance(in Mpc) is (method H01)"
print popt[0]*24*3600*10**14*3.24*10**(-23)
print "error"
print pcov[0,0]**0.5*24*3600*10**14*3.24*10**(-23)
#print "e" 
#plt.plot(xpar_H01,y, ls = "none", marker = 'o', color ="blue")

#plt.plot(xpar_H01,yfit)
#plt.ylim(plt.ylim()[::-1])
#plt.show()
popt, pcov = curve_fit(line,unumpy.nominal_values(xpar_D05),ypar,sigma=unumpy.std_devs(xpar_D05))
yfit = line(unumpy.nominal_values(xpar_D05),*popt)
print "a =", popt[0], "+/-", pcov[0,0]**0.5
print "b =", popt[1], "+/-", pcov[1,1]**0.5
#1m = 3.24*10**-23 Mpc 
print "distance(in Mpc) is (method D05)"
print popt[0]*24*3600*10**14*3.24*10**(-23)
print "error"
print pcov[0,0]**0.5*24*3600*10**14*3.24*10**(-23)
#print "e" 
plt.plot(unumpy.nominal_values(xpar_D05)*10**3,ypar, ls = "none", marker = 'o', color ="blue")
plt.errorbar(unumpy.nominal_values(xpar_D05)*10**3, ypar, xerr = unumpy.std_devs(xpar_D05)*10**3 , ls = "none", marker = 'o', color ="blue")
plt.plot(unumpy.nominal_values(xpar_D05)*10**3,yfit)
plt.xlabel(r'$\theta$/v$_{ph}$ (10$^{-17}$ rad m$^{-1}$ s)',fontsize=15)
plt.ylabel('Days since discovery',fontsize=15)
plt.savefig('EPM_2016gfy_D05_BVI.pdf')
#plt.ylim(plt.ylim()[::-1])
#plt.ylim((4,41))
#plt.xlim((0.05,0.36))
plt.show()
