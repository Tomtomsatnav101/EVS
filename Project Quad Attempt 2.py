# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 10:05:36 2021

@author: tomsw
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import time
from scipy.interpolate import CubicSpline

#Constants

h = 0.7
#R = 8/h
A_s = 1                              #0.81/92276.53381005568 #0.81/(114790)
n_s = 0.96
omega_c = 0.23
omega_b = 0.046
omega_r = 8.6e-5
omega_m = omega_c + omega_b
D_0 = 0.76
c = 3e5
H_0 = 100*h
G = 6.67*10**(-11)
rho_crit = (3*H_0)/(8*np.pi*G)
rho_m = omega_m * rho_crit
curly_R = 1
sqrt_a = 0.9
delta_c = sqrt_a * 1.686

div = 100
omega_l = 1 - omega_r - omega_m
f_sky = 1
a = 0.707


alpha = 1-0.328*np.log(431*omega_m*(h**2))*(omega_b/omega_m) + 0.38*np.log(22.3*omega_m*(h**2))*(omega_b/omega_m)**2
s = (44.5*np.log(9.83/(omega_m*(h**2))))/((1+10*((omega_b*(h**2))**0.75))**0.5)

def W(k):
    if k <=0.05:
        return 1.0 - (k**2)/10.0 + (k**4)/280.0 - (k**6)/15120.0 + (k**8)/1330560.0
    else:
        return 3*(np.sin(k)/(k**3) - np.cos(k)/(k**2))

def k_EH(k):
    return ((c*k*(omega_r**0.5))/(H_0*omega_m)) / (alpha + (1-alpha)/(1+(0.43*k*s)**4))

def T(k):
    if k<=0.05:
        return 1.0-0.5072815*(k**2)+0.9178786265*(k**4)-1.962195729*(k**6)+4.283146237*(k**8)
    else:
        return ((np.log(1+(0.124*k)**2))/((0.124*k)**2)) * ((1+((1.257*k)**2) + ((0.4452*k)**4) + ((0.2197*k)**6))/(1+((1.606*k)**2) + ((0.8568*k)**4) + ((0.3927*k)**6)))**0.5

def A(k,z):
    return (2/(3*omega_m)) * (k*c/H_0)**2 * T(k_EH(k)) * dee(z)

def P(k,z):
    return (A(k,z))**2 * (c*k/H_0)**(n_s-1)


def f(x,R,z):
    return ((W(R*x))**2)*P(x,z)*(x**(-1))

def g(x,R,z):
    return ((W(R/x))**2)*P(1/x, z)*(x**(-1))

def quad_8():
    R = 8/h
    I1 = integrate.quad(f,1,np.inf,args=(R,0))
    I2 = integrate.quad(g,1,np.inf, args=(R,0))
    I3 = I1[0]+I2[0]
    I = A_s*np.sqrt(4*np.pi*I3)
    return I

def sigma_R(R, z):
    I1 = integrate.quad(f,1,np.inf,args=(R,z))
    I2 = integrate.quad(g,1,np.inf, args=(R,z))
    I3 = I1[0]+I2[0]
    I = A_s*np.sqrt(4*np.pi*I3)
    return I

def calc_radius(M):
    return (M/(1.16*(10**12)*omega_m*h*h))**(1/3)

def calc_mass(R):
    return 1.16*(10**12)*omega_m*((R*h)**3)/h

def Press_Schechter(v):
    return np.sqrt(2/np.pi) * v * np.exp(-0.5*(v**2))

def Sheth_Tormen(v):
    return 0.322 * np.sqrt(2*a/np.pi) * np.exp((-a*v**2)/2)*(1+(a*(v**3))**(-0.3))

def beta(z):
    return 0.589*(1+z)**0.2
    
def phi(z):
    return -0.729*(1+z)**(-0.08)

def eta(z):
    return -0.243*(1+z)**0.27

def gamma(z):
    return 0.864*(1+z)**(-0.01)

def Tinker(v,z):
    return 0.368*(1+(beta(z)*v)**(-2*phi(z)))*(v**(2*eta(z)+1))*(np.exp((-gamma(z)*(v**2))/2))

def calc_dsigma_dU(U,z):    # took h inside brackets
    delta=10**(-3)          #      <-------dsigma_R/dR------->                                              <-------------------dR/dM term------------------------->           <--dM/du-->
    return ((sigma_R(delta+calc_radius((10**U)/h),z)-sigma_R(-delta+calc_radius((10**U)/h),z))/(2*delta)) * ((1/3)*((10**U)**(-2/3)) * (1.16*(10**12)*omega_m*h*h)**(-1/3)) *(10**U)*np.log(10)

def dn_dU(U,z):# 
    #temp = (-1/sigma_R(calc_radius(10**U),z))
    #print("test="+str(temp))
    #                                                                                                                                                      <-----dM/du----->
    return curly_R * Press_Schechter(delta_c/sigma_R(calc_radius((10**U)/h),z)) * (rho_m/(10**U)) * (-1/sigma_R(calc_radius((10**U)/h),z)) * calc_dsigma_dU(U,z) #* (10**U)*np.log(10)
#                  took h inside bracket                                                                                                                         -Dont need as returns dsigma/du 
def mean_num_den(U, M, z):
    I, err = integrate.quad(dn_dU, U, M, args=(z))  #M = upper bound
    return I


def number(z, U):
    dVdz = dV_dz(z)
    mnd = mean_num_den(U,np.inf, z)
    return dVdz*mnd#*(10**U)  ########THINK THIS MIGHT FIX IT

#sets value of A_s
def normalise():
        global A_s
        A_s= 0.801/quad_8()

def dee(z):
    zcube = (1+z)**3
    omega_z = (omega_m*zcube)/(omega_l + omega_r*((1+z)**4) + zcube*omega_m)    
    Oxz = omega_l/(omega_l + omega_r*((1+z)**4) + zcube*omega_m)
    gz = 2.5*omega_z/(omega_z**(4/7)-Oxz+(1+0.5*omega_z)*(1+Oxz/70))    
    Dplus = gz/(1+z)
    return Dplus

def H(z):
    return H_0*((omega_m*(1+z)**3 + omega_l)**0.5)/c

def invH(z):
    return c/(H_0*((omega_m*(1+z)**3 + omega_l)**0.5))   ##Factor of c maybe

def dV_dz(z):
    I,err = integrate.quad(invH,0,z)
    return (f_sky*4*np.pi/H(z))*(I**2)

def plot_n():
    zz = np.linspace(0,1.1,12)
    M = np.linspace(12, 16, 5)
    N = np.zeros((5, 11))

    for i in range(5):
        print("Done "+str(i))
        for j in range(11):
            print("Donee "+str(j))
            N[i][j] = integrate.quad(number,zz[j], zz[j+1], args=(M[i]))[0]

    N0 = N[0][:]
    N1 = N[1][:]
    N2 = N[2][:]
    N3 = N[3][:]
    N4 = N[4][:]
    
    plt.figure(0)
    plt.semilogy(zz[0:11], N0, '+' , label="M="+str(M[0]))
    plt.semilogy(zz[0:11], N1, '+' ,label="M="+str(M[1]))
    plt.semilogy(zz[0:11], N2, '+' ,label="M="+str(M[2]))
    plt.semilogy(zz[0:11], N3, '+' ,label="M="+str(M[3]))
    plt.semilogy(zz[0:11], N4, '+' ,label="M="+str(M[4]))
    plt.xlabel("Redshift")
    plt.ylabel("n(>M, z)")
    plt.title('Number of Clusters larger than M for given redshift')
    plt.legend(loc="lower left",fontsize = 'small')
    plt.xlim([-0.05, 1.05])

def P0_dU(z,U):
    dndu = dn_dU(U,z)
    #print("dndu = %s" %dndu)
    n = number(z, U) # = n(>m,z)*dV/dz
    #print("n = %s" %n)
    #print("exp(-n) = %s" %np.exp(-n))
    #dVdz = dV_dz(z)
    #print("dVdz = %s" %dVdz)
    return -np.exp(-n)*dndu

def m_0(z):
    U = np.linspace(15.069,15.071,21) #15.324, 15.326
    M = np.zeros(21)
    z_mean=z+0.05
    V = integrate.quad(dV_dz,z,z+0.1)[0]
    for i in range(21):
        M[i] = -np.log(V) + np.log(dn_dU(U[i], z_mean)) - mean_num_den(U[i], np.inf, z_mean)*V

    plt.figure(1)
    plt.plot(U[0:21], M,'+')
    plt.title("p_G(m)")
    plt.xlabel("log_10(M/Solar Masses)")
    plt.ylabel("log(dP/dm)")
    
def dP_0_dm():
    #%%
    #Use this function to manually search for m_0
    z = 1.1
    steps=21
    U = np.linspace(14.988, 14.99, steps)
    N = np.zeros(steps)

    V = integrate.quad(dV_dz,z,z+0.1)[0]
    print(V)

    for i in range(steps):
        #print("Done "+str(i))
        mnd = mean_num_den(U[i],np.inf, z+0.05)
        #print("mnd = ", mnd)
        #print(-np.exp(-mnd*V))
        #print(dn_dU(U[i], z))
        N[i] = np.exp(-mnd*V)*V*dn_dU(U[i], z+0.05)  #Should have a - sign
        #N[j][0] = integrate.quad(P0_dU,0,0.1, args=(U[j]))[0]


    #N0 = N[:][0]
    #N1 = N[:][1]
    #N2 = N[:][2]
    #N3 = N[:][3]
    #N4 = N[:][4]
    
    plt.figure(505)
    plt.semilogy(U[0:steps], N, '+' , label="z="+str(z))
    # plt.plot(U[0:5], N1, '+' ,label="z="+str(z[1]))
    # plt.plot(U[0:5], N2, '+' ,label="z="+str(z[2]))
    # plt.plot(U[0:5], N3, '+' ,label="z="+str(z[3]))
    # plt.plot(U[0:5], N4, '+' ,label="z="+str(z[4]))
    plt.xlabel("M")
    plt.ylabel("dP0/dm")
    plt.title('?')
    plt.legend(loc="lower left",fontsize = 'small')
    #plt.xlim([-0.05, 1.05])

#%%
def EVS_para(u0,z):
    V = integrate.quad(dV_dz,z,z+0.1)[0]
    z_mean = z+0.05

    gamma = mean_num_den(u0, np.inf, z_mean) *V -1
    beta = ((1+gamma)**(1+gamma)) / (dn_dU(u0, z_mean) * V * u0 *np.log(10))
    alpha = u0 - beta/gamma *(((1+gamma)**(-gamma))-1)
    return gamma, alpha, beta
    


def inv_Weibull(x, gamma, alpha,beta):
    inv_F = alpha + (beta/gamma)*(((-np.log(x))**(-gamma))-1)
    return inv_F
    
    
###############Program########################################################################################
#%%
t = time.time()
normalise()                     #Calibrates A_s and outputs known value of sigma_8
print("sigma_8 = %s" %sigma_R(calc_radius(calc_mass(8/h)),0))
print("D(0) = %s" %dee(0))              #Confirms the given value of D(0)

U = np.linspace(14,16,div)

z=0

# sig_R = np.zeros(div)
# for i in range(div):
#     sig_R[i] = sigma_R(calc_radius(10**U[i]),z)   

# plt.figure(0) 
# plt.plot(U, sig_R, label="z="+str(z))            #Calculates and plots sigma_R against M
# plt.xlabel("M/ Solar Masses")
# plt.ylabel("Sigma_R")
# plt.legend(loc="upper right")

# dsig_R_dU = np.zeros(div)
# for i in range(div):
#     dsig_R_dU[i] = calc_dsigma_dU(U[i], z)

# plt.figure(1)               #Calculates and plots dsigma_R/dM against M
# plt.plot(U,dsig_R_dM, label="z="+str(z))
# plt.xlabel("M/ Solar Masses")
# plt.ylabel("dSigma_R/dM")
# plt.legend(loc="lower right")

# dndU = np.zeros(div)      #JUST COMMENTED THIS OUT
# for i in range(div):
#     dndU[i] = dn_dU(U[i], z)

# plt.figure(0)
# plt.plot(U, np.log10(dndU), label="z="+str(z))               #Calculates and plots dn/dM against M
# plt.xlabel("M/Solar Masses")
# plt.ylabel("dn/dM")
# plt.legend(loc="lower left")      #JUST COMMENTED THIS OUT

# num_den = np.zeros(div)
# for i in range(div):
#     num_den[i] = mean_num_den(U[i],np.inf, z)
    
# plt.figure(3)
# plt.plot(U,np.log10(num_den),'+', label="z="+str(z))     #Calculates and plots n(>m) against M
# plt.xlabel("log_10(M/Solar Masses)")
# plt.ylabel("log_10(Number density)")
# plt.legend(loc="lower left")

# num = np.zeros(div)
# for i in range(div):
#     #print("i = %s" %i)
#     num[i] = number(U[i],z)

# plt.figure(4)
# plt.plot(U, np.log10(num), '+', label="z="+str(z))
# plt.xlabel("log_10(M/Solar Masses)")
# plt.ylabel("log_10(n(>M))")
# plt.legend(loc="lower left")

# U1 = 13     ##JUSTCOMMENTED THIS OUT
# U2 = 14 
# U3 = 15 

# z=np.linspace(0,5,div)

# dndU1 = np.zeros(div)
# for i in range(div):
#     dndU1[i] = dn_dU(U1, z[i])
    
# print("U1 done")
# dndU2 = np.zeros(div)
# for i in range(div):
#     dndU2[i] = dn_dU(U2, z[i])    

# print("U2 done")
# dndU3 = np.zeros(div)
# for i in range(div):
#     dndU3[i] = dn_dU(U3, z[i])

# print("U3 done")
# plt.figure(1)
# plt.semilogy(z, dndU1, label="M="+str(U1))               #Calculates and plots dn/dM against M
# plt.semilogy(z, dndU2, label="M="+str(U2))
# plt.semilogy(z, dndU3, label="M="+str(U3))
# plt.xlabel("Z")
# plt.ylabel("dn/dU")
# plt.legend(loc="upper right")
# plt.xlim([0,5])
#plt.ylim([1e-10, 1e-3])  #JUST COMMENTED THIS OUT


# z = np.linspace(0, 10, 101)    #JUST COMMENTED THIS
# M = np.linspace(12, 16, 5)

# num_den = np.zeros((4,101))

# for i in range(3):
#     print("Done "+str(i))
#     for j in range(101):
#         num_den[i][j] = mean_num_den(M[i], M[i+1], z[j])

# n1 = num_den[0][:]
# n2 = num_den[1][:]
# n3 = num_den[2][:]
# #n4 = num_den[3][:]

# plt.figure(0)
# plt.semilogy(z, n1, label="M="+str(M[0])+" to "+str(M[1]))
# plt.semilogy(z, n2, label="M="+str(M[1])+" to "+str(M[2]))
# plt.semilogy(z, n3, label="M="+str(M[2])+" to "+str(M[3]))
# #plt.semilogy(z, n4, label="M="+str(M[3])+" to "+str(M[4]))
# plt.xlabel("Z")
# plt.ylabel("n(M_1, M_2, Z)")
# plt.legend(loc="lower left")

# plt.figure(1)
# plt.semilogy(z[0:50], n1[0:50], label="M="+str(M[0])+" to "+str(M[1]))
# plt.semilogy(z[0:50], n2[0:50], label="M="+str(M[1])+" to "+str(M[2]))
# plt.semilogy(z[0:50], n3[0:50], label="M="+str(M[2])+" to "+str(M[3]))
# #plt.semilogy(z[0:50], n4[0:50], label="M="+str(M[3])+" to "+str(M[4]))
# plt.xlabel("Z")
# plt.ylabel("n(M_1, M_2, Z)")
#plt.legend(loc="lower left")     #JUST COMMENTED THIS END


#Sheth-Torman
u0st = [15.325, 15.45, 15.465, 15.448, 15.415, 15.374, 15.328, 15.279, 15.227, 15.175, 15.123, 15.07]
#Press-Schecter
u0ps = [15.342, 15.431, 15.431, 15.404, 15.364, 15.317, 15.256, 15.213, 15.157, 15.102, 15.045, 14.989]
z_mean = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15]

plt.figure(1)
plt.plot(z_mean, u0st,label="Sheth-Torman")
plt.plot(z_mean,u0ps, label="Press-Schecter")
plt.xlabel("z")
plt.ylabel("m_0")
plt.title("Mass at which p_G peaks for each redshift bin")
plt.legend(loc="lower left")

gamma1 = np.zeros(12)
alpha1 = np.zeros(12)
beta1 = np.zeros(12)

for i in range(12):
    temp = EVS_para(u0ps[i], z_mean[i]-0.05)    
    gamma1[i] = temp[0]
    alpha1[i] = temp[1]
    beta1[i] = temp[2]

plt.figure(2)
plt.plot(u0ps-alpha1)
plt.title("Difference between m0 and alpha")

#%%
x = np.linspace(14.988,14.99,21)
y=np.zeros(21)
f1 = np.zeros(21)
#F1 = np.zeros(21)
z_bin = 11

for i in range(21):
    y[i] = (x[i]-alpha1[z_bin])/(beta1[z_bin])
    f1[i] = -(1+gamma1[z_bin]*y[i])**(-1/gamma1[z_bin]) - np.log(beta1[z_bin]) - (1+1/gamma1[z_bin])*np.log(1+gamma1[z_bin]*y[i])
    #F1[i] = np.exp(-(1+gamma1[z_bin]*y[i])**(-1/gamma1[z_bin]))
    #print(f1[i])
    
plt.figure(3)
plt.plot(x,f1)
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Weibull Distribution")
#%%
#Need to find m0 for 0<z<1

U = np.linspace(15.549,15.551,21) #15.324, 15.326
M = np.zeros(21)
z_mean=0.5
V = integrate.quad(dV_dz,0,1)[0]
for i in range(21):
    M[i] = -np.log(V) + np.log(dn_dU(U[i], z_mean)) - mean_num_den(U[i], np.inf, z_mean)*V

plt.figure(4)
plt.plot(U[0:21], M,'+')
plt.title("p_G(m)")
plt.xlabel("log_10(M/Solar Masses)")
plt.ylabel("log(dP/dm)")


u0 = 15.55
V = integrate.quad(dV_dz,0,1)[0]
z_mean = 0.5

uni_gamma = mean_num_den(u0, np.inf, z_mean) *V -1
uni_beta = ((1+uni_gamma)**(1+uni_gamma)) / (dn_dU(u0, z_mean) * V * u0 *np.log(10))
uni_alpha = u0 - uni_beta/uni_gamma *(((1+uni_gamma)**(-uni_gamma))-1)

print(uni_gamma)
print(uni_alpha)
print(uni_beta)

x = np.linspace(15.549,15.551,41)
y=np.zeros(41)
f1 = np.zeros(41)
F1 = np.zeros(41)

for i in range(41):
    y[i] = (x[i]-uni_alpha)/(uni_beta)
    f1[i] = -(1+uni_gamma*y[i])**(-1/uni_gamma) - np.log(uni_beta) - (1+1/uni_gamma)*np.log(1+uni_gamma*y[i])
    #F1[i] = np.exp(-(1+gamma1[z_bin]*y[i])**(-1/gamma1[z_bin]))
    print(f1[i])
    
plt.figure(5)
plt.plot(x,f1)
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Weibull Distribution")
#%%
#Sheth-Torman
Mmaxst = [15.325, 15.45, 15.465, 15.448, 15.415, 15.374, 15.328, 15.279, 15.227, 15.175, 15.123, 15.07]
#Press-Schecter
Mmaxps= [15.342, 15.431, 15.431, 15.404, 15.364, 15.317, 15.256, 15.213, 15.157, 15.102, 15.045, 14.989]
z_mean = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15]

plt.figure(500)
cs_mst = CubicSpline(z_mean, Mmaxst)
cs_mps = CubicSpline(z_mean, Mmaxps)
xs = np.arange(0.05, 1.15, 0.01)

plt.plot(xs, cs_mst(xs), label="Sheth-Torman")
plt.plot(xs, cs_mps(xs), label="Press-Schechter")
plt.ylim(15,15.5)
plt.xlabel("Redshift Bin")
plt.ylabel("M_max")
plt.title("Most likely maximum mass in each redshift bin")
plt.legend(loc="lower left")


Umaxst = np.zeros(12)
Umaxps = np.zeros(12)
for i in range(12):
    Umaxst[i] = 10**Mmaxst[i]
    Umaxps[i] = 10**Mmaxps[i]

plt.figure(501)
cs_ust = CubicSpline(z_mean, Umaxst)
cs_ups = CubicSpline(z_mean, Umaxps)

plt.plot(xs, cs_ust(xs) ,label="Sheth-Torman")
plt.plot(xs, cs_ups(xs) ,label="Press-Schechter")
plt.ylim(10**15,10**15.5)
plt.xlabel("Redshift Bin")
plt.ylabel("M_max (10^15 Solar Masses)")
plt.title("Most likely maximum mass in each redshift bin")
plt.legend(loc="lower left")
#%%

F_low = np.zeros(z_bin+1)
F_high = np.zeros(z_bin+1)
for i in range(z_bin+1):
    F_low[i] = inv_Weibull(0.05,gamma1[i], alpha1[i], beta1[i])
    F_high[i] = inv_Weibull(0.999,gamma1[i], alpha1[i], beta1[i])

cs_Fhigh_ps = CubicSpline(z_mean, F_high)
cs_Flow_ps = CubicSpline(z_mean, F_low)
xs = np.arange(0.05, 1.15, 0.01)

plt.figure(600)
plt.plot(xs,cs_Flow_ps(xs), label="Lower bound")
plt.plot(xs,cs_Fhigh_ps(xs), label="Upper bound")
plt.plot(xs, cs_mps(xs), label="Press-Schechter")
plt.title('Press-Schechter')
plt.xlabel("Redshift Bin")
plt.ylabel("M_max")
plt.legend(loc="lower left")

plt.figure(601)
plt.plot(z_mean,F_low, label="Lower bound")
plt.plot(z_mean,F_high, label="Upper bound")
plt.plot(z_mean, Mmaxps, label="Press-Schechter")
plt.title('Press-Schechter')
plt.xlabel("Redshift Bin")
plt.ylabel("log(M_max/Solar Masses)")
plt.legend(loc="upper right")

data = [0.869, 0.819, 1.968, 1.26, 1.083]
error = [0.063, 0.013, 0.048, 0.02, 0.024]

log_data = np.log10(data)+15
log_error = np.zeros([2,5])

for i in range(5):
    log_error[0,i] = -(np.log10(data[i]-error[i]) + 15 - log_data[i])
    log_error[1,i] = np.log10(data[i]+error[i]) + 15 - log_data[i]

plt.errorbar(z_mean[0:5], log_data, yerr=log_error, marker="*", linestyle="none")

#%%
#Planck data

lin_F_low = np.zeros(z_bin+1)
lin_F_high = np.zeros(z_bin+1)
for i in range(z_bin+1):
    lin_F_low[i] = 10**F_low[i]/(10**15)
    lin_F_high[i] = 10**F_high[i]/(10**15)

lin_cs_Fhigh_ps = CubicSpline(z_mean, lin_F_high)
lin_cs_Flow_ps = CubicSpline(z_mean, lin_F_low)

Umaxps = np.zeros(12)
for i in range(12):
    Umaxps[i] = 10**Mmaxps[i]/(10**15)
    

lin_cs_mps = CubicSpline(z_mean, Umaxps)
xs = np.arange(0.05, 1.15, 0.01)

plt.figure(602)

data_2 = [0.8771, 0.8859, 1.6116, 1.4693 ,1.2250, 1.1487, 1.0727, 0.9481, 1.0754, 0.6774]
u_error = [0.0186, 0.0323, 0.0297, 0.0392, 0.0525, 0.0535, 0.0630, 0.0669, 0.0478, 0.0487]
l_error =[0.0209, 0.0320, 0.0292, 0.0417, 0.0550, 0.0548, 0.0664, 0.0530, 0.0472, 0.0535]
error_2 = np.zeros([2,10])
for i in range(10):
    error_2[0,i] = l_error[i]
    error_2[1,i] = u_error[i]

plt.plot(xs,lin_cs_Flow_ps(xs), label="Lower bound")
plt.plot(xs,lin_cs_Fhigh_ps(xs), label="Upper bound")
plt.plot(xs, lin_cs_mps(xs), label="Press-Schechter")
plt.scatter(z_mean, Umaxps, marker="x")
plt.scatter(z_mean, lin_F_high, marker="x")
plt.scatter(z_mean, lin_F_low, marker="x")

plt.errorbar(z_mean[0:10], data_2, yerr=error_2, marker="_", linestyle="none", capsize=6)
plt.xlabel("Redshift bin")
plt.ylabel("Mass (10^15 solar masses)")
plt.title("Press-Schechter")
plt.legend(loc="upper right")

#Zoom on error bars
#plt.xlim([0,0.2])
#plt.ylim([14.8, 15])

#Zoom on lines
#plt.ylim([15.3, 15.5])
#plt.xlim([0, 0.4])

width = F_high-F_low
print(width)

delta = Mmaxps - F_low
print(delta)

#%%

elapsed = time.time()-t
print("%s seconds" %elapsed)