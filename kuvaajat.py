#%% 

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import constants

# global constants
R = constants.gas_constant
T_0 = 373.15 # boiling point in kelvins at regular air pressure
P_0 = 753 # ambient air pressure in mmHg
water_mol_mass = 18.01528 # grams per mol
mass_container = 34 # grams


# ----- IMPORTING DATA ----- #
mass, teho_in, teho_out = np.genfromtxt("data1uusi.txt", comments='#', unpack=True, missing_values='\n')
temp, psi = np.genfromtxt("data2uusi.txt", comments='#', missing_values='\n', unpack=True)

# mass difference in grams
mass = mass - mass_container
# conversion to kg
mass = mass*0.001 

# teho difference, in kWh
teho = teho_out - teho_in

# converting kWh to si-units 
# watt = Joule per second
# 1 hour = 3600 seconds
# 3600 * measured value = watt seconds from which we get joules
# I don't remember all the details but the following is in appropriate units lol
teho = teho*3600

# the actual pressure in mmHg
psi = P_0 - psi
# converting mmHg to pascal
# 1 mmHg = 101325/760 Pa
P_0 = P_0*(101325/760) # updating P_0 to 
psi = psi*(101325/760)
#print(psi)

# converting temps from celcius to kelvin
temp = temp + 273.15
#print(temp)

# ---- INITIAL GUESS ---- #
psiError = np.zeros_like(psi)+(0.5*(101325/760))
alkuarvaus = np.array( [1000.0, -10.0])

# ----- FUNCTION DEFINITIONS ----- #
# function for linear fitting
def linear(x, a, b):
    y = a*x + b
    return y

#--- exponential function used for fitting ---
# T_0 and P_0 are constants that have to be given when called
# e.g. room temperature and or air pressure
# a and b are fit parameters
# T and P are measured variables
# T being temperature
# T_0 is measured room temperature
# P_0 is measured room air pressure
def exponential(T, a, b):
    P = P_0*np.exp(a*(1/T_0 - 1/T))*(T/T_0)**b
    return P


# ---- LINEAR FITTING ---- #
popt = optimize.curve_fit(linear, xdata=mass, ydata=teho)
parameters = popt[0]
popterr = np.sqrt(np.diag(popt[1]))
print('optimal parameteres for the linear fit are ', parameters, '\n and the errors are ', popterr)

x1 = np.linspace(0, 0.001*50, 2000)

linearfit = linear(x1, parameters[0], parameters[1])


# ----- EXPONENTIAL FITTING ----- #
optparams = optimize.curve_fit(exponential, temp, psi, p0=alkuarvaus, sigma=psiError, absolute_sigma=True)
params2 = optparams[0]
# For error estimates we will also need the off-diagonal covariant 
# matrix elements as the variables are dependent on each other.
popterr2 = np.sqrt(np.diag(optparams[1]))
print('optimal parameteres for the exponential fit are ', params2, '\n and the errors are ', popterr2)

x2 = np.linspace(270, 380, 1000)
expfit = exponential(x2, params2[0], params2[1])

# latent heat based on parameters
l_23 = R*(params2[0] + params2[1]*T_0)
print('Estimated latent heat based on calculated parameters at T = 373 K is ', l_23)

# ----- ERRORS FOR GRAPHING ----- #
# fitYErr = np.sqrt( pCov[0,0] + fitX*(2*pCov[0,1] + pCov[1,1]*fitX ) )
#yerr = np.sqrt(mass*popterr[0] + popterr[1])
#print(linearfit)
#yerr = mass.std() * np.sqrt(1/len(mass) +
#                          (mass - mass.mean())**2 / np.sum((mass - mass.mean())**2))

# %%
# ----- ADDITIONAL ERROR CALCULATIONS ----- #
da=popterr2[0]
db=popterr2[1]
#da=np.sqrt(8.64740931e-01)
#db=np.sqrt(6.68646211e-06)
#sigma=-4.27355624e+01

errors = optparams[1]
sigma = errors[0,1]
# print(sigma)

# using the error propagation law
def error(t):
    dL = np.sqrt((R*da)**2 + (R*t*db)**2 + 2*R**2*t*sigma)
    return dL

# set temperature
T = np.array([273.15, 300, 340, 373.15])
l_23 = R*(params2[0] + params2[1]*T)
# note: printing values with or without molar mass of water
print('latent heat value is: ', l_23/water_mol_mass)
print('latent heat error is: ', error(T)/water_mol_mass)

# %%
# --- K VALUE CALCULATION --- #

# power loss
p_0 = parameters[1]/300.0*1000 # this is in W
p_0_err = popterr[1]/300*1000
d = 0.47*0.1 # meter
# divided by 300 due to 5 min = 300 s
print('p_0 is ', parameters[1]/300.0)
print('p_0 error is ', popterr[1]/300)

### --- radius of the glass ball --- ###
V = 3 # dm^3

R = (3/4 * V/np.pi)**(1/3)*0.1 #  meters
print(R)

### --- calculating the k value for the insulator --- ###
k = (p_0 * d)/(4*np.pi*(R+d)**2*78.1)
print(k)

# k value error
# here I use the law of error propagation
# R+d is used as a single variable with 
# their respective errors summed
# hardcoded :D
# 1st line p_0, 2nd d, 3rd temp, and finally R+d
k_err = np.sqrt( ((d)/(4*np.pi*(R+d)**2*78.1) * p_0_err)**2 \
                + ((p_0)/(4*np.pi*(R+d)**2*78.1)*0.005)**2 \
                + ((p_0 * d)/(4*np.pi*(R+d)**2*78.1**2)*0.5)**2 \
                + ((2*p_0 * d)/(4*np.pi*(R+d)**3*78.1)*0.005)**2) 

# %%
# ----- GRAPHING ----- #
plt.rcParams.update({
    "text.usetex": True
})

fig, ax1 = plt.subplots(dpi=600)
ax1.scatter(mass,  teho, s=5, label='mittauspisteet')
ax1.plot(x1, linearfit, '-', color='orange')
#ax1.fill_between(x1, linearfit - yerr, linearfit + yerr, color='blue', alpha=0.2)
ax1.set_xlabel(r'Tiivistyneen veden massa $m$ [kg]')
ax1.set_ylabel(r'Energia jouleina $Q$ [kJ]')
#fig.legend()
plt.savefig("ekakuvaaja.pdf")
plt.savefig("ekakuvaaja.png")

fig2, ax2 = plt.subplots(dpi=600)
ax2.scatter(temp, psi, s=5)
ax2.plot(x2, expfit, '-', color='orange')
ax2.set_ylabel(r'Ilmanpaine $P$ [Pa]')
ax2.set_xlabel(r'Lämpötila $T$ [K]')
plt.savefig("tokakuvaaja.pdf")
plt.savefig("tokakuvaaja.png")

#plt.show()
# %%
