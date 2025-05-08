#%% 

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import constants

# constants
T_0 = 373.15 # boiling point in kelvins at regular air pressure
P_0 = 773*(101325/760) # ambient air pressure in pascals
water_mol_mass = 18.01528 # grams per mol


# ----- IMPORTING DATA -------
teho, mass = np.genfromtxt("data1.txt", comments='#', unpack=True, missing_values='\n')
psi, temp = np.genfromtxt("data2.txt", comments='#', missing_values='\n', unpack=True)

mass = mass*0.001 # conversion to kg

# converting kWh to si-units
# kilo watt = 1000*watt
# watt = Joule per second
# 1 hour = 3600 seconds
# 1000 * 3600 * measured value = sth watt seconds from which we get joules

teho = teho*3600

# converting mmHg to pascal
# 1 mmHg = 101325/760 Pa

psi = psi*(101325/760)

# ---- ADDING YERROR ---- #
psiError = np.zeros_like(psi)+(0.5*(101325/760))
alkuarvaus = np.array( [1000.0, -10.0])

# ------ DEFINING FUNCTIONS --------
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


# ---- LINEAR FITTING ------
popt = optimize.curve_fit(linear, xdata=mass, ydata=teho)
parameters = popt[0]
popterr = np.sqrt(np.diag(popt[1]))
print('optimal parameteres for the linear fit are ', parameters, '\n and the errors are ', popterr)

x1 = np.linspace(0, 0.001*50, 2000)

linearfit = linear(x1, parameters[0], parameters[1])


# ----- EXPONENTIAL FITTING -------
optparams = optimize.curve_fit(exponential, temp, psi, p0=alkuarvaus, sigma=psiError, absolute_sigma=True)
params2 = optparams[0]
# For error estimates we need the off-diagonal covariant matrix elements as the variables are dependent. 
# Here we can use np.fliplr() in order to get off-diag elements to diagonal.
# flip = np.fliplr(optparams[1])
# popterr2 = np.sqrt(np.abs(np.diag(flip)))
# The actual error will sum both the original diagonals and the off-diagonals together.
popterr2 = np.sqrt(np.diag(optparams[1]))
print('optimal parameteres for the exponential fit are ', params2, '\n and the errors are ', popterr2)
# print('flipped matrix ', flip)
# print(np.diag(optparams[1]))
# print(np.linalg.cond(popt[1])) # these two for testing how good the fits are

x2 = np.linspace(270, 380, 1000)
expfit = exponential(x2, params2[0], params2[1])

print('printing the full fit matrix: ', optparams)

# latent heat based on parameters
l_23 = constants.gas_constant*(params2[0] + params2[1]*T_0)
print('Estimated latent heat based on calculated parameters is ', l_23)

# ----- ERRORS FOR GRAPHING ----- #
# fitYErr = np.sqrt( pCov[0,0] + fitX*(2*pCov[0,1] + pCov[1,1]*fitX ) )
#yerr = np.sqrt(mass*popterr[0] + popterr[1])
#print(linearfit)
#yerr = mass.std() * np.sqrt(1/len(mass) +
#                          (mass - mass.mean())**2 / np.sum((mass - mass.mean())**2))

# %%
# ----- GRAPHING -------
plt.rcParams.update({
    "text.usetex": True
})

fig, ax1 = plt.subplots(dpi=600)
ax1.scatter(mass,  teho, s=5, label='mittauspisteet')
ax1.plot(x1, linearfit, '-', color='orange')
#ax1.fill_between(x1, linearfit - yerr, linearfit + yerr, color='blue', alpha=0.2)
ax1.set_xlabel(r'Tiivistyneen veden massa $m$ [kg]')
ax1.set_ylabel(r'Energia jouleina $Q$ [J]')
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
# ----- ADDITIONAL ERROR CALCULATIONS -----
R=constants.gas_constant
da=popterr2[0]
db=popterr2[1]
#da=np.sqrt(8.64740931e-01)
#db=np.sqrt(6.68646211e-06)
sigma=-4.27355624e+01

def error(t):
    dL = np.sqrt((R*da)**2 + (R*t*db)**2 + 2*R**2*t*sigma)
    return dL

# print('list of temperatures: ',temp)
T = 373.15
print('molar latent heat error is: ', error(T)/water_mol_mass)
# print('latent heat value: ', R*(params2[0]+params2[1]*temp))
l_23 = constants.gas_constant*(params2[0] + params2[1]*T)
print(l_23/(water_mol_mass))

# %%

print('p_0 is ', parameters[1]/300.0)
print('p_0 error is ', popterr[1]/300)

# radius of the glass ball
V = 3 # dm^3

R = (3/4 * V/np.pi)**(1/3)

# %%
