import numpy as np
import matplotlib.pyplot as plt
from iminuit import Minuit
import scipy.stats as stats
import random


#################### LOADING DATA

data = np.genfromtxt("Bigrun.phsp", dtype='unicode')

events = data[:,3].astype(int)

dim = len(np.unique(events))

allparticles = data[:,4]
allenergies = data[:,8].astype(float)
alltimes = data[:,9].astype(float)

ignoreindices = []
for n in range(len(events)-1):
    if events[n] == events[n+1]:
        ignoreindices.append(n+1)

particles = np.delete(allparticles, ignoreindices)
energies = np.delete(allenergies, ignoreindices)
times = np.delete(alltimes, ignoreindices)


##################### GUESSING ENERGY FROM FLIGHT TIME

velocigies = (1 / np.sqrt(1. - ((0.5025 / times)/299792458)**2) -1) * 938.27

plt.figure(1)
plt.minorticks_on()
plt.grid(True, which='both')
plt.grid(linestyle='dashed', linewidth=0.25, which="minor")
plt.scatter(energies[particles=="proton"], velocigies[particles=="proton"], s=2, alpha=0.5, color="red", label="Protons")
plt.scatter(energies[particles=="neutron"], velocigies[particles=="neutron"], s=2, alpha=0.5, color="blue", label="Neutrons")
plt.title("Geant4 energy vs energy from time of flight, 100MeV incident neutrons")
plt.legend()
plt.xlabel("Geant4 energy (MeV)")
plt.ylabel("Flight energies from TOF (MeV)")
#plt.show()



######################### CALCULATING INCIDENT ENERGY FROM ANGLE

a = np.pow(velocigies/938.27, 2) * np.cos(np.atan(5/50))**2  - np.cos(np.atan(5/50))**2 - np.pow(velocigies/938.27, 2)
b = 2 * velocigies/938.27
c = - 1 * np.cos(np.atan(5/50))**2 * (np.pow(velocigies/938.27, 2) -1)

incidents = (-b - np.sqrt(np.pow(b,2) - 4 * a * c)) / (2*a) * 938.27 - 938.27


totals, edges = np.histogram(incidents, bins=np.linspace(0, 100, 100))
protons, edges = np.histogram(incidents[particles=="proton"], bins=np.linspace(0, 100, 100))
neutrons, edges = np.histogram(incidents[particles=="neutron"], bins=np.linspace(0, 100, 100))


plt.figure(2)
plt.minorticks_on()
plt.grid(True, which='both')
plt.grid(linestyle='dashed', linewidth=0.25, which="minor")
plt.stairs(totals, edges, label="Total", color="black")
plt.stairs(protons, edges, label="Protons", color="red")
plt.stairs(neutrons, edges, label="Neutrons", color="blue")
plt.title("Calculated incident energy")




########################## ML FITTING

def totallogl(mean, deviation): # The data "means" is hardcoded into the function. Usually this is not a proble, since data rarely changes after the experiment :p
    return -np.sum(- 0.5 * np.log(np.pow(deviation, 2)) - np.pow(incidents[incidents > 80] - mean, 2) / (2 * np.pow(deviation, 2)))

neutronincidents = incidents[particles == "neutron"]

def neutronlogl(mean, deviation): # The data "means" is hardcoded into the function. Usually this is not a proble, since data rarely changes after the experiment :p
    return -np.sum(- 0.5 * np.log(np.pow(deviation, 2)) - np.pow(neutronincidents[neutronincidents > 90] - mean, 2) / (2 * np.pow(deviation, 2)))


print("Maximum likelihood on total data, not taking uncertainties into account")
totalminu = Minuit(totallogl, mean=10, deviation= 1)
totalminu.errordef = 0.5 #Because of log-likelihood
totalminu.limits = [(None, None), (0, None)]

totalminu.migrad()
totalminu.hesse() # Finite difference (assume parabolic) covariance matrix
totalminu.minos() # Profile likelihood method
#print(minu) # Printing all out
print(repr(totalminu.values[:])) # Accessing the fitted parameters
#print(repr(minu.covariance)) # Accessing the covariance matrix from Hesse
#print(repr(totalminu.merrors[0].lower), repr(totalminu.merrors[0].upper), repr(totalminu.merrors[1].lower), repr(totalminu.merrors[1].upper)) # Accessing the assymmetrical errors from MINOS. Haven't found a better way yet :/

plt.plot(np.linspace(0, 100, 1000), len(incidents[incidents > 80]) * (edges[1] - edges[0]) * stats.norm.pdf(np.linspace(0, 100, 1000), loc=totalminu.values[0], scale=totalminu.values[1]), color="black", label="Fit")

print("Maximum likelihood on neutron data, not taking uncertainties into account")
neutronminu = Minuit(neutronlogl, mean=10, deviation= 1)
neutronminu.errordef = 0.5 #Because of log-likelihood
neutronminu.limits = [(None, None), (0, None)]

neutronminu.migrad()
neutronminu.hesse() # Finite difference (assume parabolic) covariance matrix
neutronminu.minos() # Profile likelihood method
#print(minu) # Printing all out
print(repr(neutronminu.values[:])) # Accessing the fitted parameters
#print(repr(minu.covariance)) # Accessing the covariance matrix from Hesse
#print(repr(neutronminu.merrors[0].lower), repr(neutronminu.merrors[0].upper), repr(neutronminu.merrors[1].lower), repr(neutronminu.merrors[1].upper)) # Accessing the assymmetrical errors from MINOS. Haven't found a better way yet :/

plt.plot(np.linspace(0, 100, 1000), len(neutronincidents[neutronincidents > 90]) * (edges[1] - edges[0]) * stats.norm.pdf(np.linspace(0, 100, 1000), loc=neutronminu.values[0], scale=neutronminu.values[1]), color="blue", label="Fit")
plt.legend()
#plt.show()




#################################### BOOTSTRAPPING!

#particles = np.delete(allparticles, ignoreindices)
#energies = np.delete(allenergies, ignoreindices)
#times = np.delete(alltimes, ignoreindices)



bootnumber = 1000000

protontimes = times[particles=="proton"]
bootprotontimes = np.zeros(bootnumber)
for n in range(bootnumber):
    index = random.randint(0,len(protontimes) -1)
    bootprotontimes[n] = protontimes[index] + np.random.normal(0, 100E-12) # 100ps time uncertainty
bootprotonvelocigies = (1 / np.sqrt(1. - ((0.5025 / bootprotontimes)/299792458)**2) -1) * 938.27
a = np.pow(bootprotonvelocigies/938.27, 2) * np.cos(np.atan(5/50))**2  - np.cos(np.atan(5/50))**2 - np.pow(bootprotonvelocigies/938.27, 2)
b = 2 * bootprotonvelocigies/938.27
c = - 1 * np.cos(np.atan(5/50))**2 * (np.pow(bootprotonvelocigies/938.27, 2) -1)
bootprotonincidents = (-b - np.sqrt(np.pow(b,2) - 4 * a * c)) / (2*a) * 938.27 - 938.27
bootprotons, edges = np.histogram(bootprotonincidents, bins=np.linspace(0, 125, 250))

neutrontimes = times[particles=="neutron"]
bootneutrontimes = np.zeros(bootnumber)
for n in range(bootnumber):
    index = random.randint(0,len(neutrontimes) -1)
    bootneutrontimes[n] = neutrontimes[index] + np.random.normal(0, 100E-12) # 100ps time uncertainty
bootneutronvelocigies = (1 / np.sqrt(1. - ((0.5025 / bootneutrontimes)/299792458)**2) -1) * 938.27
a = np.pow(bootneutronvelocigies/938.27, 2) * np.cos(np.atan(5/50))**2  - np.cos(np.atan(5/50))**2 - np.pow(bootneutronvelocigies/938.27, 2)
b = 2 * bootneutronvelocigies/938.27
c = - 1 * np.cos(np.atan(5/50))**2 * (np.pow(bootneutronvelocigies/938.27, 2) -1)
bootneutronincidents = (-b - np.sqrt(np.pow(b,2) - 4 * a * c)) / (2*a) * 938.27 - 938.27
bootneutrons, edges = np.histogram(bootneutronincidents, bins=np.linspace(0, 125, 250))


plt.figure(3)
plt.minorticks_on()
plt.grid(True, which='both')
plt.grid(linestyle='dashed', linewidth=0.25, which="minor")
#plt.stairs(totals, edges, label="Total", color="black")
plt.stairs(bootprotons, edges, label="Protons", color="red")
plt.stairs(bootneutrons, edges, label="Neutrons", color="blue")
#plt.stairs(neutrons, edges, label="Neutrons", color="blue")
plt.title("Bootstrapped incident energy")
plt.legend()




#### ML FITTING ON THE BOOTSTRAP


def bootprotonlogl(mean, deviation): # The data "means" is hardcoded into the function. Usually this is not a proble, since data rarely changes after the experiment :p
    return -np.sum(- 0.5 * np.log(np.pow(deviation, 2)) - np.pow(bootprotonincidents[bootprotonincidents > 70] - mean, 2) / (2 * np.pow(deviation, 2)))

def bootneutronlogl(mean, deviation): # The data "means" is hardcoded into the function. Usually this is not a proble, since data rarely changes after the experiment :p
    return -np.sum(- 0.5 * np.log(np.pow(deviation, 2)) - np.pow(bootneutronincidents[bootneutronincidents > 80] - mean, 2) / (2 * np.pow(deviation, 2)))


print("Maximum likelihood on bootstrapped protons")
totalminu = Minuit(bootprotonlogl, mean=10, deviation= 1)
totalminu.errordef = 0.5 #Because of log-likelihood
totalminu.limits = [(None, None), (0, None)]

totalminu.migrad()
totalminu.hesse() # Finite difference (assume parabolic) covariance matrix
totalminu.minos() # Profile likelihood method
#print(minu) # Printing all out
print(repr(totalminu.values[:])) # Accessing the fitted parameters
#print(repr(minu.covariance)) # Accessing the covariance matrix from Hesse
#print(repr(totalminu.merrors[0].lower), repr(totalminu.merrors[0].upper), repr(totalminu.merrors[1].lower), repr(totalminu.merrors[1].upper)) # Accessing the assymmetrical errors from MINOS. Haven't found a better way yet :/

plt.plot(np.linspace(0, 125, 1000), len(bootprotonincidents[bootprotonincidents > 70]) * (edges[1] - edges[0]) * stats.norm.pdf(np.linspace(0, 125, 1000), loc=totalminu.values[0], scale=totalminu.values[1]), color="red", label="Fit", linestyle="dashed")

print("Maximum likelihood on bootstrapped neutrons")
neutronminu = Minuit(bootneutronlogl, mean=10, deviation= 1)
neutronminu.errordef = 0.5 #Because of log-likelihood
neutronminu.limits = [(None, None), (0, None)]

neutronminu.migrad()
neutronminu.hesse() # Finite difference (assume parabolic) covariance matrix
neutronminu.minos() # Profile likelihood method
#print(minu) # Printing all out
print(repr(neutronminu.values[:])) # Accessing the fitted parameters
#print(repr(minu.covariance)) # Accessing the covariance matrix from Hesse
#print(repr(neutronminu.merrors[0].lower), repr(neutronminu.merrors[0].upper), repr(neutronminu.merrors[1].lower), repr(neutronminu.merrors[1].upper)) # Accessing the assymmetrical errors from MINOS. Haven't found a better way yet :/

plt.plot(np.linspace(0, 125, 1000), len(bootneutronincidents[bootneutronincidents > 80]) * (edges[1] - edges[0]) * stats.norm.pdf(np.linspace(0, 125, 1000), loc=neutronminu.values[0], scale=neutronminu.values[1]), color="blue", label="Fit", linestyle="dashed")






plt.legend()
plt.show()

