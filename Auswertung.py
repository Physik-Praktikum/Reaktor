import numpy as np
import sympy as sy
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as opt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import os.path
from scipy import integrate



def differentieller_Fit(x,a,b):
	return a*(np.sin(np.pi*x/b))**2
def integrale_Kennlinie(x,a,b):
	return ((a*b)/(np.pi*beta))*( ((np.pi*x)/(2*b)) - 0.25*np.sin((2*np.pi*x/b)))

lambdas = np.array([0.0124,0.0305,0.111,0.301,1.14,3.01]) # Zerfallsbreiten der 6 Gruppen verzögerter Neutronen in 1/s
betas = np.array([0.021,0.14,0.126,0.253,0.074,0.027]) # Anteil der verzögerten Neutronen an allen Spaltneutronen
z1 = np.array([799,1601,2401,3203,3994])
z2 = np.array([707,1806,2682,3422,4000])
#z2 = np.array([4000,3422,2682,1806,707])
beta = betas.sum() 
a = betas/(beta) # beta_i/beta für die 6 Gruppen verzögerter Neutronen
l = 0.0001 # Lebensdauer der prompten Neutronen in s

T2 = np.array([175.99,88.83,60.11,63,104.22]) # gemessene Verdopplungszeiten in s
T2_error = np.array([1.43,0.43,0.42,0.81,0.74]) # statistische Unsicherheiten der Verdopplungszeiten in s
Stab1 = np.array([0,799,1601,2401,3203,3994]) #Stabstellungen von Stab 1 in units
Stab2 = np.array([4000,3422,2682,1806,707,0]) #Stabstellungen von Stab 2 in units

Ts = T2/np.log(2) # Stabile Reaktorperioden in s
Ts_error = T2_error/np.log(2) # deren Unsicherheiten in s

rho = np.ones_like(Ts)
for i in range(len(rho)):
	rho[i] = l/Ts[i] + (betas/(1 + lambdas*Ts[i])).sum()
rho_norm = rho/beta
rho_norm_error = np.ones_like(rho_norm)

for i in range(len(rho_norm_error)):
	x = ((betas*lambdas)/((1+lambdas*Ts[i])**2)).sum()
	x = x*Ts_error[i]/beta
	x = x + (l*Ts_error[i])/(beta*Ts[i]**2)
	rho_norm_error[i] = x

Text1 = "Fitansatz:  "
Text1 = Text1 + r"$y = a\cdot\sin^2\left(\frac{\pi\cdot x}{b}\right)$"


p_opt,kov = opt(differentieller_Fit,z1/np.amax(z1),rho_norm/np.amax(rho_norm))#,sigma=rho_norm_error)
z1_fit = np.linspace(np.amin(z1/np.amax(z1)),np.amax(z1/np.amax(z1)),1000)
rho1_fit = differentieller_Fit(z1_fit,*p_opt)
x = np.linspace(0,4000)
diff1 = differentieller_Fit(x,0.120,5440)
int1 = integrate.cumtrapz(diff1, x, initial=0) # kumulierte Intergralwerte
#int1 = integrale_Kennlinie(x,0.120,5440)

Text2 = r"$a = ($"
Text2 = Text2 + "{0:.2f}".format(p_opt[0])
Text2 = Text2 + r"$\pm$"
Text2 = Text2 + "{0:.2f}".format(np.sqrt(kov[0][0]))
Text2 = Text2 + r"$)$"


Text3 = r"$b = ($"
Text3 = Text3 + "{0:.2f}".format(p_opt[1])
Text3 = Text3 + r"$\pm$"
Text3 = Text3 + "{0:.2f}".format(np.sqrt(kov[1][1]))
Text3 = Text3 + r"$)$"


p_opt,kov = opt(differentieller_Fit,z2/np.amax(z2),rho_norm/np.amax(rho_norm))#,sigma=rho_norm_error)
z2_fit = np.linspace(np.amin(z2/np.amax(z2)),np.amax(z2/np.amax(z2)),1000)
rho2_fit = differentieller_Fit(z2_fit,*p_opt)
diff2 = differentieller_Fit(x,0.116,5600)
int2 = integrate.cumtrapz(diff2, x, initial=0) # kumulierte Intergralwerte
#int2 = integrale_Kennlinie(x,0.116,5600)


Text4 = r"$a = ($"
Text4 = Text4 + "{0:.2f}".format(p_opt[0])
Text4 = Text4 + r"$\pm$"
Text4 = Text4 + "{0:.2f}".format(np.sqrt(kov[0][0]))
Text4 = Text4 + r"$)$"


Text5 = r"$b = ($"
Text5 = Text5 + "{0:.2f}".format(p_opt[1])
Text5 = Text5 + r"$\pm$"
Text5 = Text5 + "{0:.2f}".format(np.sqrt(kov[1][1]))
Text5 = Text5 + r"$)$"

#Datei = open("/home/tom/GitHub_Repositories/Reaktor/diff_Kennlinie.txt","a")
#for i in range(len(rho_norm)):
#	print("{0:.0f} & {1:.0f} & {2:.5f} \pm {3:.5f}".format(z1[i],z2[i],rho_norm[i],rho_norm_error[i]),file=Datei)
#Datei.close()



fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.errorbar(z1/np.amax(z1),rho_norm/np.amax(rho_norm),rho_norm_error/np.amax(rho_norm),fmt="o",capsize=5,color="red",label="Messung Stab 1")
ax1.plot(z1_fit,rho1_fit,color="red",label="Fitfunktion Stab 1")
ax1.errorbar(z2/np.amax(z2),rho_norm/np.amax(rho_norm),rho_norm_error/np.amax(rho_norm),fmt="o",capsize=5,color="green",label="Messung Stab 2")
ax1.plot(z2_fit,rho2_fit,color="green",label="Fitfunktion Stab 2")
ax1.tick_params(axis='both', which='major', width=1.5,length=10,direction="in", labelsize=14)
ax1.tick_params(axis='both', which='minor', width=1.5, length=4,direction="in", labelsize=8)
ax1.tick_params(which="both", top=True,labeltop=True,right=True,labelright=True)
ax1.xaxis.set_major_locator(MultipleLocator(0.1))
ax1.xaxis.set_minor_locator(MultipleLocator(0.02))
ax1.yaxis.set_major_locator(MultipleLocator(0.2))
ax1.yaxis.set_minor_locator(MultipleLocator(0.04))
ax1.set_xlim(min(np.amin(z1/np.amax(z1)),np.amin(z2/np.amax(z2))),max(np.amax(z1/np.amax(z1)),np.amax(z2/np.amax(z2))))
ax1.set_ylim(min(np.amin(rho_norm/np.amax(rho_norm)), np.amin(rho1_fit), np.amin(rho2_fit)),max(np.amax(rho_norm/np.amax(rho_norm)), np.amax(rho1_fit), np.amax(rho2_fit))) 
ax1.set_xlabel(r"$\frac{z}{z_{\mathrm{max.}}}$",fontsize=20)
ax1.set_ylabel(r"$\left.\frac{\mathrm{d}\rho'}{\mathrm{d}z}\right/ \frac{\mathrm{d}\rho'_{\mathrm{max}}}{\mathrm{d}z}$",fontsize=20)
ax1.annotate(Text1,(0.5,0.7),fontsize=20,color="black")
ax1.annotate(Text2,(0.5,0.6),fontsize=17,color="red")
ax1.annotate(Text3,(0.5,0.55),fontsize=17,color="red")
ax1.annotate(Text4,(0.5,0.45),fontsize=17,color="green")
ax1.annotate(Text5,(0.5,0.4),fontsize=17,color="green")
ax1.grid(True)
ax1.legend()


fig2 = plt.figure()
ax11 = fig2.add_subplot(121)
ax12 = fig2.add_subplot(122)
ax11.plot(x,diff1,color="red",label="Stab 1")
ax11.plot(x,diff2,color="green",label="Stab 2")
ax12.plot(x,int1/10,color="red",label="Stab 1")
ax12.plot(x,int2/10,color="green",label="Stab 2")
#ax11.set_title("Differentielle Kennlinien",fontsize=14)
#ax12.set_title("Integrale Kennlinien",fontsize=14)
ax11.set_xlabel(r"$z/\mathrm{units}$",fontsize=14)
ax12.set_xlabel(r"$z/\mathrm{units}$",fontsize=14)
ax11.set_ylabel(r"$\left.\frac{\mathrm{d}\rho'}{\mathrm{d}z}\right/\frac{\mathrm{ct}}{\mathrm{units}}$",fontsize=14)
ax12.set_ylabel(r"$\rho'/\mathrm{ct}$",fontsize=14)
ax11.grid(True)
ax12.grid(True)
ax11.legend(fontsize=14)
ax12.legend(fontsize=14)
ax11.tick_params(axis='both', which='major', width=1.5,length=10,direction="in", labelsize=14)
ax11.tick_params(axis='both', which='minor', width=1.5, length=4,direction="in", labelsize=8)
ax11.tick_params(which="both", top=True,labeltop=True,right=True,labelright=True)
ax11.xaxis.set_major_locator(MultipleLocator(500))
ax11.xaxis.set_minor_locator(MultipleLocator(100))
ax11.yaxis.set_major_locator(MultipleLocator(0.01))
ax11.yaxis.set_minor_locator(MultipleLocator(0.002))
ax12.tick_params(axis='both', which='major', width=1.5,length=10,direction="in", labelsize=14)
ax12.tick_params(axis='both', which='minor', width=1.5, length=4,direction="in", labelsize=8)
ax12.tick_params(which="both", top=True,labeltop=True,right=True,labelright=True)
ax12.xaxis.set_major_locator(MultipleLocator(500))
ax12.xaxis.set_minor_locator(MultipleLocator(100))
ax12.yaxis.set_major_locator(MultipleLocator(2.5))
ax12.yaxis.set_minor_locator(MultipleLocator(0.5))


plt.show()
