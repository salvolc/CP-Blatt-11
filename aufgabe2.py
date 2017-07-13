import numpy as np
import itertools 
from matplotlib import pyplot as plt 
import matplotlib.cm as cm
#from mpl_toolkits.mplot3d import Axes3D
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')



















































































# nb = 64
# xA = np.genfromtxt("Aufgabe2bxAN.txt")
# yA = np.genfromtxt("Aufgabe2byAN.txt")

# zg = np.genfromtxt("Aufgabe2a.txt")

# fig, ax1 = plt.subplots()
# ax1.hist(zg,nb,lw=1,linestyle='solid',edgecolor='black',color='blue')
# ax1.tick_params('y', colors='b')
# ax1.set_xlabel('z')
# ax1.set_ylabel('N', color='b')
# plt.title("Gaußverteilte Zahlen mit Box-Muller")

# ax2 = ax1.twinx()
# ax2.set_position([0.125,0.076,0.775,0.803])
# ax2.set_frame_on(False)
# ax2.plot(xA,yA,lw=1,color='r')
# ax2.tick_params('y', colors='r')
# ax2.set_xlabel('x')
# ax2.set_ylabel('f(x)', color='r')
# plt.savefig("Aufgabe2a.pdf")
# plt.close()

# zwg = np.genfromtxt("Aufgabe2b.txt")



# fig, ax1 = plt.subplots()
# ax1.hist(zwg,nb,lw=1,linestyle='solid',edgecolor='black',color='blue')
# ax1.tick_params('y', colors='b')
# ax1.set_xlabel('z')
# ax1.set_ylabel('N', color='b')
# plt.title("Gaußverteilte Zahlen mit zentralen Grenzwertsatz")

# ax2 = ax1.twinx()
# ax2.set_position([0.125,0.076,0.775,0.803])
# ax2.set_frame_on(False)
# ax2.plot(xA,yA,lw=1,color='r')
# ax2.tick_params('y', colors='r')
# ax2.set_xlabel('x')
# ax2.set_ylabel('f(x)', color='r')
# plt.savefig("Aufgabe2b.pdf")
# plt.close()

# ###################################################
# ###################################################
# ###################################################

# zNEU = np.genfromtxt("Aufgabe2c.txt")

# fig, ax1 = plt.subplots()
# ax1.hist(zNEU,nb,lw=1,linestyle='solid',edgecolor='black',color='blue')
# ax1.tick_params('y', colors='b')
# ax1.set_xlabel('z')
# ax1.set_ylabel('N', color='b')
# #plt.title("Gaußverteilte Zahlen mit Neumann")

# ax2 = ax1.twinx()
# ax2.set_position([0.125,0.076,0.775,0.803])
# ax2.set_frame_on(False)
# ax2.plot(xA,yA,lw=1,color='r')
# ax2.tick_params('y', colors='r')
# ax2.set_xlabel('x')
# ax2.set_ylabel('f(x)', color='r')
# plt.savefig("Aufgabe2c.pdf")
# plt.close()



# zVNEU = np.genfromtxt("Aufgabe2c2.txt")
# a = np.array([0.01,0.1,0.5,1.0,1.5,10])

# fig, axes = plt.subplots(2,3,squeeze=False)#,sharey=True)
# fig.set_size_inches(14,7)

# axes[0,0].hist(zVNEU[0],nb,lw=1,linestyle='solid',edgecolor='black',color='blue',normed=True)
# axes[0,0].tick_params('y', colors='b')
# axes[0,0].set_xlabel('z')
# axes[0,0].set_ylabel('N bzw. f(x)', color='b')
# axes[0,0].plot(xA,yA,lw=1,color='r')
# axes[0,0].set_title(r"Gaußverteilte Zahlen mit Neumann und $\alpha = 0.01$")


# plt.title("Gaußverteilte Zahlen mit Neumann")
# axes[0,1].hist(zVNEU[1],nb,lw=1,linestyle='solid',edgecolor='black',color='blue',normed=True)
# axes[0,1].tick_params('y', colors='b')
# axes[0,1].set_xlabel('z')
# axes[0,1].set_ylabel('N bzw. f(x)', color='b')
# axes[0,1].plot(xA,yA,lw=1,color='r')
# axes[0,1].set_title(r"Gaußverteilte Zahlen mit Neumann und $\alpha = 0.1$")



# axes[0,2].hist(zVNEU[2],nb,lw=1,linestyle='solid',edgecolor='black',color='blue',normed=True)
# axes[0,2].tick_params('y', colors='b')
# axes[0,2].set_xlabel('z')
# axes[0,2].set_ylabel('N bzw. f(x)', color='b')
# axes[0,2].plot(xA,yA,lw=1,color='r')
# axes[0,2].set_title(r"Gaußverteilte Zahlen mit Neumann und $\alpha = 1.0$")


# axes[1,0].hist(zVNEU[3],nb,lw=1,linestyle='solid',edgecolor='black',color='blue',normed=True)
# axes[1,0].tick_params('y', colors='b')
# axes[1,0].set_xlabel('z')
# axes[1,0].set_ylabel('N bzw. f(x)', color='b')
# axes[1,0].plot(xA,yA,lw=1,color='r')
# axes[1,0].set_title(r"Gaußverteilte Zahlen mit Neumann und $\alpha = 1.5$")



# axes[1,1].hist(zVNEU[4],nb,lw=1,linestyle='solid',edgecolor='black',color='blue',normed=True)
# axes[1,1].tick_params('y', colors='b')
# axes[1,1].set_xlabel('z')
# axes[1,1].set_ylabel('N bzw. f(x)', color='b')
# axes[1,1].plot(xA,yA,lw=1,color='r')
# axes[1,1].set_title(r"Gaußverteilte Zahlen mit Neumann und $\alpha = 10$")



# axes[1,2].hist(zVNEU[5],nb,lw=1,linestyle='solid',edgecolor='black',color='blue',normed=True)
# axes[1,2].tick_params('y', colors='b')
# axes[1,2].set_xlabel('z')
# axes[1,2].set_ylabel('N bzw. f(x)', color='b')
# axes[1,2].plot(xA,yA,lw=1,color='r')
# axes[1,2].set_title(r"Gaußverteilte Zahlen mit Neumann und $\alpha = 100$")

# plt.tight_layout()

# plt.savefig("Aufgabe2c2.pdf")
# plt.close()


# zinv = np.genfromtxt("Aufgabe2d.txt")
# xANd = np.genfromtxt("Aufgabe2dxAN.txt")
# yANd = np.genfromtxt("Aufgabe2dyAN.txt")
# zinv = zinv[np.absolute(zinv) < 100]


# fig, ax1 = plt.subplots()
# ax1.hist(zinv,60,lw=1,linestyle='solid',edgecolor='black',color='blue')
# ax1.tick_params('y', colors='b')
# ax1.set_xlabel('z')
# ax1.set_ylabel('N', color='b')
# plt.title("Zufallszahlen mit Inversionsmethode")


# ax2 = ax1.twinx()
# ax2.set_position([0.125,0.076,0.775,0.803])
# ax2.set_frame_on(False)
# ax2.plot(xANd,yANd,lw=1,color='r')
# ax2.tick_params('y', colors='r')
# ax2.set_xlabel('x')
# ax2.set_ylabel('f(x)', color='r')
# plt.savefig("Aufgabe2d.pdf")
# plt.close()


# zinvSMD = np.genfromtxt("Aufgabe2dSMD.txt")
# zinvSMD = zinvSMD[np.absolute(zinvSMD) < 100]

# plt.hist(zinvSMD,1500,lw=0.5,linestyle='solid',edgecolor='black',color='blue',normed=True,stacked=True)
# plt.plot(xANd,yANd,color="red")
# plt.xlabel("z")
# plt.ylabel('N bzw. f(x)')
# plt.xlim(-4.5,4.5)
# plt.title("Zufallszahlen mit Inversionsmethode aus SMD")
# plt.savefig("Aufgabe2dSMD.pdf")