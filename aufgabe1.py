import numpy as np
import itertools 
from matplotlib import pyplot as plt 
import matplotlib.cm as cm






















































































# r1 = np.genfromtxt("Aufgabe1i.txt")


# counts,bin_edges,a = plt.hist(r1,20,lw=1,linestyle='solid',edgecolor='black',color='blue')
# bin_centres = (bin_edges[:-1] + bin_edges[1:])/2.
# err = np.sqrt(counts)
# plt.errorbar(bin_centres, counts, yerr=err, fmt='.', color='k')

# plt.xlabel("r")
# plt.ylabel("N")
# plt.title("Linear kongruenter Zufallszahlen i)")

# plt.savefig("Aufg1i.pdf")
# plt.close()



# r2 = np.genfromtxt("Aufgabe1ii.txt")

# from matplotlib import gridspec
# f = plt.figure(figsize=(8, 6)) 
# gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
# ax1 = plt.subplot(gs[0])
# ax2 = plt.subplot(gs[1])

# counts,bin_edges,a = ax1.hist(r2,20,lw=1,linestyle='solid',edgecolor='black',color='blue')
# bin_centres = (bin_edges[:-1] + bin_edges[1:])/2.
# err = np.sqrt(counts)
# ax1.errorbar(bin_centres, counts, yerr=err, fmt='.', color='k')
# ax1.set_title('Linear kongruenter Zufallszahlen ii)')
# ax2.axhline(y = np.sum(counts)/20, linewidth=4, color='r')
# ax2.errorbar(bin_centres, counts, yerr=err, fmt='o', color='k')
# ax2.set_ylim(min(counts)-max(err),max(counts)+max(err))
# f.subplots_adjust(hspace=0)
# plt.tight_layout()
# plt.savefig("Aufg1ii.pdf")
# plt.close()




# rM = np.genfromtxt("Aufgabe1Mars.txt")


# from matplotlib import gridspec
# f = plt.figure(figsize=(8, 6)) 
# gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
# ax1 = plt.subplot(gs[0])
# ax2 = plt.subplot(gs[1])

# counts,bin_edges,a = ax1.hist(rM,20,lw=1,linestyle='solid',edgecolor='black',color='blue')
# bin_centres = (bin_edges[:-1] + bin_edges[1:])/2.
# err = np.sqrt(counts)
# ax1.errorbar(bin_centres, counts, yerr=err, fmt='.', color='k')
# ax1.set_title('Marsenne Twister 19937 Zufallszahlen')
# ax2.axhline(y = np.sum(counts)/20, linewidth=4, color='r')
# ax2.errorbar(bin_centres, counts, yerr=err, fmt='o', color='k')
# ax2.set_ylim(min(counts)-max(err),max(counts)+max(err))
# f.subplots_adjust(hspace=0)
# plt.tight_layout()
# plt.savefig("Aufg1M.pdf")
# plt.close()


# #################################################################
# #################################################################
# #################################################################
# #################################################################
# #################################################################
# x = np.genfromtxt("Aufgabe1xAN.txt")
# y = np.genfromtxt("Aufgabe1yAN.txt")



# plt.plot(x,y,lw=1)
# plt.xlabel("r")
# plt.ylabel("N")

# plt.savefig("Aufg1cAN.pdf")
# plt.close()

# nb = 20
# bnorm = False

# z = np.genfromtxt("Aufgabe1c1.txt")
# z = (z[np.logical_not(np.isnan(z))])

# fig, ax1 = plt.subplots()
# ax1.hist(z,nb,lw=1,linestyle='solid',edgecolor='black',color='blue',normed=bnorm)
# ax1.tick_params('y', colors='b')
# ax1.set_xlabel('z')
# ax1.set_ylabel('N', color='b')
# plt.title("Neumann mit linear kongruenter Zufallszahlen i)")

# ax2 = ax1.twinx()
# ax2.set_position([0.125,0.076,0.775,0.803])
# ax2.set_frame_on(False)
# ax2.plot(x,y,lw=1,color='r')
# ax2.tick_params('y', colors='r')
# ax2.set_xlabel('x')
# ax2.set_ylabel('p(x)', color='r')
# plt.savefig("Aufg1c1.pdf")
# plt.close()



# z = np.genfromtxt("Aufgabe1c2.txt")
# z = (z[np.logical_not(np.isnan(z))])

# fig, ax1 = plt.subplots()
# ax1.hist(z,nb,lw=1,linestyle='solid',edgecolor='black',color='blue',normed=bnorm)
# ax1.tick_params('y', colors='b')
# ax1.set_xlabel('z')
# ax1.set_ylabel('N', color='b')
# plt.title("Neumann mit linear kongruenter Zufallszahlen ii)")

# ax2 = ax1.twinx()
# ax2.set_position([0.125,0.076,0.775,0.803])
# ax2.set_frame_on(False)
# ax2.plot(x,y,lw=1,color='r')
# ax2.tick_params('y', colors='r')
# ax2.set_xlabel('x')
# ax2.set_ylabel('p(x)', color='r')


# plt.savefig("Aufg1c2.pdf")
# plt.close()



# z = np.genfromtxt("Aufgabe1c3.txt")
# z = (z[np.logical_not(np.isnan(z))])

# fig, ax1 = plt.subplots()

# ax1.hist(z,nb,lw=1,linestyle='solid',edgecolor='black',color='blue',normed=bnorm)
# ax1.tick_params('y', colors='b')
# ax1.set_xlabel('z')
# ax1.set_ylabel('N', color='b')
# plt.title("Neumann mit Mersenne Twister Zufallszahlen")

# ax2 = ax1.twinx()
# ax2.set_position([0.125,0.076,0.775,0.803])
# ax2.set_frame_on(False)
# ax2.plot(x,y,lw=1,color='r')
# ax2.tick_params('y', colors='r')
# ax2.set_xlabel('x')
# ax2.set_ylabel('p(x)', color='r')
# plt.savefig("Aufg1c3.pdf")
# plt.close()


# #########################################################
# #########################################################
# #########################################################
# #########################################################

# #n = 100000
# dpi = 300
# plt.plot(r1[:-1],r1[1:],'b.',markersize=2)
# plt.xlabel(r"$z_n$")
# plt.ylabel(r"$z_{n+1}$")

# plt.savefig("Aufg1di.png",dpi=dpi)
# plt.close()


# plt.plot(r2[:-1],r2[1:],'b.',markersize=0.2)
# plt.xlabel(r"$z_n$")
# plt.ylabel(r"$z_{n+1}$")

# plt.savefig("Aufg1dii.png",dpi=dpi)
# plt.close()


# plt.plot(rM[:-1],rM[1:],'b.',markersize=0.2)
# plt.xlabel(r"$z_n$")
# plt.ylabel(r"$z_{n+1}$")

# plt.savefig("Aufg1dM.png",dpi=dpi)
# plt.close()


# ###############################################
# ###############################################
# ###############################################

# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# ax.set_xlabel(r'$z_i$')
# ax.set_ylabel(r'$z_{i+1}$')
# ax.set_zlabel(r'$z_{i+2}$')
# ax.scatter(r1[:-2],r1[1:-1],r1[2:],s=4)

# ax.view_init(27, -60)
# plt.savefig("Aufg1ei.png",dpi=dpi)
# plt.close()

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# ax.set_xlabel(r'$z_i$')
# ax.set_ylabel(r'$z_{i+1}$')
# ax.set_zlabel(r'$z_{i+2}$')
# ax.scatter(r2[:-2],r2[1:-1],r2[2:],s=2)

# ax.view_init(20, -128)
# plt.savefig("Aufg1eii.png",dpi=dpi)
# plt.close()

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# ax.set_xlabel(r'$z_i$')
# ax.set_ylabel(r'$z_{i+1}$')
# ax.set_zlabel(r'$z_{i+2}$')
# ax.scatter(rM[:-2],rM[1:-1],rM[2:],s=2)

# ax.view_init(27, -60)
# plt.savefig("Aufg1eM.png",dpi=dpi)
# plt.close()


