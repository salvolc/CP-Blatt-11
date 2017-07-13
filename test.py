import numpy as np
import itertools 
from matplotlib import pyplot as plt 
import matplotlib.cm as cm

pi = 3.141

def f(y):
	return -1/np.sqrt(4*pi*y**3 - 4 * pi**2 * y**4)

def p(x):
	return 1/pi * 1/(1+x**2)

x = np.linspace(0.01,1/pi-0.01,10000)

plt.plot(x,np.absolute(f(x))/300)
plt.plot(-x,np.absolute(f(x))/300)

y = np.linspace(-1,1,1000)
plt.plot(y/100,p(y))

plt.savefig("test2.pdf")