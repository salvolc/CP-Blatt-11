import numpy as np
from matplotlib import pyplot as plt
import scipy.integrate as integrate
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
import matplotlib.cm as cm


n_particles = 1
dim = 2
maxframe = 1200
secs = maxframe/60
h = 0.01
sts = secs/h

fig = plt.figure()
lim = 8

ax = plt.axes(xlim=(-100,100), ylim=(0,1))

line, = ax.plot([], [], lw=2)
plotlays, plotcols = [n_particles], cm.rainbow(np.linspace(0, 1, n_particles))
lines = []



for index in range(n_particles):
    lobj = ax.plot([],[],'-',lw=1,color=plotcols[index])[0]
    lines.append(lobj)

def init():
    for line in lines:
        line.set_data([],[])
    return lines 


xlist = []
ylist = []
 
print("\nAnimation: Einlesen der Daten..")
y = np.genfromtxt("aufg1d.txt",unpack=False)
x = np.linspace(-100,100,200)

steps = int(len(y)/n_particles)
for i in range(n_particles):
	xlist.append(x)
	ylist.append(y[i*steps:(i+1)*steps])


print("Animation: Animiere nun! \nAnimation: Dies könnte ein Weilchen dauern")

def animate(i):
	maxval = int(i*(steps/maxframe))
	for lnum,line in enumerate(lines):
		line.set_data(x, ylist[lnum][maxval])

	if(maxval > steps):
		ax.plot([-100,-1/2*10,-1/2*10,1/2*10,1/2*10,100],[0,0,10,10,0,0],'r',label="Barriere")

	print("Animation: Frame ",i," von ",maxframe, " also ",(i/maxframe)*100,"%",end="\r")
	return lines

anim = animation.FuncAnimation(fig, animate, init_func=init,frames=maxframe, interval=10, blit=True)

anim.save('animation.mp4', fps=60, extra_args=['-vcodec', 'libx264'])
print("Animation: Fertig \n")