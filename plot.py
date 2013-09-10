#!/usr/bin/env python
#
# Visualization of planet movement via numpy and matplotlib
#

import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

# load position data and reshape it
Ndays = 0
N = 0
Nplanets = 0
Tstep = 0
V = np.array([])
W = np.array([])
X = np.array([])

filehandle = open("run.dat", "r")
section = ""
it = 0

for line in filehandle:
	# skip comments
	if line[0] == "#":
		continue
	# determine section	
	if line[0] in ("*"):
		section = line[1]
		it = 0
		continue
	# parameters
	if section == "P":
		temp = line.split()
		Ndays = int(temp[0])
		N = int(temp[1])
		Nplanets = int(temp[2])
		Tstep = float(temp[3])
		# resize array structures
		V.resize((Ndays,N,3))
		W.resize((Ndays,3))
		X.resize((Ndays,N,3))
		continue
	# abs(v)
	if section == "V":
		temp = line.split()
		V[it] = np.array(temp).reshape((N,3))
		it += 1
		continue
	# energy
	if section == "W":
		temp = line.split()
		W[it,:] = temp	
		it += 1
		continue
	# positions
	if section == "X":
		temp = line.split()
		X[it] = np.array(temp).reshape((N,3))
		it += 1
		continue

filehandle.close()

# X shape [Ndays,N,3] => [N,Ndays,3]
X = np.swapaxes(X, 0,1)

# calculate abs(V) and store in V[i,j,3]
V.resize((Ndays,N,4));
for i in range(Ndays):
	for j in range(N):
		V[i,j,3] = np.linalg.norm( V[i,j,0:3]  )

# use real planet names and make up names for probes
n = ["sun", "mercury", "venus", "earth", "mars", "jupiter", "saturn", "uranus", "neptune", "pluto"]
n_probe = ["probe" + str(i) for i in range(N-Nplanets)]
n += n_probe

print str(N) + " objects found."

# set 3d figure
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d", axisbg='black')

# choose different colors for each trajectory
color_planets = plt.cm.jet(np.linspace(0, 1, Nplanets))
color_probes = plt.cm.binary(np.linspace(0, 0.6, N-Nplanets))
colors = np.vstack((color_planets, color_probes))
#colors = ("#ffcc00', 

# set trajectories, only lines get labels
lines  = sum( [ax.plot([], [], [], '-', color=c, label=l) for c, l in zip(colors, n)], [] )
points = sum( [ax.plot([], [], [], 'o', color=c) for c in colors], [] )

# set plot layout 
limit = (-1,1)			# axes limits
#limit = (-20,20)
ax.set_xlim(limit)
ax.set_ylim(limit)
ax.set_zlim(limit)
ax.axis("off") 				# disable axes
# put legend to the right of the plot
#ax.legend(loc="center left", bbox_to_anchor=(1., 0.5), prop={"size":12})	
ax.legend(loc="center right", bbox_to_anchor=(1., 0.5), prop={"size":12})	
counter = plt.figtext(0.1, 0.9, "-", color="white")	# prepare text window

# set point-of-view: theta and phi in degrees
ax.view_init(90, 0) # 60, 20

# function is called each frame 
def anim_sequence(i):
	# fast forward
	#i *= 5
	# set counter
	counter.set_text("t = " + str(i) + " d")
	# set trajectories
	for line, point, pos in zip(lines, points, X):
		x, y, z = pos[:i].T

		line.set_data(x, y)
		line.set_3d_properties(z)

		point.set_data(x[-1:], y[-1:])
		point.set_3d_properties(z[-1:])
		
		fig.canvas.draw()

	return lines + points


# start animator
#anim = animation.FuncAnimation(fig, anim_sequence, frames=Ndays, interval=30, blit=True)

# or simply draw the last frame to show full data
anim_sequence(Ndays)

# save animation as mp4, ffmpeg needed
#anim.save("test.mp4", fps=30), extra_args=["-vcodec", "libx264"])

# show plot with tight layout
plt.tight_layout()
plt.show()
