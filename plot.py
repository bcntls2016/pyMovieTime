#!/users/p1039/coppens/_meuk/python-2.7.12-intel/bin/python
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
from params import *
from plotsettings import *

### SOME GLOBAL SETTINGS 
plt.rcParams['text.color'] = 'lightgrey'
plt.rcParams['axes.labelcolor'] = 'lightgrey'
plt.rcParams['axes.facecolor'] = "black"
plt.rcParams['axes.edgecolor'] = '#A9A9A9'
plt.rcParams['grid.color'] = '#A9A9A9'
plt.rcParams['xtick.color'] = '#A9A9A9'
plt.rcParams['ytick.color'] = '#A9A9A9'
#plt.rc('text', usetex=True)
dpiout = 300

### LABEL PLACEMENT
timex = -xmax + 2.5
timey = ymax - 4.5
fframex = -xmax + 2.5
fframey = -ymax + 2.5

posx = xmax - 22.5
posy = timey
speedvxx = posx
speedvxy = posy - 3.25
speedvyx = speedvxx
speedvyy = speedvxy - 3.25
energyx = speedvyx
energyy = speedvyy - 3.25
particlex = energyx
particley = energyy - 3.25

### FILENAMES AND LABELS
infilename = sys.argv[1]
outfile = 'density.png'
timestring = "$t=" + str(time) + "\,\mathrm{ps}$"

### READ RAW DATA
X = []						# Declare empty list
X.append([])				# Declare new empty list as a first entry
Y = []
Y.append([])
D = []
D.append([])
Cx = []
Cx.append([])
Cy = []
Cy.append([])
blockID = 0					# Define data block counter and set to zero
for line in open(infilename):		# Open file and start reading line by line
	if not line.strip():	# Check if empty line
		X.append([])		# Add a new empty list to the list
		Y.append([])
		D.append([])
		Cx.append([])
		Cy.append([])
		blockID += 1			# Increment current data block by 1
	else:						# Read data in current block into the nested list
		columns = line.split()	# Split line into columns
		X[blockID].append(float(columns[0]))	# Extend the list in the current block by
		Y[blockID].append(float(columns[1]))	# one entry and store the data in it
		D[blockID].append(float(columns[2]))
		Cx[blockID].append(float(columns[3]))
		Cy[blockID].append(float(columns[4]))
X = np.asarray(X)		# Convert to numpy array
X = X.T					# Transpose the matrix
Y = np.asarray(Y)
Y = Y.T
D = np.asarray(D)
D = D.T
Mask = D < rhoscaling*rho0
#D = np.rot90(D, k=2)
D = np.flipud(D)
Cx = np.asarray(Cx)
Cx = Cx.T
Cx[Mask] = None
Cy = np.asarray(Cy)
Cy = Cy.T
Cy[Mask] = None

### PLOT DATA
if DrawDensity:
	im = plt.imshow(D,
					extent=[-xmax, xmax, -ymax, ymax],
					cmap=plt.cm.CMRmap,
					vmin=0,
					vmax=denmaxvalue,
					)
	im.set_interpolation('lanczos')
	plt.colorbar(im, fraction=0.0377, pad=0.02)

if DrawImpurity:
	circle=plt.Circle((ximp, yimp),
				radius=1.5,
				color='#6D8E2B')
	plt.axes().add_patch(circle)

if DrawCirculation:				
	plt.streamplot(X, Y, Cx, Cy,
					color="#76D6FF",
					linewidth=0.9,
					arrowstyle='fancy',
					arrowsize=3,
					density=2.5)

if FirstFrame:
	plt.axes().text(fframex, fframey, FirstFrameTitle,
					color='lightgrey',
					fontsize=10,
					family='Arial')

if IncludePosition:
	X=ximp-xcom
	Y=yimp-ycom
	R=np.sqrt(X*X+Y*Y)
	positionstring = "$r_{I}=" + str(round(R,1)) + "\,\mathrm{\AA}$"
	plt.axes().text(posx, posy, positionstring, 
					color='lightgrey',
					fontsize=11)

if IncludeSpeed:
	X=xlabel[1]
	Y=ylabel[1]
	vxstring = "$v_{" + X + "}=" + str(round(vximp,1)) + "\,\mathrm{ms^{-1}}$"
	vystring = "$v_{" + Y + "}=" + str(round(vyimp,1)) + "\,\mathrm{ms^{-1}}$"
	plt.axes().text(speedvxx, speedvxy, vxstring, 
					color='lightgrey',
					fontsize=11)
	plt.axes().text(speedvyx, speedvyy, vystring, 
					color='lightgrey',
					fontsize=11)

if IncludeEnergy:
	energystring = "$E_{k}=" + str(round(ekin,1)) + "\,\mathrm{K}$"
	plt.axes().text(energyx, energyy, energystring, 
					color='lightgrey',
					fontsize=11)

if IncludeNParticles:
	particlestring = "$N_{He}=" + str(nparticles) + "\,\\#$"
	plt.axes().text(particlex, particley, particlestring, 
					color='lightgrey',
					fontsize=11)

plt.minorticks_on()
plt.axes().set_aspect('equal')
plt.axes().set_xlim([-xmax,xmax])
plt.axes().set_ylim([-ymax,ymax])
plt.xlabel(xlabel, fontsize=14)
plt.ylabel(ylabel, fontsize=14)
plt.axes().text(timex,timey, timestring, color='white', fontsize=14, family='Arial')
plt.savefig(outfile, dpi=dpiout, facecolor='black', bbox_inches='tight')