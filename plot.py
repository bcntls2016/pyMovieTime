#!/usr/bin/python
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
from plotsettings import *
vars = __import__(sys.argv[2].replace('.py',''))

def getStateColor( state ):
    if state == 'gs':
        return ColorGS;
    elif State == 'es1':
        return ColorES1;
    elif State == 'es2':
        return ColorES2;
    elif State == 'ion':
        return ColorIon;
    else:
        return ColorGS;

### SOME GLOBAL SETTINGS
plt.rcParams['text.color'] = 'lightgrey'
plt.rcParams['axes.labelcolor'] = 'lightgrey'
plt.rcParams['axes.facecolor'] = "black"
plt.rcParams['axes.edgecolor'] = '#A9A9A9'
plt.rcParams['grid.color'] = '#A9A9A9'
plt.rcParams['xtick.color'] = '#A9A9A9'
plt.rcParams['ytick.color'] = '#A9A9A9'
dpiout = 300

### LABEL PLACEMENT
timex = -vars.xmax + 2.5
timey = vars.ymax - 4.5
fframex = -vars.xmax + 2.5
fframey = -vars.ymax + 2.5

posx = vars.xmax - 22.5
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
ID = sys.argv[1].split('.')[1]
outfile = "den." + ID + ".png"
timestring = "$t=" + str(vars.time) + "\,\mathrm{ps}$"

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
			extent=[-vars.xmax, vars.xmax, -vars.ymax, vars.ymax],
			cmap=plt.cm.CMRmap,
			vmin=0,
			vmax=denmaxvalue)
	im.set_interpolation('lanczos')
	cbar=plt.colorbar(im, ticks=[0, 0.0109, 0.0218, 0.0327], pad=0.02)
	cbar.ax.set_yticklabels([r'$0$', r'$\frac{1}{2}\rho_0$', r'$\rho_0$', r'$\frac{3}{2}\rho_0$']) 
if DrawImpurity:
        circle=plt.Circle((vars.ximp, vars.yimp),
			radius=1.25,
                        color=getStateColor(State))
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
	X=vars.ximp-vars.xcom
	Y=vars.yimp-vars.ycom
	R=np.sqrt(X*X+Y*Y)
	positionstring = "$r_{I}=" + str(round(R,1)) + "\,\mathrm{\AA}$"
	plt.axes().text(posx, posy, positionstring,
			color='lightgrey',
			fontsize=11)

if IncludeSpeed:
	X=xlabel[1]
	Y=ylabel[1]
	vxstring = "$v_{" + X + "}=" + str(round(vars.vximp,1)) + "\,\mathrm{ms^{-1}}$"
	vystring = "$v_{" + Y + "}=" + str(round(vars.vyimp,1)) + "\,\mathrm{ms^{-1}}$"
	plt.axes().text(speedvxx, speedvxy, vxstring,
			color='lightgrey',
			fontsize=11)
	plt.axes().text(speedvyx, speedvyy, vystring,
			color='lightgrey',
			fontsize=11)

if IncludeEnergy:
	energystring = "$E_{k}=" + str(round(vars.ekin,1)) + "\,\mathrm{K}$"
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
plt.axes().set_xlim([-vars.xmax,vars.xmax])
plt.axes().set_ylim([-vars.ymax,vars.ymax])
plt.xlabel(xlabel, fontsize=14)
plt.ylabel(ylabel, fontsize=14)
plt.axes().text(timex,timey, timestring, color='white', fontsize=14, family='Arial')
plt.savefig(outfile, dpi=dpiout, facecolor='black', bbox_inches='tight')
