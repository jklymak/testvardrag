from numpy import *
from scipy import *
from pylab import *
import numpy.matlib as matlib
from shutil import copy
from os import mkdir
import os
import matplotlib.pyplot as plt

def lininc(n,Dx,dx0):
  a=(Dx-n*dx0)*2./n/(n+1)
  dx = dx0+arange(1.,n+1.,1.)*a
  return dx
H = 200.
N0=5.2e-3
u0 = 0.01

outdir='../runs/quadPkg/'
try:
  mkdir(outdir)
except:
  print( outdir+' Exists')
try:
  mkdir(outdir+'/figs')
except:
  print(outdir+'/figs Exists')
copy('./gendata.py',outdir)

# These must match ../code/SIZE.h
ny = 2*40
nx = 2*40
nz = 25

dx = np.ones(nx) * 2000
dy = np.ones(ny) * 2000

with open(outdir+"/delXvar.bin", "wb") as f:
	dx.tofile(f)
with open(outdir+"/delYvar.bin", "wb") as f:
	dy.tofile(f)

x = np.cumsum(dx)
x = x - np.mean(x)
y = np.cumsum(dy)
y = y - np.mean(y)

# topo
topo = np.zeros((ny, nx))
# make a wall along south bdy:
topo = topo - H
topo[0, :] = 0

with open(outdir+"/topo.bin", "wb") as f:
	topo.tofile(f)

# dz:
# dz is from the surface down (right?).  Its saved as positive.

dz=zeros(nz)+H/nz

with open(outdir+"/delZvar.bin", "wb") as f:
	dz.tofile(f)
z = -np.cumsum(dz)

# temperature profile...
g=9.8
alpha = 2e-4
T0 = 28+cumsum(N0**2/g/alpha*(-dz))

with open(outdir+"/TRef.bin", "wb") as f:
	T0.tofile(f)

# save T0 over whole domain
TT0 = np.zeros((nz, ny, nx)) + T0[:, np.newaxis, np.newaxis]
with open(outdir+"/T0.bin", "wb") as f:
	TT0.tofile(f)

# save U0 over whole domain
U0 = np.zeros((nz, ny, nx)) + u0
with open(outdir+"/U0.bin", "wb") as f:
	U0.tofile(f)

# make drag co-efficients:
qdrag = 1.0e-2 * np.ones((ny, nx))
ldrag = 0.0 * np.ones((ny, nx))

#X, Y = np.meshgrid(x, y)
#R2 = X**2 + Y**2
#ldrag[R2<10000**2] = 0.005
fig, ax = plt.subplots()
pc = ax.pcolormesh(ldrag)
fig.colorbar(pc)
with open(outdir+"/DraguQuad.bin", "wb") as f:
	qdrag.tofile(f)
with open(outdir+"/DragvQuad.bin", "wb") as f:
	qdrag.tofile(f)
with open(outdir+"/DraguLin.bin", "wb") as f:
	ldrag.tofile(f)
with open(outdir+"/DragvLin.bin", "wb") as f:
	ldrag.tofile(f)

# make a file that spreads the stress...  Note sum fac*dz = 1

dragfac = np.zeros((nz, ny, nx))
Ndrag = 1
for i in range(nx):
    for j in range(ny):
        ind = np.where(z>=topo[j, i])[0]
        if len(ind) > Ndrag:
            print(ind)
            dragfac[range(ind[-Ndrag], ind[-1]+1), j, i] = 1
            sum = np.sum(dragfac[:, j, i] * dz[:])
            dragfac[range(ind[-Ndrag], ind[-1]+1), j, i] =  dragfac[ind[-Ndrag]:ind[-1]+1, j, i] / sum

print(dragfac[:, 4, 4])

with open(outdir+"/BotDragFac.bin", "wb") as f:
	dragfac.tofile(f)



## Copy some other files
import shutil
shutil.copy('data', outdir+'/data')
shutil.copy('eedata', outdir)
shutil.copy('data.diagnostics', outdir)
shutil.copy('data.var_bot_drag', outdir)
shutil.copy('data.pkg', outdir+'/data.pkg')
# also store these.  They are small and helpful to document what we did
for nm in {'input','code','build_options','analysis'}:
    to_path = outdir+'/'+nm
    if os.path.exists(to_path):
        shutil.rmtree(to_path)
    shutil.copytree('../'+nm, outdir+'/'+nm)
