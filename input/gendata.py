from numpy import *
from scipy import *
from pylab import *
import numpy.matlib as matlib
from shutil import copy
from os import mkdir
import os

def lininc(n,Dx,dx0):
  a=(Dx-n*dx0)*2./n/(n+1)
  dx = dx0+arange(1.,n+1.,1.)*a
  return dx
H = 2000.
N0=5.2e-3
u0 = 0.1


outdir='../runs/QuadPkg'
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

# topo
topo = np.zeros((ny, nx))
# make a wall along south bdy:
topo[0, :] = 0
topo = topo - H

with open(outdir+"/topo.bin", "wb") as f:
	topo.tofile(f)

# dz:
# dz is from the surface down (right?).  Its saved as positive.

dz=zeros(nz)+H/nz

with open(outdir+"/delZvar.bin", "wb") as f:
	dz.tofile(f)

# temperature profile...
g=9.8
alpha = 2e-4
T0 = 28+cumsum(N0**2/g/alpha*(-dz))

with open(outdir+"/TRef.bin", "wb") as f:
	T0.tofile(f)

# save T0 over whole domain
TT0 = matlib.repmat(T0, ny, nx)
with open(outdir+"/T0.bin", "wb") as f:
	TT0.tofile(f)

# save U0 over whole domain
U0 = np.zeros(nz, ny, nx) + u0
with open(outdir+"/U0.bin", "wb") as f:
	U0.tofile(f)

# make drag co-efficients:
qdrag = 0.02 * np.ones((ny, nx))
ldrag = 0.0 * np.ones((ny, nx))
with open(outdir+"/DraguQuad.bin", "wb") as f:
	qdrag.tofile(f)
with open(outdir+"/DragvQuad.bin", "wb") as f:
	qdrag.tofile(f)
with open(outdir+"/DraguLin.bin", "wb") as f:
	ldrag.tofile(f)
with open(outdir+"/DragvLin.bin", "wb") as f:
	ldrag.tofile(f)





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
