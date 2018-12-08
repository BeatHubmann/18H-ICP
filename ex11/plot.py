import mpl_toolkits.mplot3d
from numpy import *
from pylab import *

u = loadtxt("x_cholesky.txt")
N = int(sqrt(u.shape[0]))
print("N=%d" % N)
u.shape = (N,N)

fig = figure()
ax = fig.gca(projection="3d")
x,y = mgrid[0:1:N*1j, 0:1:N*1j]
ax.plot_surface(x,y,u)
show()

pcolormesh(x,y,u)
colorbar()
show()

u = loadtxt("x_jacobi.txt")
N = int(sqrt(u.shape[0]))
print("N=%d" % N)
u.shape = (N,N)

fig = figure()
ax = fig.gca(projection="3d")
x,y = mgrid[0:1:N*1j, 0:1:N*1j]
ax.plot_surface(x,y,u)
show()

pcolormesh(x,y,u)
colorbar()
show()

u = loadtxt("x_gauss_seidel.txt")
N = int(sqrt(u.shape[0]))
print("N=%d" % N)
u.shape = (N,N)

fig = figure()
ax = fig.gca(projection="3d")
x,y = mgrid[0:1:N*1j, 0:1:N*1j]
ax.plot_surface(x,y,u)
show()

pcolormesh(x,y,u)
colorbar()
show()

u = loadtxt("x_diff_chol_jaco.txt")
N = int(sqrt(u.shape[0]))
print("N=%d" % N)
u.shape = (N,N)

fig = figure()
ax = fig.gca(projection="3d")
x,y = mgrid[0:1:N*1j, 0:1:N*1j]
ax.plot_surface(x,y,u)
show()

pcolormesh(x,y,u)
colorbar()
show()

u = loadtxt("x_diff_chol_gaus.txt")
N = int(sqrt(u.shape[0]))
print("N=%d" % N)
u.shape = (N,N)

fig = figure()
ax = fig.gca(projection="3d")
x,y = mgrid[0:1:N*1j, 0:1:N*1j]
ax.plot_surface(x,y,u)
show()

pcolormesh(x,y,u)
colorbar()
show()