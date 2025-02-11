import numpy as np

f =  'RARR'
posarr = np.loadtxt(f,delimiter=',')

x=posarr[:,0]
y=posarr[:,1]
z=posarr[:,2]

f =  'LOUISPL14_dipoles.dat'
dipolearr = np.loadtxt(f,delimiter=',',skiprows=1)
u=dipolearr[:,0]
v=dipolearr[:,1]
w=dipolearr[:,2]
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(projection='3d')


ax.scatter(x,y,z)
plt.show()


fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(projection='3d')


ax.quiver(x,y,z,u,v,w,length=0.5, normalize=True,pivot='middle')
plt.show()