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
print(u)
print(v)
print(w)


print(np.sqrt(u[0]**2+v[0]**2+w[0]**2), 'MAGNITUDE')
# u = [1,1,1,1,1,1]
# v = [0,0,0,0,0,0]
# w = [0,0,0,0,0,0]
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(projection='3d')


ax.scatter(x,y,z)
plt.show()


fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(projection='3d')


ax.quiver(x,y,z,u,v,w,length=0.5, normalize=False,pivot='middle')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
limits = 6
ax.set_xlim(-limits,limits)
ax.set_ylim(-limits,limits)
ax.set_zlim(-limits,limits)


plt.show()
plt.tight_layout()
# print(u)
# print(z)