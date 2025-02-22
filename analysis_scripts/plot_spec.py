import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
f = sys.argv[1]
print(f)
pl_file = f + '_pl.csv'
abs_file = f + '_abs.csv'
# with open(pl_file) as file:
#     lines = [line.rstrip() for line in file]
# exciton_energy = float(lines[1].split()[-1])
# hw = float(lines[2].split()[-1])
# max_vibs = int(lines[3].split()[-1])
df= pd.read_csv(pl_file,skiprows=3)
df = df.rename(columns=lambda x: x.strip())
# fig,[[ax1,ax2],[ax3,ax4]] = plt.subplots(ncols=2,nrows=2,figsize=(10,10))
# energies = []
# for i in range(0,max_vibs+1):
#     energies.append(exciton_energy-i*hw)
# ax1.plot(df['Energy'],df['PLX'])
# for i,e in enumerate(energies):
#     indx = min(range(len(df['Energy'])), key=lambda i: abs(df['Energy'][i]-e))
#     height = df['PLX'][indx]
#     ax1.vlines(e,ymax=height,ymin=0, ls='dashed')
#     ax1.text(e,height+(max(df['PLX'])/40),f'0-{i}',horizontalalignment='center', verticalalignment='center')
# ax1.set_xlabel('Wavenumber (cm$^{-1}$)')
# ax1.set_ylabel('PL Intensity (a.u.)')
# ax2.plot(df['Energy'],df['PLY'])
# for i,e in enumerate(energies):
#     indx = min(range(len(df['Energy'])), key=lambda i: abs(df['Energy'][i]-e))
#     height = df['PLY'][indx]
#     ax2.vlines(e,ymax=height,ymin=0, ls='dashed')
#     ax2.text(e,height+(max(df['PLY'])/40),f'0-{i}',horizontalalignment='center', verticalalignment='center')
# ax2.set_xlabel('Wavenumber (cm$^{-1}$)')
# ax2.set_ylabel('PL Intensity (a.u.)')


dfabs= pd.read_csv(abs_file,skiprows=3)
dfabs = dfabs.rename(columns=lambda x: x.strip())

# ax3.plot(dfabs['Energy'],dfabs['ABSX'])
# ax4.plot(dfabs['Energy'],dfabs['ABSY'])

# plt.tight_layout()

# plt.show()


fig,ax = plt.subplots()

df['Wavelength'] = (1/df['Energy'])*1E7
dfabs['Wavelength'] = (1/dfabs['Energy'])*1E7

ax.plot(dfabs['Wavelength'],dfabs['ABSX'],color='blue',label='x-polarized')
ax.plot(dfabs['Wavelength'],dfabs['ABSY'],color='blue',linestyle='dashed',label='y-polarized')
ax.plot(dfabs['Wavelength'],dfabs['ABSZ'],color='blue',linestyle='dotted',label='z-polarized')

ax2=ax.twinx()
ax2.plot(df['Wavelength'],df['PLX'],color='red')
ax2.plot(df['Wavelength'],df['PLY'],color='red',linestyle='dashed')
ax2.plot(df['Wavelength'],df['PLZ'],color='red',linestyle='dotted')

ax.set_xlim(min(dfabs['Wavelength']),dfabs['Wavelength'].max()+(dfabs['Wavelength'].max()-dfabs['Wavelength'].min()))
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('Absorption')
ax2.set_ylabel('Photoluminescence')
fig.legend(loc='upper right')
ax.set_xlim(200,600)
plt.tight_layout()
plt.savefig(f+'_plabs.png')
plt.show()



fig,ax = plt.subplots()

PLAVG = np.sqrt(df['PLX']**2+df['PLY']**2+df['PLZ']**2)
ax.plot(df['Wavelength'],PLAVG,color='red')
# ax.plot(df['Wavelength'],df['PLZ'],color='blue')
ax.set_xlim(min(dfabs['Wavelength']),dfabs['Wavelength'].max()+(dfabs['Wavelength'].max()-dfabs['Wavelength'].min()))
plt.tight_layout()
plt.show()





import numpy as np
dipolef = f + '_dipoles.dat'
posf = f+'_positions.dat'
posarr = np.loadtxt(posf,delimiter=',',skiprows=1)

x=posarr[:,0]
y=posarr[:,1]
z=posarr[:,2]

dipolearr = np.loadtxt(dipolef,delimiter=',',skiprows=1)
u=dipolearr[:,0]
v=dipolearr[:,1]
w=dipolearr[:,2]

import matplotlib.pyplot as plt

# fig = plt.figure(figsize=(12, 12))
# ax = fig.add_subplot(projection='3d')


# ax.scatter(x,y,z)
# plt.show()

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(projection='3d')


ax.quiver(x,y,z,u,v,w,length=0.5, normalize=False,pivot='middle')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
limits = 10
ax.set_xlim(-limits,limits)
ax.set_ylim(-limits,limits)
ax.set_zlim(-limits,limits)
plt.tight_layout()

plt.savefig(f+'_dipoles.png')
plt.show()



# dipolef = f + '_dipoles.dat'
# dipoles = np.loadtxt(dipolef,delimiter=',',skiprows=1)
# with open(dipolef) as file:
#     lines = [line.rstrip() for line in file]
# latticex=int(lines[0].split()[-1].split(',')[0])
# latticey=int(lines[0].split()[-1].split(',')[1])
# latticez=int(lines[0].split()[-1].split(',')[2])
# lattice_index_arr = np.zeros((latticex,latticey,latticez))
# counter = 0
# for i in range(0,latticex):
#     for j in range(0,latticey):
#         for k in range(0,latticez):
#             lattice_index_arr[i,j,k]=counter
#             counter+=1
# x = []
# y = []
# z = []
# u= []
# v = []
# w = []
# for i in range(0,latticex):
#     for j in range(0,latticey):
#         for k in range(0,latticez):
#             x.append(i)
#             y.append(j)
#             z.append(k)
#             u.append(dipoles[int(lattice_index_arr[i,j,k])][0])
#             v.append(dipoles[int(lattice_index_arr[i,j,k])][1])
#             w.append(dipoles[int(lattice_index_arr[i,j,k])][2])
# ax = plt.figure().add_subplot(projection='3d')
# ax.quiver(x, y, z, u, v, w, length=0.5, normalize=True,pivot='middle')
# ax.set_xlim(-1,max(x)+1)
# ax.set_ylim(-1,max(y)+1)
# ax.set_zlim(-1,max(z)+1)
# plt.savefig(f+'_dipoles.png')

# plt.show()

cd_file = f + '_cd.csv'


df= pd.read_csv(cd_file,skiprows=3)
df = df.rename(columns=lambda x: x.strip())
print(df.columns)
fig,ax = plt.subplots()

df['Wavelength'] = (1/df['Energy'])*1E7
# ax.plot(df['Energy'],df['CD'])
# ax.set_xlabel('Wavenumber (cm$^{-1}$)')
ax.plot(df['Wavelength'],df['CD'])
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('CD')
plt.savefig(f+'_cd.png')
plt.show()

cpl_file = f + '_cpl.csv'


df= pd.read_csv(cpl_file,skiprows=4)
df = df.rename(columns=lambda x: x.strip())
fig,ax = plt.subplots()

df['Wavelength'] = (1/df['Energy'])*1E7
ax.plot(df['Energy'],df['GLUM'])
ax.set_xlabel('Wavenumber (cm$^{-1}$)')
ax.set_ylabel('g-factor (arbitrary scale)')
plt.savefig(f+'_cpl_cm^-1.png')

plt.show()