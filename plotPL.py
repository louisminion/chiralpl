import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
f = sys.argv[1]
print(f)
with open(f) as file:
    lines = [line.rstrip() for line in file]
exciton_energy = float(lines[1].split()[-1])
hw = float(lines[2].split()[-1])
max_vibs = int(lines[3].split()[-1])
df= pd.read_csv(f,skiprows=4)
df = df.rename(columns=lambda x: x.strip())
fig,[ax,ax1,ax2] = plt.subplots(ncols=3,figsize=(15,5))
avgPL = (df['PLX']+df['PLY'])/2
ax.plot(df['Energy'],avgPL)
# ax.set_ylabel('Eigenvalue')
# sticks = [exciton_energy, exciton_energy-1, exciton_energy-2]
energies = []
for i in range(0,max_vibs+1):
    energies.append(exciton_energy-i*hw)
# energies = [exciton_energy, exciton_energy-hw, exciton_energy-2*hw, exciton_energy-3*hw]
for i,e in enumerate(energies):
    indx = min(range(len(df['Energy'])), key=lambda i: abs(df['Energy'][i]-e))
    height = avgPL[indx]
    ax.vlines(e,ymax=height,ymin=0, ls='dashed')
    ax.text(e,height+(max(avgPL)/40),f'0-{i}',horizontalalignment='center', verticalalignment='center')
# ax.stem(sticks, [1,1,1])
# ax.set_xlabel('Number')
ax.set_xlabel('Wavenumber (cm$^{-1}$)')
ax.set_ylabel('PL Intensity (a.u.)')


ax1.plot(df['Energy'],df['PLX'])

for i,e in enumerate(energies):
    indx = min(range(len(df['Energy'])), key=lambda i: abs(df['Energy'][i]-e))
    height = df['PLX'][indx]
    ax1.vlines(e,ymax=height,ymin=0, ls='dashed')
    ax1.text(e,height+(max(df['PLX'])/40),f'0-{i}',horizontalalignment='center', verticalalignment='center')
# ax.stem(sticks, [1,1,1])
# ax.set_xlabel('Number')
ax1.set_xlabel('Wavenumber (cm$^{-1}$)')
ax1.set_ylabel('PL Intensity (a.u.)')


ax2.plot(df['Energy'],df['PLY'])

for i,e in enumerate(energies):
    indx = min(range(len(df['Energy'])), key=lambda i: abs(df['Energy'][i]-e))
    height = df['PLY'][indx]
    ax2.vlines(e,ymax=height,ymin=0, ls='dashed')
    ax2.text(e,height+(max(df['PLY'])/40),f'0-{i}',horizontalalignment='center', verticalalignment='center')
# ax.stem(sticks, [1,1,1])
# ax.set_xlabel('Number')
ax2.set_xlabel('Wavenumber (cm$^{-1}$)')
ax2.set_ylabel('PL Intensity (a.u.)')

plt.tight_layout()


plt.show()