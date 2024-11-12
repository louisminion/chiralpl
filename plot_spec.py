import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
f = sys.argv[1]
print(f)
pl_file = f + '_pl.csv'
abs_file = f + '_abs.csv'
with open(pl_file) as file:
    lines = [line.rstrip() for line in file]
exciton_energy = float(lines[1].split()[-1])
hw = float(lines[2].split()[-1])
max_vibs = int(lines[3].split()[-1])
df= pd.read_csv(pl_file,skiprows=4)
df = df.rename(columns=lambda x: x.strip())
fig,[[ax1,ax2],[ax3,ax4]] = plt.subplots(ncols=2,nrows=2,figsize=(10,10))
energies = []
for i in range(0,max_vibs+1):
    energies.append(exciton_energy-i*hw)
ax1.plot(df['Energy'],df['PLX'])
for i,e in enumerate(energies):
    indx = min(range(len(df['Energy'])), key=lambda i: abs(df['Energy'][i]-e))
    height = df['PLX'][indx]
    ax1.vlines(e,ymax=height,ymin=0, ls='dashed')
    ax1.text(e,height+(max(df['PLX'])/40),f'0-{i}',horizontalalignment='center', verticalalignment='center')
ax1.set_xlabel('Wavenumber (cm$^{-1}$)')
ax1.set_ylabel('PL Intensity (a.u.)')
ax2.plot(df['Energy'],df['PLY'])
for i,e in enumerate(energies):
    indx = min(range(len(df['Energy'])), key=lambda i: abs(df['Energy'][i]-e))
    height = df['PLY'][indx]
    ax2.vlines(e,ymax=height,ymin=0, ls='dashed')
    ax2.text(e,height+(max(df['PLY'])/40),f'0-{i}',horizontalalignment='center', verticalalignment='center')
ax2.set_xlabel('Wavenumber (cm$^{-1}$)')
ax2.set_ylabel('PL Intensity (a.u.)')


dfabs= pd.read_csv(abs_file,skiprows=3)
dfabs = dfabs.rename(columns=lambda x: x.strip())

ax3.plot(dfabs['Energy'],dfabs['ABSX'])
ax4.plot(dfabs['Energy'],dfabs['ABSY'])

plt.tight_layout()

plt.show()


fig,ax = plt.subplots()

df['Wavelength'] = (1/df['Energy'])*1E7
dfabs['Wavelength'] = (1/dfabs['Energy'])*1E7

ax.plot(dfabs['Wavelength'],dfabs['ABSX'],color='blue',label='x-polarized')
ax.plot(dfabs['Wavelength'],dfabs['ABSY'],color='blue',linestyle='dashed',label='y-polarized')

ax2=ax.twinx()
ax2.plot(df['Wavelength'],df['PLX'],color='red')
ax2.plot(df['Wavelength'],df['PLY'],color='red',linestyle='dashed')

ax.set_xlim(min(dfabs['Wavelength']),dfabs['Wavelength'].max()+(dfabs['Wavelength'].max()-dfabs['Wavelength'].min()))
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('Absorption')
ax2.set_ylabel('Photoluminescence')
fig.legend(loc='upper right')
plt.tight_layout()
plt.show()