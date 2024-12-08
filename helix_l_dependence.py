# Test length of helix
def loadInputFile(fname):
    with open(fname) as file:
        lines = [line.rstrip() for line in file]
    input_vars = {}
    for line in lines:
        if line.split(' ')[0] == '#':
            continue
        input_vars[line.split(' ')[0]] = line.split(' ')[1]
    return input_vars
def createInputFile(varsdict,fname):
    with open(fname,mode='w') as f:
        for v in varsdict.keys():
            f.write(f'{v} {varsdict[v]}\n')
import subprocess
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation



input_vars = loadInputFile('EXAMPLE.inp')
for i in range(1, 25):
    input_vars['LATTICE_DIMZ'] = i
    input_vars['NUMBER_CONFIGURATIONS'] = 100
    input_vars['INPUT_NAME'] = f'LOUISPL{i}'
    createInputFile(input_vars,'EXAMPLE.inp')
    subprocess.run(["powershell",f"python chiralpl.py EXAMPLE.inp"]) 

nframes = 24
plt.subplots_adjust(top=1, bottom=0, left=0, right=1)
def animate(i):
    i = i +1
    im = plt.imread(f'simulations/LOUISPL{i}/'+'LOUISPL'+str(i)+'_cd.png')
    plt.imshow(im)

anim = FuncAnimation(plt.gcf(), animate, frames=nframes)
anim.save('output.gif', writer='imagemagick')
