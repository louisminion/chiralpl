import subprocess
import os
import shutil
import sys
import argparse



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


parser = argparse.ArgumentParser(
                    prog='CHIRALPL',
                    description='What the program does',
                    epilog='Text at the bottom of help')
parser.add_argument("input_file", help="Input file specifying simulation parameters.")
parser.add_argument('--version', action='version', version='%(prog)s 0.1')
# parser.add_argument('-v','--verbose', action='store_true')
# parser.add_argument('-n','--no_plot', action='store_true')
args = parser.parse_args()
with open(args.input_file) as f:
    lines = [line.rstrip() for line in f]

for l in lines:
    if lines[0].split()[0] == 'INPUT_NAME':
        name = lines[0].split()[1]

print(name)


input_vars = loadInputFile(args.input_file)


# if not os.path.exists('../simulations'):
#     os.mkdir('../simulations')

# if not os.path.exists(f'../simulations/{name}'):
#     os.mkdir(f'../simulations/{name}')

# shutil.copy(args.input_file, f'../simulations/{name}/{name}.inp')
# subprocess.run(["powershell",f".\/bin/chiralpl.exe {args.input_file} | tee {name}.out"]) 
# if args.verbose:
#     print('Cleaning up files')
# shutil.move(f"{name}.out",f"../simulations/{name}/{name}.out")
# shutil.move(f"{name}_dipoles.dat",f"../simulations/{name}/{name}_dipoles.dat")
# shutil.move(f"{name}_abs.csv",f"../simulations/{name}/{name}_abs.csv")
# shutil.move(f"{name}_pl.csv",f"../simulations/{name}/{name}_pl.csv")
# shutil.move(f"{name}_cpl.csv",f"../simulations/{name}/{name}_cpl.csv")


# os.chdir(f'simulations/{name}')

# def runSimulation(threads, input_vars, name,)
import time
times = []
no_configs = []
for i in range(1,100):
    input_vars['NUMBER_CONFIGURATIONS'] = i
    createInputFile(input_vars,args.input_file)
    stime = time.time_ns()
    subprocess.run(["powershell",f"../bin/chiralpl.exe {args.input_file} | tee {name}.out"]) 
    ftime = time.time_ns()
    times.append((ftime-stime)/(1E9))
    no_configs.append(i)
print('Cleanup')
print('Plotting')
os.remove(f'{name}_cd.csv')
os.remove(f'{name}_pl.csv')
os.remove(f'{name}_abs.csv')
os.remove(f'{name}_dipoles.dat')
os.remove(f'{name}.out')
import matplotlib.pyplot as plt

fig,ax = plt.subplots()

ax.scatter(no_configs,times)
ax.set_ylabel('Time (s)')
ax.set_xlabel('Number of configurations')

fig.savefig(f'{name}_timing_over_configs.png')
plt.show()
