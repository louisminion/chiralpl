import numpy as np
import matplotlib.pyplot as plt
import sys
f = sys.argv[1]
print(f)
evals = np.loadtxt(f)

fig,ax = plt.subplots()

ax.stem(evals)
ax.set_ylabel('Eigenvalue')
ax.set_xlabel('Number')
plt.show()