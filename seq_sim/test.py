
import numpy as np
import matplotlib.pyplot as plt

fl = 0.2
fd = 0.51
ld = 0.21
fb = 0.29
lb = 0.03

mfed = []

for i in range(100):
    random = np.random.random()
    if random < fl:
        mfed.append(0)
    elif random < fl+fd:
        mfed.append(-np.random.exponential(ld)+1)
    else:
        mfed.append(np.random.exponential(lb)+1)
plt.hist(mfed, bins = np.arange(-0.05,1.5,0.1))