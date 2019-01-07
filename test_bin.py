import numpy as np
import matplotlib.pyplot as plt

n = np.asarray([80, 4, 0, 2, 3, 4])
bins = np.arange(0, 121, 20)

x = np.random.randint(0, 20, 80)
y = np.random.randint(20,40, 4)
z = np.random.randint(60,80, 2)
zz = np.random.randint(80,100, 3)
xx = np.random.randint(100, 121, 4)
x = np.concatenate((x,y,z,zz,xx))

lim = x.min() #np.mean(vel_diff)
mask = x >= lim
while len(np.where(mask == True)[0]) > len(x) * .1  and lim < x.max():
    
    lim += 1
    mask = x >= lim
#plt.hist(x)
#plt.hist(x, range = [x.min(), lim])
plt.hist(np.clip(x, x.min(), lim) , alpha = .6)
plt.show()
"""
mask = n > np.sum(n) * 5/100
new_bin = np.append( bins[0],( list(bins[1:][mask]), [bins[1:][~mask].min(), bins.max()] ))

#new_bin = [bins[0], bins[1:][mask], bins[1:][~mask].min(), bins.max()]

print(new_bin)
plt.hist(x, bins)
plt.hist(x, new_bin, alpha = .6)
plt.show()"""
