import numpy as np
import matplotlib.pyplot as plt 
from py.elf import line_models as lm

coeff = (1,1,1)
X = np.arange(-10,10.1)
p = lm.Polynomial(*[p for p in coeff])
v = lm.polynomial(*[p for p in coeff], wave = X)

plt.plot(X,p(X))
plt.plot(X,v,'.r')
plt.show()

