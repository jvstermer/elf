import numpy as np
from elf import line_models as lm
import matplotlib.pyplot as plt

x = np.arange(0, 20, .5)
a = 1
b = 2

y = lm.polynomial(a,b, wave= x)
plt.plot(x, a +b*x,label='bx+a')
plt.plot(x,y, '.',label = 'pol')
plt.legend()
plt.show()

