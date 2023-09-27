import numpy as np
import matplotlib.pyplot as plt


x = np.linspace(0, 1000, 1000)
y = x * 5


fig, ax = plt.subplots()
ax.plot(x, y)
plt.show()
