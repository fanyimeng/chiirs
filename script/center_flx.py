import numpy as np
import matplotlib.pyplot as plt

x = np.arange(0, 3200)
y = np.arange(0, 3200)
arr = np.ones((y.size, x.size))

cx = 1200.
cy = 1600.
r = 50.

# The two lines below could be merged, but I stored the mask
# for code clarity.
mask = (x[np.newaxis,:]-cx)**2 + (y[:,np.newaxis]-cy)**2 < r**2
print(arr[mask].sum())
arr[mask] = 2

# This plot shows that only within the circle the value is set to 123.
plt.figure(figsize=(6, 6))
plt.pcolormesh(x, y, arr)
plt.colorbar()
plt.show()