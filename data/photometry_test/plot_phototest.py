import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d


data = np.load('photo_test_result.npy').transpose()
data_abs = np.load('photo_test_result_abs.npy').transpose()
# print(data.shape)
# newdata = []
# for i in range(data.shape[1]):
#     if str(i)[-1] == '0':
#         newdata.append(data[:, i])
# newdata = np.array(newdata)
plt.plot(data[0], data[2] / data_abs[2])
plt.savefig('phototest.pdf')
plt.clf()
plt.plot(data[0], data[1])
f = interp1d(data[1], data[0])
f2 = interp1d(data[0], data_abs[2]/data[2])
print(f2(f(2)))
plt.savefig('phototest1.pdf')
plt.clf()