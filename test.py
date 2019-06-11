from sympy import *
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import numpy as np

n=10
plt.figure(figsize=(9,6))
n=8
X = np.arange(n)+1
#X是1,2,3,4,5,6,7,8,柱的个数
# numpy.random.uniform(low=0.0, high=1.0, size=None), normal
#uniform均匀分布的随机数，normal是正态分布的随机数，0.5-1均匀分布的数，一共有n个
Y1 = np.random.uniform(0.5,1.0,n)
Y2 = np.random.uniform(0.5,1.0,n)
plt.bar(X,Y1,width = 0.35,facecolor = 'lightskyblue',edgecolor = 'white')
plt.ylim(0,+1.25)
plt.show()