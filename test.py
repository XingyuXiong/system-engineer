import numpy
import matplotlib.pyplot as plt

#x = numpy.linspace(0, 2, 100)
x=numpy.array([1,2,3,4])
#plt.figure(figsize=(8, 4))

for i in range(3):
    plt.plot(x, x**i, label='%d'%i)

# 设置展示区间
plt.xlim(-1, 3)
plt.ylim(-1, 10)

# 设置展示信息
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.title('simple_plot')
plt.legend(loc='best')

# 展示图片
plt.show()
# 关闭图片
plt.close()