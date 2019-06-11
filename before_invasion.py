#所有模型均来自于论文_
import numpy as np
import sympy as sy
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import random


def init(q,death,n,extern=0):
    p=t=m=c2=m2=np.zeros(n)
    c=np.ones(n)
    random.seed()
    c=np.array([q*(1-q)**(n-i-1) for i in range(n)])
    c2=np.array(extern)#外来物种扩散率，一般只考虑只有一个种群
    m=[death for i in range(n)]
    return n,p,t,c,m,c2


def before_invasion(args):
    '''
    以下为未受外来物种入侵的多物种竞争共存的集合种群模型，假设pi与时间无关，但是保留参数p初值
    所有参数均为数组np.array，pi表示物种i占据生境斑块的比例; ci表示物种i的扩散率; mi表示物种i的灭绝率; n表示物种种数。
    直接解各个物种比例稳定值方程组，此处pi为符号数组
    '''
    n,p,t,c,m=args[0:5]
    pi=sy.symbols('p1:%d'%(n+1))
    #将方程组分离为增广矩阵形式，即常数项单独分离开来
    #注意下面下标是错位的，从0开始
    para=[[0 for i in range(n)] for i in range(n)]#常数项矩阵
    for i in range(n):
        for j in range(n):
            if j<i:
                para[i][j]=1+(c[j]/c[i])#int后期略去
            elif i==j:
                para[i][j]=1
    para=sy.Matrix(para)
    para2=[[(1-m[i]/c[i])] for i in range(n)]#系数项矩阵
    para2=sy.Matrix(para2)
    print(para,para2)
    result=sy.linsolve((para,para2),pi)
    result=list(list(result)[0])
    print(result)
    X=np.arange(n)+1
    Y=np.array(result)
    plt.figure(figsize=(n,6))
    plt.bar(X,Y,width=0.30,facecolor = 'lightskyblue',edgecolor = 'white')
    sum=np.array(result).sum()
    print(sum)
    plt.ylim(0,float(sum))
    plt.show()


def cut_in_line(args):
    '''
    以下为外来物种插队竞争模型，保持三点基本假设
    1 外来物种入侵前后所有物种的竞争力和扩散率排序规律不变，值可以变化（为了使竞争力之和趋近于1）
    2 外来物种入侵前后，所有物种的比例稳定值之和恒定
    3 物种进化过程中遵循 “趋利避害”的原则
    参数与上面相同，同时方程形式也与入侵前模型相同，只是物种数为n+1
    #此时外来物种与本地物种并无区别
    '''
    n,p,t,c,m,c2=args
    death=m[0]
    c=list(c)
    c.append(c2)
    c=np.array(c)
    c.sort()
    n+=1
    m=np.array([death for i in range(n)])
    print(c)
    before_invasion((n,p,t,c,m))


def show(plt):
    pass