import numpy as np
import sympy as sy


N=4


def init():
    n=N
    t=m=np.zeros(n)
    c=np.ones(n)
    return t,c,m,n


def before_invasion(args):
    '''
    以上为未受外来物种入侵的多物种竞争共存的集合种群模型
    所有参数均为数组np.array，pi表示物种i占据生境斑块的比例; ci表示物种i的扩散率; mi表示物种i的灭绝率; n表示物种种数。
    直接解各个物种比例稳定值方程组，此处p为符号数组
    '''
    t,c,m,n=args
    pi=sy.symbols('p1:5')
    #将方程组分离为增广矩阵形式，即常数项单独分离开来
    ##求pj从0到i-1的和，是以i为下标的数组，以下矩阵会用到
    #注意下面下标是错位的，从0开始
    para=[[0 for i in range(n)] for i in range(n)]#常数项矩阵
    for i in range(n):
        for j in range(n):
            if j<i:
                para[i][j]=1+int(c[j]/c[i])#int后期略去
            elif i==j:
                para[i][j]=1
    para=sy.Matrix(para)
    para2=[]#系数项矩阵
    for i in range(n):
        para2.append(m[i]/c[i]-1)
    print(para,para2)
