#所有模型均来自于论文_
import numpy as np
import sympy as sy
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import random


TIME_LIMIT=20
FLAG=0
DATA=[]


def init(p,q,death,n,extern=0):
    t=0
    p=np.array(p)
    m=c2=m2=np.zeros(n)
    c=np.ones(n)
    random.seed()
    c=np.array([q*(1-q)**(n-i-1) for i in range(n)])
    c2=extern#外来物种扩散率，一般只考虑只有一个种群
    m=[death for i in range(n)]
    return n,p,q,t,c,m,c2

#before_invasion和cut_in_line都是用来解稳定物种比例方程组，但是好像不符合change函数给出的预期
#即论文给出的稳定比例表达式不是微分方程组真正的稳定解
def before_invasion(args):
    '''
    以下为未受外来物种入侵的多物种竞争共存的集合种群模型，假设pi与时间无关，但是保留参数p初值
    所有参数均为数组np.array，pi表示物种i占据生境斑块的比例; ci表示物种i的扩散率; mi表示物种i的灭绝率; n表示物种种数。
    直接解各个物种比例稳定值方程组，此处pi为符号数组
    '''
    n,p,q,t,c,m=args[0:5]
    pi=sy.symbols('p1:%d'%(n+1))
    #将方程组分离为增广矩阵形式，即常数项单独分离开来
    #注意下面下标是错位的，从0开始
    para=[[0 for i in range(n)] for i in range(n)]#常数项矩阵
    for i in range(n):
        for j in range(n):
            if j<i:
                para[i][j]=1+(c[j]/c[i])
            elif i==j:
                para[i][j]=1
    para=sy.Matrix(para)
    para2=[[(1-m[i]/c[i])] for i in range(n)]#系数项矩阵
    para2=sy.Matrix(para2)
    result=sy.linsolve((para,para2),pi)
    result=list(list(result)[0])
    #print(result)
    X=np.arange(n)+1
    Y=np.array(result)
    plt.figure(figsize=(n,6))
    plt.bar(X,Y,width=0.30,facecolor = 'lightskyblue',edgecolor = 'white')
    sum=np.array(result).sum()
    #print(sum)
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
    n,p,q,t,c,m,c2=args
    death=m[0]
    c=list(c)
    c.append(c2)
    c=np.array(c)
    c.sort()
    n+=1
    m=np.array([death for i in range(n)])
    #print(c)
    before_invasion((n,p,t,c,m))


def change(args,flag):
    '''
    三种模型，未受外来物种入侵，插队竞争，等位竞争，flag为0，1，2
    '''
    FLAG=flag
    i=0
    pex=0
    args=args
    while(i<TIME_LIMIT):
        args,pex=cal(args,pex,flag)
        i+=1
    show(args)


def cal(args,pex,flag):
    n,p,q,t,c,m,c2=args
    pnext=[0 for i in range(n)]
    if flag==0:
        for i in range(n):
            pnext[i]=p[i]+c[i]*p[i]*(1-p[:i].sum())-m[i]*p[i]-(c[:(i-1)]*p[:(i-1)]).sum()*p[i]
    if flag==1:
        for i in range(n):
            pnext[i]=p[i]+c[i]*p[i]*(1-p[:i].sum())-m[i]*p[i]-(c[:(i-1)]*p[:(i-1)]).sum()*p[i]
            pexnext=pex+c[i]*pex*(1-p[:i].sum())-m[i]*pex-(c[:(i-1)]*p[:(i-1)]).sum()*pex
    if flag==2:
        #外来物种等位竞争共存模型，外来物种专一与本地物种竞争实现共存
        #首先找出外来物种与之竞争的本地物种k，这里假设k是竞争力最接近外来物种但小于它的本地物种
        k=n-1
        for i in range(n):
            if c[i]>c2:
                k=max(0,i-1)
                break
        for i in range(n):
            if i<k:
                pnext[i]=p[i]+c[i]*p[i]*(1-p[:i].sum())-m[i]*p[i]-(c[:(i-1)]*p[:(i-1)]).sum()*p[i]
            elif i==k:
                pnext[i]=p[i]+c[i]*p[i]*(1-p[:i].sum()-pex)-m[i]*p[i]-(c[:(i-1)]*p[:(i-1)]).sum()*p[i]
                pexnext=pex+c[i]*pex*(1-p[:i].sum()-pex)-m[i]*pex-(c[:(i-1)]*p[:(i-1)]).sum()*pex
            else:
                pnext[i]=p[i]+c[i]*p[i]*(1-p[:i].sum()-pex)-m[i]*p[i]-(c[:(i-1)]*p[:(i-1)]).sum()*p[i]-c[k]*pex*p[i]
    #检查种群分布是否有小于0的情况，如果有，变成0
    for i in range(n):
        if pnext[i]<0:
            pnext[i]=0
    save(p,pex)
    p=np.array(pnext)
    t+=1
    args=n,p,q,t,c,m,c2
    return args,pex


def save(p,pex):
    datarow=list(p)
    datarow.append(pex)
    DATA.append(datarow)


def show(args):
    n,p,q,t,c,m,c2=args
    print(c)
    x=np.array([i for i in range(len(DATA))])
    for i in range(len(DATA[0])):
        y=np.array([j[i] for j in DATA])
        if i<len(DATA[0])-1:
            plt.plot(x,y,label='%s:%s'%(i+1,c[i]))
        elif FLAG:
            #如果FLAG不为0说明该模型含有外来物种
            plt.plot(x,y,label='入侵种群')
    plt.legend(loc='best')#用来显示不同标签
    plt.savefig(fname='change_flag%s_q%s_m%s.png'%(FLAG,q,m[0]), dpi=500)
    plt.show()
    plt.close()