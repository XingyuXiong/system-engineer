import numpy as np
import sympy as sy


def init():
    n=10
    p=t=c=m=np.zeros(n)
    return p,t,c,m,n


def before_invasion(args):
    '''
    以上为未受外来物种入侵的多物种竞争共存的集合种群模型
    所有参数均为数组np.array，pi表示物种i占据生境斑块的比例; ci表示物种i的扩散率; mi表示物种i的灭绝率; n表示物种种数。
    '''
    p,t,c,m,n=args
    derive=sy.diff(p,t,1)
    block=np.array([1-np.sum(p[:i]) for i in range(n)])#剩余的生境斑块
    spread=p*t*block#物种扩散部分
    instinct=m*p#灭绝
    unknown=np.array([np.sum(c*p[:i])*p for i in range(n)])#这部分意思没看懂
    return -derive+spread-instinct-unknown