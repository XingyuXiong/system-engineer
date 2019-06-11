import sympy as sy
import before_invasion as bi


if __name__=='__main__':
    '''
    以下为未受外来物种入侵的多种群竞争共存模型，有三个参数
    第一个参数是竞争力分布的参数q，假设竞争力为从小到大的等比数列q*(1-q)**(n-i),i为第i个种群；
    第二个参数是灭绝率death或者说m；
    第三个参数是种群数n
    绘图结果为各个种群的稳定生境斑块比例柱状分布图
    请谨慎调参，如果出现稳定生境斑块比例为负数的情况将无法绘图，结果我已打印在控制台中
    '''
    #bi.before_invasion(bi.init(0.5,0.04,4))
    bi.cut_in_line(bi.init(0.5,0.04,4,0.7))
    #后续模型待添加