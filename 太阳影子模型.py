#模型建立

import math

#自建坐标系下，以拍摄点与经线相切且向东方向为x轴
#以拍摄点与纬线相切且向北方向为y轴

def getLength(n,longiNow,latiNow,stanTime): 
    
    #直射点纬度
    longDir = math.asin(0.39775*math.sin(180/186*n*math.pi/180))

    #当地时间
    localTime = stanTime - (120 - longiNow) * 4 / 60

    #拍摄地与太阳直射点经度差
    delta = (12 - localTime)*15*math.pi/180

    x =  -1*math.cos(longDir)*math.sin(delta)*h/(math.cos(latiNow)*math.cos(longDir)*math.cos(delta) + math.sin(latiNow)*math.sin(longDir))

    y =  (math.sin(latiNow)*math.cos(longDir)*math.cos(delta) - math.cos(latiNow)*math.sin(longDir))*h/(math.cos(latiNow)*math.cos(longDir)*math.cos(delta) + math.sin(latiNow)*math.sin(longDir))

    #影子长度
    length = math.sqrt(x**2+y**2)
    
    print(length)



#可以利用getLength函数 做出各种变化曲线
#一个例子
#影子长度随北京时间变化
i = 0
while i < 24:
    #天数
    n = 10
    #杆的高度
    h = 1
    #拍摄地经度
    longiNow = 70 * math.pi/180  
    #拍摄地纬度
    latiNow = 27 * math.pi/180
    #北京时间
    stanTime = i
    getLength(n,longiNow,latiNow,stanTime)
    i+=0.5

#优化算法
#采用L-M最优化算法
#首先梯度下降法和高斯牛顿法都是最优化方法。其区别之处在于，
#梯度下降法在寻找目标函数极小值时，是沿着反梯度方向进行寻找的。梯度的定义就是指向标量场增长最快的方向，在寻找极小值时，先随便定初始点（x0，y0）然后进行迭代不断寻找直到梯度的模达到预设的要求。但是梯度下降法的缺点之处在于：在远离极小值的地方下降很快，而在靠近极小值的地方下降很慢，靠近的时候可能成zig-zag下降。
#而高斯牛顿法是一种非线性最小二乘最优化方法。其利用了目标函数的泰勒展开式把非线性函数的最小二乘化问题化为每次迭代的线性函数的最小二乘化问题。高斯牛顿法的缺点在于：若初始点距离极小值点过远，迭代步长过大会导致迭代下一代的函数值不一定小于上一代的函数值。
#LM算法在高斯牛顿法中加入了因子μ，当μ大时相当于梯度下降法，μ小时相当于高斯牛顿法。在使用Levenberg-Marquart时，先设置一个比较小的μ值，当发现目标函数反而增大时，将μ增大使用梯度下降法快速寻找，然后再将μ减小使用牛顿法进行寻找。

##算法实现参考matlab文档
##LM_al.m
    
