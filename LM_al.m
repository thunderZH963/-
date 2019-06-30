%待训练数据集  可修改 
data_1=[8.90 8.93 8.97 9.00 9.03 9.07 9.10 9.13 9.17 9.20 9.23 9.27 9.30 9.33 9.37 9.40 9.43 9.47 9.50 9.53];      % 北京时间
obs_1=[2.37 2.34 2.32 2.29 2.25 2.22 2.19 2.17 2.14 2.11 2.08 2.05 2.03 2.00 1.97 1.94 1.92 1.89 1.86 1.84];    %影子长度

%作图
figure(1)
plot(data_1,obs_1);    %data_1作为横坐标  obs_1作为纵坐标
xlabel('北京时间/h');
ylabel('影子长度/m');

%n取114，代表距离春分114天
n = 114;

%直射点纬度
%利用黄经推导 参考文献  蒋洪力.太阳直射点纬度的数学推导和分析.学通报. 2007 年
latiDir = asin(0.39775*sin(180/186*n*pi/180));

syms stanTime longiNow
%当地时间
%stanTime代表北京时间
localTime = stanTime - (120 - longiNow / pi * 180) * 4 / 60;

%拍摄地与太阳直射点经度差
delta = (12 - localTime)*15*pi/180;

%杆高2m 
%latiNow 代表当地纬度
syms h latiNow

%得到影子坐标
%关于坐标系的说明参考.py文件
%x
x =  -1*cos(latiDir)*sin(delta)*h/(cos(latiNow)*cos(latiDir)*cos(delta) + sin(latiNow)*sin(latiDir));
%y
y =  (sin(latiNow)*cos(latiDir)*cos(delta) - cos(latiNow)*sin(latiDir))*h/(cos(latiNow)*cos(latiDir)*cos(delta) + sin(latiNow)*sin(latiDir));

%h==2
x = x/h*2;
y = y/h*2;

%坐标变换   此处不需要
%turnAngle
%syms turnAngle
%x_new = x*cos(turnAngle) + y * sin(turnAngle)
%y_new = y*cos(turnAngle) - x * sin(turnAngle)

%影子长度
length = sqrt(x^2+y^2);

%加入拍摄角度干扰因子
%turnAngle
syms turnAngle
length_new = length*cos(turnAngle)

%存储逐次的参数
a1 = [];
a2 = [];
a3 = [];
%初始猜想  
turnAngle_0 = 0.2*pi/180;
longiNow_0 = 100*pi/180;
latiNow_0 = 30*pi/180;
a = inline(length_new);
y_init = a(latiNow_0,longiNow_0,data_1,turnAngle_0); 
Ndata = 20;     %len(data_1)
Nparams = 3;    %3个参数
a1 = [a1,turnAngle_0/pi*180];
a2 = [a2,longiNow_0/pi*180];
a3 = [a3,latiNow_0/pi*180];

%雅可比矩阵
Jsym = jacobian(length_new,[turnAngle,longiNow,latiNow]);

%迭代次数
n_iters = 500;

%阻尼系数
lamda = 0.01;

updateJ=1;

%参与迭代的参数设置
turnAngle_est = turnAngle_0;
longiNow_est = longiNow_0;
latiNow_est = latiNow_0;

% 迭代开始
for it = 1:n_iters
    if updateJ == 1
        % 根据当前估计值，计算雅克比矩阵
        J=zeros(Ndata,Nparams);
        for i=1:20      %根据数据集修改
            b=inline(Jsym);
            J(i,:)= b(latiNow_est,longiNow_est,data_1(i),turnAngle_est);
        end
        % 根据当前参数，得到函数值
        y_est = a(latiNow_est,longiNow_est,data_1,turnAngle_est);
        % 计算误差
        d=obs_1-y_est;
        % 计算海塞矩阵
        H=J'*J;
        % 若是第一次迭代，计算误差
        if it==1
            e=dot(d,d);
        end
    end
    % 根据阻尼系数lamda混合得到H矩阵
    H_lm=H+(lamda*eye(Nparams,Nparams));
    
    % 计算步长dp，并根据步长计算新的可能的\参数估计值
    
    g = J'*d(:);
    dp=H_lm\g;
    turnAngle_lm = turnAngle_est+dp(1);
    longiNow_lm =  longiNow_est +dp(2);
    latiNow_lm =  latiNow_est + dp(3);
    %turnAngle_lm= turnAngle_lm + 180* fix(turnAngle_lm/180);
    %longiNow_lm = longiNow_lm - 180 * fix(longiNow_lm/180);
    %latiNow_lm = latiNow_lm - 90 * fix(latiNow_lm/90);
    % 计算新的可能估计值对应的y和计算残差e
    y_est_lm = a(latiNow_lm,longiNow_lm,data_1,turnAngle_lm);
    d_lm=obs_1-y_est_lm;
    e_lm=dot(d_lm,d_lm);
    
    % 根据误差，决定如何更新参数和阻尼系数
    if e_lm  < e
        lamda=lamda/10;
        turnAngle_est= turnAngle_lm;
        longiNow_est=  longiNow_lm;
        latiNow_est =  latiNow_lm;
        a1 = [a1,turnAngle_est/pi*180];
        a2 = [a2,-1*longiNow_est/pi*180];
        a3 = [a3,latiNow_est/pi*180];
        e=e_lm;
        disp(e);
        updateJ=1;
    else
        updateJ=0;
        lamda=lamda*10;
    end
end

%作图 迭代次数
[temp1,temp2] = size(a1);
k = [1:temp2];
figure(2)
plot(k,a1);
xlabel('迭代次数');
ylabel('角度');
figure(3)
plot(k,a2);
xlabel('迭代次数');
ylabel('经度');
figure(4)
plot(k,a3);
xlabel('迭代次数');
ylabel('纬度');

%显示优化的结果
fprintf("预测纬度为%f\n",a3(temp2));
fprintf("预测经度为%f\n",a2(temp2));

