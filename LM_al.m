%��ѵ�����ݼ�  ���޸� 
data_1=[8.90 8.93 8.97 9.00 9.03 9.07 9.10 9.13 9.17 9.20 9.23 9.27 9.30 9.33 9.37 9.40 9.43 9.47 9.50 9.53];      % ����ʱ��
obs_1=[2.37 2.34 2.32 2.29 2.25 2.22 2.19 2.17 2.14 2.11 2.08 2.05 2.03 2.00 1.97 1.94 1.92 1.89 1.86 1.84];    %Ӱ�ӳ���

%��ͼ
figure(1)
plot(data_1,obs_1);    %data_1��Ϊ������  obs_1��Ϊ������
xlabel('����ʱ��/h');
ylabel('Ӱ�ӳ���/m');

%nȡ114��������봺��114��
n = 114;

%ֱ���γ��
%���ûƾ��Ƶ� �ο�����  ������.̫��ֱ���γ�ȵ���ѧ�Ƶ��ͷ���.ѧͨ��. 2007 ��
latiDir = asin(0.39775*sin(180/186*n*pi/180));

syms stanTime longiNow
%����ʱ��
%stanTime������ʱ��
localTime = stanTime - (120 - longiNow / pi * 180) * 4 / 60;

%�������̫��ֱ��㾭�Ȳ�
delta = (12 - localTime)*15*pi/180;

%�˸�2m 
%latiNow ������γ��
syms h latiNow

%�õ�Ӱ������
%��������ϵ��˵���ο�.py�ļ�
%x
x =  -1*cos(latiDir)*sin(delta)*h/(cos(latiNow)*cos(latiDir)*cos(delta) + sin(latiNow)*sin(latiDir));
%y
y =  (sin(latiNow)*cos(latiDir)*cos(delta) - cos(latiNow)*sin(latiDir))*h/(cos(latiNow)*cos(latiDir)*cos(delta) + sin(latiNow)*sin(latiDir));

%h==2
x = x/h*2;
y = y/h*2;

%����任   �˴�����Ҫ
%turnAngle
%syms turnAngle
%x_new = x*cos(turnAngle) + y * sin(turnAngle)
%y_new = y*cos(turnAngle) - x * sin(turnAngle)

%Ӱ�ӳ���
length = sqrt(x^2+y^2);

%��������Ƕȸ�������
%turnAngle
syms turnAngle
length_new = length*cos(turnAngle)

%�洢��εĲ���
a1 = [];
a2 = [];
a3 = [];
%��ʼ����  
turnAngle_0 = 0.2*pi/180;
longiNow_0 = 100*pi/180;
latiNow_0 = 30*pi/180;
a = inline(length_new);
y_init = a(latiNow_0,longiNow_0,data_1,turnAngle_0); 
Ndata = 20;     %len(data_1)
Nparams = 3;    %3������
a1 = [a1,turnAngle_0/pi*180];
a2 = [a2,longiNow_0/pi*180];
a3 = [a3,latiNow_0/pi*180];

%�ſɱȾ���
Jsym = jacobian(length_new,[turnAngle,longiNow,latiNow]);

%��������
n_iters = 500;

%����ϵ��
lamda = 0.01;

updateJ=1;

%��������Ĳ�������
turnAngle_est = turnAngle_0;
longiNow_est = longiNow_0;
latiNow_est = latiNow_0;

% ������ʼ
for it = 1:n_iters
    if updateJ == 1
        % ���ݵ�ǰ����ֵ�������ſ˱Ⱦ���
        J=zeros(Ndata,Nparams);
        for i=1:20      %�������ݼ��޸�
            b=inline(Jsym);
            J(i,:)= b(latiNow_est,longiNow_est,data_1(i),turnAngle_est);
        end
        % ���ݵ�ǰ�������õ�����ֵ
        y_est = a(latiNow_est,longiNow_est,data_1,turnAngle_est);
        % �������
        d=obs_1-y_est;
        % ���㺣������
        H=J'*J;
        % ���ǵ�һ�ε������������
        if it==1
            e=dot(d,d);
        end
    end
    % ��������ϵ��lamda��ϵõ�H����
    H_lm=H+(lamda*eye(Nparams,Nparams));
    
    % ���㲽��dp�������ݲ��������µĿ��ܵ�\��������ֵ
    
    g = J'*d(:);
    dp=H_lm\g;
    turnAngle_lm = turnAngle_est+dp(1);
    longiNow_lm =  longiNow_est +dp(2);
    latiNow_lm =  latiNow_est + dp(3);
    %turnAngle_lm= turnAngle_lm + 180* fix(turnAngle_lm/180);
    %longiNow_lm = longiNow_lm - 180 * fix(longiNow_lm/180);
    %latiNow_lm = latiNow_lm - 90 * fix(latiNow_lm/90);
    % �����µĿ��ܹ���ֵ��Ӧ��y�ͼ���в�e
    y_est_lm = a(latiNow_lm,longiNow_lm,data_1,turnAngle_lm);
    d_lm=obs_1-y_est_lm;
    e_lm=dot(d_lm,d_lm);
    
    % ������������θ��²���������ϵ��
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

%��ͼ ��������
[temp1,temp2] = size(a1);
k = [1:temp2];
figure(2)
plot(k,a1);
xlabel('��������');
ylabel('�Ƕ�');
figure(3)
plot(k,a2);
xlabel('��������');
ylabel('����');
figure(4)
plot(k,a3);
xlabel('��������');
ylabel('γ��');

%��ʾ�Ż��Ľ��
fprintf("Ԥ��γ��Ϊ%f\n",a3(temp2));
fprintf("Ԥ�⾭��Ϊ%f\n",a2(temp2));

