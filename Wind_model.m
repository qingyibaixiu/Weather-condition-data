function [Psite] = Wind_model(vx)
% 
%基本参数
a=-3.3872;%
b=87.7673;
c=-505.1193;
d=871.6202;
H1=9;%测量点高度,m
H2=70;%风机轮毂高度,m
Pr_site=1.004*1e5;%实际环境气压,Pa
T_site=290.15;%实际环境温度,K
alpha=1/7;%幂律指数α;
vci=4;%切入风速,m/s
vr=12;%额定风速,m/s
vco=25;%切出风速,m/s
Pr_SL=1.01*1e5;%标准测试条件下的大气压,Pa
T_SL=293;%标准测试条件下的温度,K
biaozhunmidu=1.225;%标准测试条件下的空气密度,kg/m^3
Pr=1650;%额定功率,kW
for i=1:24
v(i)=vx(i)*(H2/H1)^alpha;
end
for i=1:24
if v(i)<vci
    Pwt(i)=0;
else if v(i)<=vr
        Pwt(i)=a*v(i)^3+b*v(i)^2+c*v(i)+d;
    else if v(i)<vco
            Pwt(i)=Pr;
        else
            Pwt(i)=0;
        end
    end
end
midu=biaozhunmidu*Pr_site*T_SL/(Pr_SL*T_site);
Psite=Pwt*midu/biaozhunmidu;           
end

