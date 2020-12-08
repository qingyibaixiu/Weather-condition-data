function [Psite] = Wind_model(vx)
% 
%��������
a=-3.3872;%
b=87.7673;
c=-505.1193;
d=871.6202;
H1=9;%������߶�,m
H2=70;%�����챸߶�,m
Pr_site=1.004*1e5;%ʵ�ʻ�����ѹ,Pa
T_site=290.15;%ʵ�ʻ����¶�,K
alpha=1/7;%����ָ����;
vci=4;%�������,m/s
vr=12;%�����,m/s
vco=25;%�г�����,m/s
Pr_SL=1.01*1e5;%��׼���������µĴ���ѹ,Pa
T_SL=293;%��׼���������µ��¶�,K
biaozhunmidu=1.225;%��׼���������µĿ����ܶ�,kg/m^3
Pr=1650;%�����,kW
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

