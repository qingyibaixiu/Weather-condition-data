function [Ppv] = PV_model(jingdu,weidu,fangweijiao,qingxiejiao,qinlang,fanzhao,m,Ppv_rate,fpv)
%光伏发电模型计算(不考虑温度对光伏阵列的影响)，本次优化调度的仿真步长为1小时
%jingdu:光伏阵列所在地的经度
%weidu:光伏阵列所在地的纬度
%fangweijiao:水平面的方位角
%qingxiejiao:光伏阵列倾斜角
%qinlang:晴朗指数
%fanzhao:反照率
%m:当日距离1月1日的天数
%Ppv-rate:光伏阵列的额定功率
%fpv:功率降额因数,通常取0.8
%*********初始设定参数**************
aveGtstc=1;%标准测试条件下的太阳辐照强度，一般取1kW/m^2
Zc=8;%光伏电池所在地的时区，是本初子午线以东的时区，中国为8
for i=1:24
    tc(i)=i-1;%光伏电池所在地时间，单位为hr
    tc2(i)=i;
end
Gsc=1.367;%太阳常数，一般取1.367kW/m^2
weidu=weidu/360*2*pi;
qingxiejiao=qingxiejiao/360*2*pi;
fangweijiao=fangweijiao/360*2*pi;
%*********模型计算公式*********
chiwei=23.45/360*2*pi*sin(2*pi*(284+m)/365);%太阳赤纬计算公式
B=2*pi*(m-1)/365;
E=3.82*(7.5e-5+1.868e-3*cos(B)-0.032077*sin(B)-0.014615*cos(2*B)-0.04089*sin(2*B));
for i=1:24
    ts(i)=tc(i)+jingdu/15-Zc+E;%太阳时间计算公式
    ts2(i)=tc2(i)+jingdu/15-Zc+E;
    xiaoshijiao(i)=15*(ts(i)-12);%小时角计算公式
    xiaoshijiao2(i)=15*(ts2(i)-12);
    xiaoshijiao1(i)=xiaoshijiao(i);
    cos_rushejiao(i)=sin(chiwei)*sin(weidu)*cos(qingxiejiao)-sin(chiwei)*cos(weidu)*sin(qingxiejiao)*cos(fangweijiao)...
    +cos(chiwei)*cos(weidu)*cos(qingxiejiao)*cos(xiaoshijiao(i)/360*2*pi)+cos(chiwei)*sin(weidu)*sin(qingxiejiao)*cos(fangweijiao)*cos(xiaoshijiao(i)/360*2*pi)...
    +cos(chiwei)*sin(qingxiejiao)*sin(fangweijiao)*sin(xiaoshijiao(i)/360*2*pi);%光伏阵列入射角计算公式
    cos_rushejiao0(i)=cos(weidu)*cos(chiwei)*cos(xiaoshijiao(i)/360*2*pi)+sin(weidu)*sin(chiwei);
    Gon=Gsc*(1+0.033*cos(2*pi*m/365));
    Go=Gon*cos_rushejiao0(i);
    aveGo(i)=12/pi*Gon*(cos(weidu)*cos(chiwei)*(sin(xiaoshijiao1(i)/360*2*pi)-sin(xiaoshijiao2(i)/360*2*pi))+pi*(xiaoshijiao1(i)-xiaoshijiao2(i))/180*sin(weidu)*sin(chiwei));
    aveG(i)=aveGo(i)*qinlang;
    if qinlang<=0.22
        aveGd(i)=(1-0.09*qinlang)*aveG(i);
    else if qinlang>0.22&qinlang<=0.8
            aveGd(i)=aveG(i)*(0.9511-0.1604*qinlang+4.388*qinlang^2-16.638*qinlang^3+12.336*qinlang^4);
        else
            aveGd(i)=0.165*aveG(i);
        end
    end
    aveGb(i)=aveG(i)-aveGd(i);
    Rb(i)=cos_rushejiao(i)/cos_rushejiao0(i);
    Ai(i)=aveGb(i)/aveGo(i);
    fp(i)=sqrt(aveGb(i)/aveG(i));
    aveGt(i)=Rb(i)*(aveGb(i)+aveGd(i)*Ai(i))+aveG(i)*fanzhao*((1-cos(qingxiejiao))/2)+aveGd(i)*(1-Ai(i))*((1+cos(qingxiejiao))/2)*(1+fp(i)*sin(qingxiejiao/2)*sin(qingxiejiao/2)*sin(qingxiejiao/2));
    Ppv(i)=fpv*Ppv_rate*(aveGt(i)/aveGtstc);
end
end

