function [Ppv] = PV_model(jingdu,weidu,fangweijiao,qingxiejiao,qinlang,fanzhao,m,Ppv_rate,fpv)
%�������ģ�ͼ���(�������¶ȶԹ�����е�Ӱ��)�������Ż����ȵķ��沽��Ϊ1Сʱ
%jingdu:����������ڵصľ���
%weidu:����������ڵص�γ��
%fangweijiao:ˮƽ��ķ�λ��
%qingxiejiao:���������б��
%qinlang:����ָ��
%fanzhao:������
%m:���վ���1��1�յ�����
%Ppv-rate:������еĶ����
%fpv:���ʽ�������,ͨ��ȡ0.8
%*********��ʼ�趨����**************
aveGtstc=1;%��׼���������µ�̫������ǿ�ȣ�һ��ȡ1kW/m^2
Zc=8;%���������ڵص�ʱ�����Ǳ����������Զ���ʱ�����й�Ϊ8
for i=1:24
    tc(i)=i-1;%���������ڵ�ʱ�䣬��λΪhr
    tc2(i)=i;
end
Gsc=1.367;%̫��������һ��ȡ1.367kW/m^2
weidu=weidu/360*2*pi;
qingxiejiao=qingxiejiao/360*2*pi;
fangweijiao=fangweijiao/360*2*pi;
%*********ģ�ͼ��㹫ʽ*********
chiwei=23.45/360*2*pi*sin(2*pi*(284+m)/365);%̫����γ���㹫ʽ
B=2*pi*(m-1)/365;
E=3.82*(7.5e-5+1.868e-3*cos(B)-0.032077*sin(B)-0.014615*cos(2*B)-0.04089*sin(2*B));
for i=1:24
    ts(i)=tc(i)+jingdu/15-Zc+E;%̫��ʱ����㹫ʽ
    ts2(i)=tc2(i)+jingdu/15-Zc+E;
    xiaoshijiao(i)=15*(ts(i)-12);%Сʱ�Ǽ��㹫ʽ
    xiaoshijiao2(i)=15*(ts2(i)-12);
    xiaoshijiao1(i)=xiaoshijiao(i);
    cos_rushejiao(i)=sin(chiwei)*sin(weidu)*cos(qingxiejiao)-sin(chiwei)*cos(weidu)*sin(qingxiejiao)*cos(fangweijiao)...
    +cos(chiwei)*cos(weidu)*cos(qingxiejiao)*cos(xiaoshijiao(i)/360*2*pi)+cos(chiwei)*sin(weidu)*sin(qingxiejiao)*cos(fangweijiao)*cos(xiaoshijiao(i)/360*2*pi)...
    +cos(chiwei)*sin(qingxiejiao)*sin(fangweijiao)*sin(xiaoshijiao(i)/360*2*pi);%�����������Ǽ��㹫ʽ
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

