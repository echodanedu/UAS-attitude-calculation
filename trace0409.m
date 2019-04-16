%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                    ����������
%
%  �������:
%  t          ����ʱ��
%  T          ���沽��
%  atti       ���������������ǣ���λ���ȣ�
%  atti_rate  ������ʡ��������ʡ���������ʣ���λ����/�룩
%  veloB      �ɻ��˶��ٶȡ���X����Y��ͷ��Z���򣨵�λ����/�룩
%  acceB      �ɻ��˶����ٶȡ���X����Y��ͷ��Z���򣨵�λ����/��/�룩
%  posi       ������������ʼλ�þ��ȡ�γ�ȡ��߶ȣ���λ���ȡ��ȡ��ף�
%  veloN      �ɻ��˶��ٶȡ���X����Y����Z���򣨵�λ����/�룩
%  posiN      ���ȡ�γ�ȡ��߶ȣ���λ���ȡ��ȡ��ף�        
%                           ������ƣ�����  ���ڣ�2007/12/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc
%%%%%%%%����ʱ������%%%%%%
t=0;
T=1;
TraceData=[];
t_stop=50.0;

%%%%%%%%%%%%%%%%����������%%%%%%%%%%%%%%%%%
atti=zeros(3,1);     %���������������ǣ���λ���ȣ�
atti_rate=zeros(3,1);%��������ʡ��������ʡ��������ʣ���λ����/�룩
atti_rate1=zeros(3,1);
delta_theta1=zeros(3,1);
delta_theta2=zeros(3,1);
veloB=zeros(3,1);    %�ɻ��˶��ٶȡ���X����Y��ͷ��Z���򣨵�λ����/�룩
acceB=zeros(3,1);    %�ɻ��˶����ٶȡ���X����Y��ͷ��Z���򣨵�λ����/��/�룩
posi=zeros(3,1);     %������������ʼλ�þ��ȡ�γ�ȡ��߶ȣ���λ���ȡ��ȡ��ף�
posi=[118.78;32.05;300.0];  

phi = zeros(1,1);                          %������
psi = zeros(1,1);                          %ƫ����
gamma = zeros(1,1);                        %�����

phi1 = zeros(1,1);                          %������
psi1 = zeros(1,1);                          %ƫ����
gamma1 = zeros(1,1);                        %�����

phi3 = zeros(1,t_stop);                           %����ĸ�����
psi3 = zeros(1,t_stop);                           %����ĺ����
gamma3 = zeros(1,t_stop);                         %����ĺ����
psi3(1) = pi/2;

Re=6378137.0;  %����뾶���ף� 
f=1/298.257;   %�������Բ��
posiN=posi;    %��ʼλ���뺽��λ��һ��
long=posiN(1,1)*pi/180.0;lati=posiN(2,1)*pi/180.0;heig=posiN(3,1);
    %������λ��

Rm=Re*(1-2*f+3*f*sin(lati)*sin(lati));
Rn=Re*(1+f*sin(lati)*sin(lati));
    %�������ʰ뾶���
atti(1,1)=0.0;
atti(2,1)=0.0;
atti(3,1)=90.0; 

Q3 =  ...                                     %��Ԫ����ֵ
   [cos(0.5*gamma3(1))*cos(0.5*psi3(1))*cos(0.5*phi3(1))+sin(0.5*gamma3(1))*sin(0.5*psi3(1))*sin(0.5*phi3(1));
    sin(0.5*gamma3(1))*cos(0.5*psi3(1))*cos(0.5*phi3(1))-cos(0.5*gamma3(1))*sin(0.5*psi3(1))*sin(0.5*phi3(1));
    cos(0.5*gamma3(1))*sin(0.5*psi3(1))*cos(0.5*phi3(1))+sin(0.5*gamma3(1))*cos(0.5*psi3(1))*sin(0.5*phi3(1));
    cos(0.5*gamma3(1))*cos(0.5*psi3(1))*sin(0.5*phi3(1))-sin(0.5*gamma3(1))*sin(0.5*psi3(1))*cos(0.5*phi3(1))];

veloB(2,1)=0.0;  %�ɻ���ʼ�˶��ٶȡ�����ͷ����λ����/�룩
while t<=t_stop
%      if(t==0 )
%       acceB(2,1)=0.0;               %��ʼ��׼ʱ��
%    elseif(t<=10)
%       acceB(2,1)=2.0;
%    elseif(t<=50)
%       acceB(2,1)=0.0;
%    elseif(t<=55)
%       acceB(2,1)=-2.0;
%    elseif(t<=60)
%       acceB(2,1)=0.0;
%       atti_rate(3,1)=-18.0;
%    elseif(t<=65)
%       acceB(2,1)=2.0;
%       atti_rate(3,1)=0.0;
%    elseif(t<=105)
%       acceB(2,1)=0.0;
%   elseif(t<=110)
%       acceB(2,1)=-2.0;
%   elseif(t<=115)
%       acceB(2,1)=0.0;
%       atti_rate(3,1)=-18.0;
%   elseif(t<=120)
%       acceB(2,1)=2.0;
%       atti_rate(3,1)=0.0; 
%   elseif (t<=125)                  
%       atti_rate(1,1)=9.0;
%    elseif (t<=240)                  
%       atti_rate(1,1)=0.0;
%    elseif (t<=245)                  
%       atti_rate(1,1)=-9.0;
%    elseif (t<=360)                  
%       atti_rate(1,1)=0.0;
%    elseif(t<=365)
%       atti_rate(1,1)=-9.0;
%    elseif(t<=480)
%       atti_rate(1,1)=0.0;
%    elseif(t<=485)
%       atti_rate(1,1)=9.0;
%    else
%       atti_rate(1,1)=0.0;
%    end

  if(t==0 )
      acceB(2,1)=0.0;               %��ʼ��׼ʱ��
  elseif (t<=5)                  
      atti_rate(2,1)=9.0;
   elseif (t<=10)                  
      atti_rate(2,1)=0.0;
   elseif (t<=15)                  
      atti_rate(2,1)=-9.0;
   elseif (t<=20)                  
      atti_rate(2,1)=0.0;
   elseif(t<=25)
      atti_rate(2,1)=-9.0;
   elseif(t<=30)
      atti_rate(2,1)=0.0;
   elseif(t<=35)
      atti_rate(2,1)=9.0;
   else
      atti_rate(2,1)=0.0;
   end
   
   t=t+T;
   t
   
   veloB(2,1)=veloB(2,1)+acceB(2,1)*T;    

   atti(1,1)=atti(1,1)+atti_rate(1,1)*T;
   atti(2,1)=atti(2,1)+atti_rate(2,1)*T;
   atti(3,1)=atti(3,1)+atti_rate(3,1)*T;
  
   roll=atti(1,1)*pi/180.0;pitch=atti(2,1)*pi/180.0;head=atti(3,1)*pi/180.0;

%   ����ϵN(����ϵ)-->B(����ϵ)
%   Cnb=[cos(roll)*cos(head)+sin(roll)*sin(pitch)*sin(head), -cos(roll)*sin(head)+sin(roll)*sin(pitch)*cos(head), -sin(roll)*cos(pitch);
%        cos(pitch)*sin(head),                               cos(pitch)*cos(head),                                sin(pitch);
%        sin(roll)*cos(head)-cos(roll)*sin(pitch)*sin(head), -sin(roll)*sin(head)-cos(roll)*sin(pitch)*cos(head), cos(roll)*cos(pitch)];
% 
%    Cbn = Cnb';
   
%    Cbn=[cos(roll)*cos(head)+sin(roll)*sin(pitch)*sin(head),  cos(pitch)*sin(head),    sin(roll)*cos(head)-cos(roll)*sin(pitch)*sin(head);
%        -cos(roll)*sin(head)+sin(roll)*sin(pitch)*cos(head),  cos(pitch)*cos(head),    -sin(roll)*sin(head)-cos(roll)*sin(pitch)*cos(head);
%        -sin(roll)*cos(pitch),                                sin(pitch),              cos(roll)*cos(pitch)];            %1,3�н����ߵ�����һ��

   
  Cnb=[cos(roll)*cos(head)-sin(roll)*sin(pitch)*sin(head), cos(roll)*sin(head)+sin(roll)*sin(pitch)*cos(head), -sin(roll)*cos(pitch);
       -cos(pitch)*sin(head),                               cos(pitch)*cos(head),                                sin(pitch);
       sin(roll)*cos(head)+cos(roll)*sin(pitch)*sin(head), sin(roll)*sin(head)-cos(roll)*sin(pitch)*cos(head), cos(roll)*cos(pitch)];
  Cbn = Cnb';                 %0412 �Ƶ�

  %������������������ٶȣ�%
  phi = asin(Cnb(2,3))*180/pi;                      %������
  psi = atan2(Cnb(2,1),Cnb(2,2))*180/pi;            %ƫ����
  gamma = atan2(-Cnb(1,3),Cbn(3,3))*180/pi;         %�����
% 
%   C =  [-sin(psi) 0 1;
%          sin(gamma)*cos(psi) cos(gamma) 0;
%          cos(gamma)*cos(psi) -sin(gamma) 0];
% %   B = [(phi-phi1);(psi-psi1);(gamma-gamma1)]; 
  B = [(gamma-gamma1);(phi-phi1);(psi-psi1)];
  
  C1 =  [ 0  cos(gamma)  -sin(gamma)*cos(phi);
          1      0          sin(phi);
          0  sin(gamma) cos(gamma)*cos(phi)];    
     
  wibb = C1*atti_rate*pi/180;                        %����������Ľ��ٶ�  
%   wibb1 = C1*atti_rate1 *pi/180;  
%   atti_rate1 = atti_rate;
  
  phi1 = phi;             %������
  psi1 = psi;             %ƫ����
  gamma1 = gamma;         %�����
  
  wnbb3 = wibb;
  
  delta_big_theta = [     0        -wnbb3(1)   -wnbb3(2) -wnbb3(3);
                     wnbb3(1)        0          wnbb3(3) -wnbb3(2);
                     wnbb3(2)  -wnbb3(3)         0        wnbb3(1);
                     wnbb3(3)   wnbb3(2)   -wnbb3(1)       0      ];
                 
  delta_theta = sqrt(wnbb3(1)^2+wnbb3(2)^2+wnbb3(3)^2);
  Q3 = (eye(4)*(1-delta_theta^2/8+delta_theta^4/384)+(1/2-delta_theta^2/48)*delta_big_theta)*Q3;
  
%     delta_theta1 = wibb';
%     delta_theta2 = wibb1';
%   
%     big_phi = [delta_theta1(1)+delta_theta2(1)+2/3*(delta_theta1(2)*delta_theta2(3)-delta_theta1(3)*delta_theta2(2));
%                 delta_theta1(2)+delta_theta2(2)+2/3*(delta_theta1(3)*delta_theta2(1)-delta_theta1(1)*delta_theta2(3));
%                 delta_theta1(3)+delta_theta2(3)+2/3*(delta_theta1(2)*delta_theta2(1)+delta_theta1(1)*delta_theta2(2))];
%            
%     big_phi_val = sqrt(big_phi(1)^2 + big_phi(2)^2 + big_phi(3)^2);
%     
%     if(big_phi_val > 0.00001)
%         Q3 =[cos(big_phi_val/2)                           -big_phi(1)/big_phi_val*sin(big_phi_val/2)     -big_phi(2)/big_phi_val*sin(big_phi_val/2)   -big_phi(3)/big_phi_val*sin(big_phi_val/2);
%             big_phi(1)/big_phi_val*sin(big_phi_val/2)    cos(big_phi_val/2)                             big_phi(3)/big_phi_val*sin(big_phi_val/2)    -big_phi(2)/big_phi_val*sin(big_phi_val/2);
%             big_phi(2)/big_phi_val*sin(big_phi_val/2)    big_phi(3)/big_phi_val*sin(big_phi_val/2)      cos(big_phi_val/2)                           -big_phi(2)/big_phi_val*sin(big_phi_val/2);
%             big_phi(3)/big_phi_val*sin(big_phi_val/2)    big_phi(2)/big_phi_val*sin(big_phi_val/2)      -big_phi(1)/big_phi_val*sin(big_phi_val/2)   cos(big_phi_val/2)]*Q3;
%     end

  Q3_val = sqrt(Q3(1)^2+Q3(2)^2+Q3(3)^2+Q3(4)^2);
  
  if Q3_val < 0.000001
      Q3 = [1;0;0;0];
  else
      Q3 = Q3/Q3_val;
  end
  
  %������̬����%
%   Cnb3 = ...                                
%          [Q3(1)^2+Q3(2)^2-Q3(3)^2-Q3(4)^2     2*(Q3(2)*Q3(3)+Q3(1)*Q3(4))       2*(Q3(2)*Q3(4)-Q3(1)*Q3(3));
%             2*(Q3(2)*Q3(3)-Q3(1)*Q3(4))     Q3(1)^2-Q3(2)^2+Q3(3)^2-Q3(4)^2     2*(Q3(3)*Q3(4)+Q3(1)*Q3(2));
%             2*(Q3(2)*Q3(4)+Q3(1)*Q3(3))       2*(Q3(3)*Q3(4)-Q3(1)*Q3(2))     Q3(1)^2-Q3(2)^2-Q3(3)^2+Q3(4)^2];
%   Cbn3 = Cnb3';                              
  
  Cbn3 = ...                                 
       [Q3(1)^2+Q3(2)^2-Q3(3)^2-Q3(4)^2     2*(Q3(2)*Q3(3)-Q3(1)*Q3(4))       2*(Q3(2)*Q3(4)+Q3(1)*Q3(3));
        2*(Q3(2)*Q3(3)+Q3(1)*Q3(4))     Q3(1)^2-Q3(2)^2+Q3(3)^2-Q3(4)^2       2*(Q3(3)*Q3(4)-Q3(1)*Q3(2));
        2*(Q3(2)*Q3(4)-Q3(1)*Q3(3))          2*(Q3(3)*Q3(4)+Q3(1)*Q3(2))     Q3(1)^2-Q3(2)^2-Q3(3)^2+Q3(4)^2];
    
  A = [Cbn;Cbn3]
  %������̬��%

% phi3(i) = asin(2*(Q3(3)*Q3(4)-Q3(1)*Q3(2)));
% gamma3(i) = atan2(2*(Q3(2)*Q3(3)+Q3(1)*Q3(4)),Q3(1)^2-Q3(2)^2+Q3(3)^2-Q3(4)^2);
% psi3(i) = atan2(2*(Q3(2)*Q3(4)+Q3(1)*Q3(3)),Q3(1)^2-Q3(2)^2-Q3(3)^2+Q3(4)^2);
 
phi3 = asin(2*(Q3(3)*Q3(4)+Q3(1)*Q3(2)));
gamma3 = atan2(-2*(Q3(2)*Q3(4)-Q3(1)*Q3(3)), Q3(1)^2-Q3(2)^2-Q3(3)^2+Q3(4)^2);
psi3 = atan2(-2*(Q3(2)*Q3(3)-Q3(1)*Q3(4)),Q3(1)^2-Q3(2)^2+Q3(3)^2-Q3(4)^2);
  
  veloN=Cbn*veloB;
  Ve=veloN(1,1);Vn=veloN(2,1);Vu=veloN(3,1);
  V=sqrt(veloN(1,1)^2+veloN(2,1)^2+veloN(3,1)^2);

  heig=heig+T*Vu;
  lati=lati+T*(Vn/(Rm+heig));
  long=long-T*(Ve/((Rn+heig)*cos(lati)));

  posiN(1,1)=long*180.0/pi;       %��λ����
  posiN(2,1)=lati*180.0/pi;       %��λ����
  posiN(3,1)=heig;

  TraceData=[TraceData;t,posiN',V,veloN',atti',wibb',acceB',gamma3,phi3,psi3];
                      %   2,3,4  6,7,8  9,10,11 12,13,14     18,   19,  20
end


% fg_1=figure('Name','Trace','NumberTitle','off');
% plot3(TraceData(:,2),TraceData(:,3),TraceData(:,4),'m',TraceData(1,2),TraceData(1,3),TraceData(1,4),'o',TraceData(490,2),TraceData(490,3),TraceData(490,4),'o');
% title('���к�������');
% xlabel('���ȣ�^{\circ}��');ylabel('γ�ȣ�^{\circ}��');zlabel('�߶ȣ�m��');
% grid;
% fg_2=figure('Name','trace','NumberTitle','off');
% subplot(3,1,1);plot(TraceData(:,1),TraceData(:,2));title('����������γ���߶ȷ���');ylabel('���ȣ�^{\circ}��');grid;
% subplot(3,1,2);plot(TraceData(:,1),TraceData(:,3));ylabel('γ�ȣ�^{\circ}��');grid;
% subplot(3,1,3);plot(TraceData(:,1),TraceData(:,4));xlabel('t(sec)');ylabel('�߶ȣ�m��');grid;
% fg_3=figure('Name','Velo','NumberTitle','off');
% plot(TraceData(:,1),TraceData(:,5));
% title('�����ٶȷ���');
% xlabel('t(sec)');ylabel('�ٶȣ�m/s��');grid;
% fg_4=figure('Name','velo','NumberTitle','off');
% subplot(3,1,1);plot(TraceData(:,1),TraceData(:,6));title('�����������������ٶȷ���');ylabel('�����ٶ�(m/s)');grid;
% subplot(3,1,2);plot(TraceData(:,1),TraceData(:,7));ylabel('�����ٶ�(m/s)');grid;
% subplot(3,1,3);plot(TraceData(:,1),TraceData(:,8));xlabel('t(sec)');ylabel('�����ٶ�(m/s)');grid;

%%%%%%%%% ��̬����  %%%%%%%%%%%
fg_5=figure('Name','attitude','NumberTitle','off');
subplot(3,1,1);plot(TraceData(:,1),TraceData(:,9),'b');title('������������������������');ylabel('���(��)');grid;
hold on  
% plot([0:1:t_stop],gamma3*180/pi,'r');
plot(TraceData(:,1),TraceData(:,18),'r');

subplot(3,1,2);plot(TraceData(:,1),TraceData(:,10));ylabel('����(��)');grid;
hold on  
% plot([0:1:t_stop],phi3*180/pi,'r');
plot(TraceData(:,1),TraceData(:,19),'r');

subplot(3,1,3);plot(TraceData(:,1),TraceData(:,11));xlabel('t(sec)');ylabel('����(��)');grid;
hold on  
% plot([0:1:t_stop],psi3*180/pi,'r');
plot(TraceData(:,1),TraceData(:,20),'r');
%%%%%%%%%%%%%�洢��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save trace0.dat TraceData -ASCII;    %�洢��������

