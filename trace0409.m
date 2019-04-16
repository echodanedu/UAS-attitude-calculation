%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                    航迹发生器
%
%  输入参数:
%  t          仿真时间
%  T          仿真步长
%  atti       横滚、俯仰、航向角（单位：度）
%  atti_rate  横滚速率、俯仰速率、航向角速率（单位：度/秒）
%  veloB      飞机运动速度——X右翼、Y机头、Z天向（单位：米/秒）
%  acceB      飞机运动加速度——X右翼、Y机头、Z天向（单位：米/秒/秒）
%  posi       航迹发生器初始位置经度、纬度、高度（单位：度、度、米）
%  veloN      飞机运动速度——X东向、Y北向、Z天向（单位：米/秒）
%  posiN      经度、纬度、高度（单位：度、度、米）        
%                           程序设计：吴玲  日期：2007/12/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc
%%%%%%%%仿真时间设置%%%%%%
t=0;
T=1;
TraceData=[];
t_stop=50.0;

%%%%%%%%%%%%%%%%航迹发生器%%%%%%%%%%%%%%%%%
atti=zeros(3,1);     %横滚、俯仰、航向角（单位：度）
atti_rate=zeros(3,1);%横滚角速率、俯仰速率、航向速率（单位：度/秒）
atti_rate1=zeros(3,1);
delta_theta1=zeros(3,1);
delta_theta2=zeros(3,1);
veloB=zeros(3,1);    %飞机运动速度——X右翼、Y机头、Z天向（单位：米/秒）
acceB=zeros(3,1);    %飞机运动加速度——X右翼、Y机头、Z天向（单位：米/秒/秒）
posi=zeros(3,1);     %航迹发生器初始位置经度、纬度、高度（单位：度、度、米）
posi=[118.78;32.05;300.0];  

phi = zeros(1,1);                          %俯仰角
psi = zeros(1,1);                          %偏航角
gamma = zeros(1,1);                        %横滚角

phi1 = zeros(1,1);                          %俯仰角
psi1 = zeros(1,1);                          %偏航角
gamma1 = zeros(1,1);                        %横滚角

phi3 = zeros(1,t_stop);                           %解算的俯仰角
psi3 = zeros(1,t_stop);                           %解算的航向角
gamma3 = zeros(1,t_stop);                         %解算的横滚角
psi3(1) = pi/2;

Re=6378137.0;  %地球半径（米） 
f=1/298.257;   %地球的椭圆率
posiN=posi;    %初始位置与航迹位置一致
long=posiN(1,1)*pi/180.0;lati=posiN(2,1)*pi/180.0;heig=posiN(3,1);
    %飞行器位置

Rm=Re*(1-2*f+3*f*sin(lati)*sin(lati));
Rn=Re*(1+f*sin(lati)*sin(lati));
    %地球曲率半径求解
atti(1,1)=0.0;
atti(2,1)=0.0;
atti(3,1)=90.0; 

Q3 =  ...                                     %四元数初值
   [cos(0.5*gamma3(1))*cos(0.5*psi3(1))*cos(0.5*phi3(1))+sin(0.5*gamma3(1))*sin(0.5*psi3(1))*sin(0.5*phi3(1));
    sin(0.5*gamma3(1))*cos(0.5*psi3(1))*cos(0.5*phi3(1))-cos(0.5*gamma3(1))*sin(0.5*psi3(1))*sin(0.5*phi3(1));
    cos(0.5*gamma3(1))*sin(0.5*psi3(1))*cos(0.5*phi3(1))+sin(0.5*gamma3(1))*cos(0.5*psi3(1))*sin(0.5*phi3(1));
    cos(0.5*gamma3(1))*cos(0.5*psi3(1))*sin(0.5*phi3(1))-sin(0.5*gamma3(1))*sin(0.5*psi3(1))*cos(0.5*phi3(1))];

veloB(2,1)=0.0;  %飞机初始运动速度——机头（单位：米/秒）
while t<=t_stop
%      if(t==0 )
%       acceB(2,1)=0.0;               %初始对准时间
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
      acceB(2,1)=0.0;               %初始对准时间
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

%   坐标系N(地理系)-->B(机体系)
%   Cnb=[cos(roll)*cos(head)+sin(roll)*sin(pitch)*sin(head), -cos(roll)*sin(head)+sin(roll)*sin(pitch)*cos(head), -sin(roll)*cos(pitch);
%        cos(pitch)*sin(head),                               cos(pitch)*cos(head),                                sin(pitch);
%        sin(roll)*cos(head)-cos(roll)*sin(pitch)*sin(head), -sin(roll)*sin(head)-cos(roll)*sin(pitch)*cos(head), cos(roll)*cos(pitch)];
% 
%    Cbn = Cnb';
   
%    Cbn=[cos(roll)*cos(head)+sin(roll)*sin(pitch)*sin(head),  cos(pitch)*sin(head),    sin(roll)*cos(head)-cos(roll)*sin(pitch)*sin(head);
%        -cos(roll)*sin(head)+sin(roll)*sin(pitch)*cos(head),  cos(pitch)*cos(head),    -sin(roll)*sin(head)-cos(roll)*sin(pitch)*cos(head);
%        -sin(roll)*cos(pitch),                                sin(pitch),              cos(roll)*cos(pitch)];            %1,3行交换惯导书上一致

   
  Cnb=[cos(roll)*cos(head)-sin(roll)*sin(pitch)*sin(head), cos(roll)*sin(head)+sin(roll)*sin(pitch)*cos(head), -sin(roll)*cos(pitch);
       -cos(pitch)*sin(head),                               cos(pitch)*cos(head),                                sin(pitch);
       sin(roll)*cos(head)+cos(roll)*sin(pitch)*sin(head), sin(roll)*sin(head)-cos(roll)*sin(pitch)*cos(head), cos(roll)*cos(pitch)];
  Cbn = Cnb';                 %0412 推倒

  %生成陀螺仪输出（角速度）%
  phi = asin(Cnb(2,3))*180/pi;                      %俯仰角
  psi = atan2(Cnb(2,1),Cnb(2,2))*180/pi;            %偏航角
  gamma = atan2(-Cnb(1,3),Cbn(3,3))*180/pi;         %横滚角
% 
%   C =  [-sin(psi) 0 1;
%          sin(gamma)*cos(psi) cos(gamma) 0;
%          cos(gamma)*cos(psi) -sin(gamma) 0];
% %   B = [(phi-phi1);(psi-psi1);(gamma-gamma1)]; 
  B = [(gamma-gamma1);(phi-phi1);(psi-psi1)];
  
  C1 =  [ 0  cos(gamma)  -sin(gamma)*cos(phi);
          1      0          sin(phi);
          0  sin(gamma) cos(gamma)*cos(phi)];    
     
  wibb = C1*atti_rate*pi/180;                        %陀螺仪输出的角速度  
%   wibb1 = C1*atti_rate1 *pi/180;  
%   atti_rate1 = atti_rate;
  
  phi1 = phi;             %俯仰角
  psi1 = psi;             %偏航角
  gamma1 = gamma;         %横滚角
  
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
  
  %计算姿态矩阵%
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
  %计算姿态角%

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

  posiN(1,1)=long*180.0/pi;       %单位：度
  posiN(2,1)=lati*180.0/pi;       %单位：度
  posiN(3,1)=heig;

  TraceData=[TraceData;t,posiN',V,veloN',atti',wibb',acceB',gamma3,phi3,psi3];
                      %   2,3,4  6,7,8  9,10,11 12,13,14     18,   19,  20
end


% fg_1=figure('Name','Trace','NumberTitle','off');
% plot3(TraceData(:,2),TraceData(:,3),TraceData(:,4),'m',TraceData(1,2),TraceData(1,3),TraceData(1,4),'o',TraceData(490,2),TraceData(490,3),TraceData(490,4),'o');
% title('飞行航迹仿真');
% xlabel('经度（^{\circ}）');ylabel('纬度（^{\circ}）');zlabel('高度（m）');
% grid;
% fg_2=figure('Name','trace','NumberTitle','off');
% subplot(3,1,1);plot(TraceData(:,1),TraceData(:,2));title('飞行器经、纬、高度仿真');ylabel('经度（^{\circ}）');grid;
% subplot(3,1,2);plot(TraceData(:,1),TraceData(:,3));ylabel('纬度（^{\circ}）');grid;
% subplot(3,1,3);plot(TraceData(:,1),TraceData(:,4));xlabel('t(sec)');ylabel('高度（m）');grid;
% fg_3=figure('Name','Velo','NumberTitle','off');
% plot(TraceData(:,1),TraceData(:,5));
% title('飞行速度仿真');
% xlabel('t(sec)');ylabel('速度（m/s）');grid;
% fg_4=figure('Name','velo','NumberTitle','off');
% subplot(3,1,1);plot(TraceData(:,1),TraceData(:,6));title('飞行器东、北、天速度仿真');ylabel('东向速度(m/s)');grid;
% subplot(3,1,2);plot(TraceData(:,1),TraceData(:,7));ylabel('北向速度(m/s)');grid;
% subplot(3,1,3);plot(TraceData(:,1),TraceData(:,8));xlabel('t(sec)');ylabel('天向速度(m/s)');grid;

%%%%%%%%% 姿态仿真  %%%%%%%%%%%
fg_5=figure('Name','attitude','NumberTitle','off');
subplot(3,1,1);plot(TraceData(:,1),TraceData(:,9),'b');title('飞行器俯仰、横滚、航向仿真');ylabel('横滚(°)');grid;
hold on  
% plot([0:1:t_stop],gamma3*180/pi,'r');
plot(TraceData(:,1),TraceData(:,18),'r');

subplot(3,1,2);plot(TraceData(:,1),TraceData(:,10));ylabel('俯仰(°)');grid;
hold on  
% plot([0:1:t_stop],phi3*180/pi,'r');
plot(TraceData(:,1),TraceData(:,19),'r');

subplot(3,1,3);plot(TraceData(:,1),TraceData(:,11));xlabel('t(sec)');ylabel('航向(°)');grid;
hold on  
% plot([0:1:t_stop],psi3*180/pi,'r');
plot(TraceData(:,1),TraceData(:,20),'r');
%%%%%%%%%%%%%存储仿真数据%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save trace0.dat TraceData -ASCII;    %存储航迹数据

