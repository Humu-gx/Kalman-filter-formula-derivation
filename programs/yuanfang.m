clear all;
clc;
%% 参数设定
kx = .01; ky = .05; % 阻尼系数
g = 9.8; % 重力
t = 15; % 仿真时间
Ts = 0.1; % 采样周期
len = fix(t/Ts); % 仿真步数
dim_observe=2;
dim_system=4;
x_noise =0.3; y_noise =0.3;
r_noise = 8; a_noise = 0.1;
Qk = diag([0; x_noise; 0; y_noise])^2;
Rk = diag([r_noise; a_noise])^2;
wx=x_noise*randn(1,len)*Ts;
wy=y_noise*randn(1,len)*Ts; % 生成干扰力噪声
X = zeros(len,4); %真实状态
Xk=zeros(4,1);
Pk=zeros(dim_system,dim_system);
X(1,:) = [0, 66, 409, 0]; % 状态模拟的初值
X_est=X;%估计状态
x_hat = [0,65,405,0]; %估计状态初值
vr=r_noise*randn(1,len);
va=a_noise*randn(1,len); % 生成量测噪声
%% 模拟水平抛物真实轨迹
for k=2:len
x = X(k-1,1); vx = X(k-1,2); y = X(k-1,3); vy = X(k-1,4);
x = x + vx*Ts;
vx = vx -kx*vx^2*Ts+wx(k-1);
y = y + vy*Ts;
vy = vy + (ky*vy^2-g)*Ts+wy(k-1);
X(k,:) = [x, vx, y, vy];
end
%% 构造量测量
for k=1:len
r = sqrt(X(k,1)^2+X(k,3)^2) + vr(k);
a = atan(X(k,1)/X(k,3))*57.3 + va(k);
Z(k,:) = [r, a];
end
%% 第一次 ekf 滤波
fX=projectile(x_hat,kx,ky,g,Ts);
A=JacobianF(x_hat,kx,ky,g);
H=JacobianH(x_hat);
hX=measurement(fX,Ts);
[Xk,Pk,Kk]=EKF(A,Qk,fX,Pk,H,Rk,Z(1,:),hX);
X_est(1,:) = Xk';
%% ekf 滤波
for k=2:len
    A=JacobianF(X_est(k-1,:),kx,ky,g);
    H=JacobianH(X_est(k-1,:));
    fX=projectile(X_est(k-1,:),kx,ky,g,Ts);
    hX=measurement(fX,Ts);
    [Xk,Pk,Kk]=EKF(A,Qk,fX,Pk,H,Rk,Z(k-1,:),hX);
    X_est(k,:) = Xk';
end
figure, hold on, grid on;
plot(X(:,1),X(:,3),'-k');
plot(Z(:,1).*sin(Z(:,2)*pi/180), Z(:,1).*cos(Z(:,2)*pi/180),'*b');
plot(X_est(:,1),X_est(:,3), 'r');
xlabel('X');
ylabel('Y');
title('EKF simulation');
legend('real', 'measurement', 'ekf estimated');
axis([-5,250,200,425]);
%% ***************************子函数*********************************
function fX=projectile(X,kx,ky,g,Ts)
    x = X(1); vx = X(2); y = X(3); vy = X(4); Ts = 0.1;
    x1= x + vx*Ts;
    vx1= vx -kx*vx^2*Ts;
    y1= y + vy*Ts;
    vy1= vy + (ky*vy^2-g)*Ts;
    fX= [x1;vx1;y1;vy1];
end
function F=JacobianF(X,kx,ky,g)%系统状态雅可比函数
    vx=X(2);
    vy=X(4);
    Ts = 0.1;
    F=zeros(4,4);
    F(1,1) = 1;
    F(1,2) = Ts;
    F(2,2) = 1-2*kx*vx*Ts;
    F(3,3) = 1;
    F(3,4) = Ts;
    F(4,4) = 1+2*ky*vy*Ts;
end
function H=JacobianH(X)%量测雅可比函数
    x=X(1);
    y=X(3);
    Ts = 0.1;
    H=zeros(2,4);
    r = sqrt(x^2+y^2);
    xy2 = 1+(x/y)^2;
    H(1,1) = x/r;
    H(1,3) = y/r;
    H(2,1) = (1/y)/xy2;
    H(2,3) = (-x/y^2)/xy2;
end
function hX=measurement(fX,Ts)
    x = fX(1);
    y = fX(3);
    Ts = 0.1;
    r = sqrt(x^2+y^2) ;
    a = atan(x/y)*57.3 ;
    hX = [r;a];
end
function [Xk,Pk,Kk]=EKF(EKF_A,EKF_Qk,EKF_x_forecast,EKF_P_last,EKF_H,EKF_Rk,EKF_Zk,EKF_h_x)
    P_forecast=EKF_A*EKF_P_last*EKF_A'+EKF_Qk;
    Kk = P_forecast*EKF_H'*(EKF_H*P_forecast*EKF_H'+EKF_Rk)^-1; %计算增益
    Xk = EKF_x_forecast+Kk*(EKF_Zk'-EKF_h_x); %校正
    Pk = (eye(4)-Kk*EKF_H)*P_forecast;
end