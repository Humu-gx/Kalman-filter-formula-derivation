%%% DrCan 状态空间的的例子

%% 数据初始化
% 确定初始参数
t = 1000;            % 运算时间
Delta_T = 1;      % 采样间隔
I = [1 0; 0 1]; 
A = [1 Delta_T; 0 1];
H = [1 0; 0 1];

% 引入噪声 协方差矩阵 
Q = [5 0; 0 4];
R = [1 0; 0 1];

% 设置初始量
Pk1 = [1 0; 0 1];
Xk1 = [ 0; 1]; % 初始距离和速度

% 设置容器
X = zeros(t/Delta_T,2);
Y = zeros(t/Delta_T,2);
Z = zeros(t/Delta_T,2);
%% 迭代计算
rng(10);            %设置随机数种子
for i = 1:t/Delta_T
    % 噪声生成
    Wk1 = Q.^0.5 * randn(2,1);
    Vk = R.^0.5 * randn(2,1);
    
    % 生成实际 真实数据 与 测量数据
    Xk = A * Xk1 + Wk1;
    Zk = H * Xk + Vk;
    
    % 预测
    Xk_p = A * Xk1;
    Pk_p = A * Pk1 * A' + Q;
    
    % 矫正
    Kk = (Pk_p * H')/(H * Pk_p * H' + R);
    Xk_head = Xk_p + Kk * (Zk - H * Xk_p);
    Pk = (I - Kk * H) * Pk_p;
    
    % 更新
    Pk1 = Pk;
    Xk1 = Xk_head;
    
    X(i,:) = Xk_head';
    Y(i,:) = Xk';
    Z(i,:) = Zk';
end

Pk1

tiledlayout('flow')
nexttile
hold on
plot(X(:,1),'r','LineWidth',1)
plot(Y(:,1),'g','LineWidth',1)
plot(Z(:,1),'b','LineWidth',1)
legend('真实值','测量值','估测值');
title('距离测量值')
xlabel('时间')
ylabel('距离')

nexttile
hold on
plot(X(:,2),'r','LineWidth',1.5)
plot(Y(:,2),'g','LineWidth',1.5)
plot(Z(:,2),'b','LineWidth',1.5)
legend('真实值','测量值','估测值');
title('速度测量值')
xlabel('时间')
ylabel('速度')












