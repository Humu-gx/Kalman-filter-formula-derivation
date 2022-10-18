%%% DrCan 状态空间的的例子

%% 数据初始化
% 确定初始参数
t = 100;           % 运算时间
Delta_T = 0.1;        % 采样间隔
I = [1 0 0; 0 1 0; 0 0 1];
A = [1 Delta_T 0; 0 1 Delta_T; 0 0 1];
H = [1 0 0; 0 1 0; 0 0 1];

% 引入噪声 协方差矩阵 
Q = [1 0 0; 0 1 0; 0 0 0.04];
R = [3 0 0; 0 0.7 0; 0 0 0.5];

% 设置初始量
Pk1 = [1 0 0; 0 1 0; 0 0 1];
Xk1 = [ 0; 0; 0.5]; % 初始距离，速度，加速度

% 设置容器
X = zeros(t/Delta_T,3);
Y = zeros(t/Delta_T,3);
Z = zeros(t/Delta_T,3);

%% 迭代计算
rng(10);            %设置随机数种子
for i = 1:t/Delta_T
    % 噪声生成
    Wk1 = Q.^0.5 * randn(3,1);
    Vk = R.^0.5 * randn(3,1);
    
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

% tiledlayout('flow')
% nexttile
figure(1)
hold on
plot(X(:,1),'r','LineWidth',1)
plot(Y(:,1),'g','LineWidth',1)
plot(Z(:,1),'b','LineWidth',1)

% nexttile
figure(2)
hold on
plot(X(:,2),'r','LineWidth',1)
plot(Y(:,2),'g','LineWidth',1)
plot(Z(:,2),'b','LineWidth',1)

% nexttile
figure(3)
hold on
plot(X(:,3),'r','LineWidth',1)
plot(Y(:,3),'g','LineWidth',1)
plot(Z(:,3),'b','LineWidth',1)
