Z = 3*(rand(100,1)*2-1)+50;
X_hat = zeros(100,1);
G_K   = zeros(100,1);
e     = zeros(100,1);
X_hat(1) = 40;
e(1)     = 7;
G_K(1)   = 0;
%% 运算部分
for k = 2:100
    G_K(k) = e(k-1)/(e(k-1)+3);
    X_hat(k) = X_hat(k-1)+G_K(k)*(Z(k)-X_hat(k-1));
    e(k) = (1-G_K(k))*e(k-1);
end
figure(1);
plot(Z);
hold on
plot(X_hat)
legend('测量值','估计值');