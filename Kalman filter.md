# Kalman filter

## 前言

最近学习了卡尔曼滤波器，所以整理再上面，之后如果有用到的话再看，这里就不再做前沿的介绍或者别的说明了，如果之后有空的话再做补充吧。

## 公式推导

#### 状态方程

跟选举

##### 真实值：

$$
\tag{1}\label{1}
X_k = A_{k}X_{k-1}+B_ku_{k-1} +w_{k-1}
$$

##### 观测值：

$$
\tag{2}
Z_k = H_{k}X_{k}+V_k
$$

其中  

$w_{k-1}$  为系统的噪声，符合高斯分布 $w_k \sim N(0,Q_k)$ 

$v_k$  为测量噪声，符合高斯分布 $v_k \sim N(0,R_k)$ 

##### 预测值：

$$
\tag{3}
{\mathbf{\hat{X}}_k^-} = {A_k}{\mathbf{\hat{X}}_{k-1}} + {B_k}{u_k}
$$

*注：* 因为无法知道 $w_{k}$  的值，所以我们无法知道真实值 $X_k$ 的具体位置，但忽略噪声影响的预测值**（先验估计值）${\mathbf{\hat{X}}_k^-}$** ，可以通过带入上一次的计算值**（后验估计值）${\hat{X}}_{k-1}$**得出。我们通过预测值，就得到了关于 真实值 $X_k$  的高斯分布预测 真实值 $X_k\sim  N({{\hat{X}}_k^-},P_k^-)$ 

其中 $P_k^-$ 为**先验估计**的协方差矩阵
$$
{P_k^-}&=& E(e_k^-,e_k^{-T})\\
e{_k^-}&=&(X_{k}-{\hat{X}_k^-})
$$


先验估计误差$e{_k^-}$
$$
\begin{aligned}
e{_k^-}&=X_{k}-{\hat{X}_k^-}\\
&= AX_{k-1}+Bu_{k-1}+w_{k-1}-(A\hat{X}_{k-1}+Bu_{k-1})\\
&= A(X_{k-1}-\hat{X}_{k-1})+w_{k-1}\\
&= Ae_{k-1}+w_{k-1}
\end{aligned}
$$
*$e^-_{k-1}$ 为**后验误差***

先验估计误差的协方差矩阵${P_k^-}$：
$$
\begin{aligned}
{P_k^-}&= E(e_k^-,e_k^{-T})\\
&= E[(Ae_{k-1}+w_{k-1})(Ae_{k-1}+w_{k-1})^T]\\
&= E[(Ae_{k-1}+w_{k-1})({e^T_{k-1}}A^T+w^T_{k-1})]\\
&= E[Ae_{k-1}{e^T_{k-1}}A^T+Ae_{k-1}w^T_{k-1}+w_{k-1}{e^T_{k-1}}A^T+w_{k-1}w^T_{k-1}]\\
&= AE(e_{k-1}{e^T_{k-1}})A^T+AE(e_{k-1}w^T_{k-1})+E(w_{k-1}{e^T_{k-1}})A^T+E(w_{k-1}w^T_{k-1})\\
\end{aligned}
$$
其中 由于 $e_{k-1}$ 与 $w_{k-1}$ 相互独立，所以  $AE(e_{k-1}w^T_{k-1})$、 $E(w_{k-1}{e^T_{k-1}})A^T$ 两项为$0$ 。$E(w_{k-1}w^T_{k-1}) = Q_{k-1}$  $E(e_{k-1}{e^T_{k-1}}) = P_{k-1}$ ，所以上式可简化为：
$$
\tag{4}
{P_k^-}= AP_{k-1}A^T+Q_{k-1}
$$
对于观测值 $Z_k = H_{k}X_{k}+V_k$  
$$
\begin{aligned}
Z_k &= H_{k}X_{k}+V_k\\
H_{k}X_{k} &= Z_k-V_k\\
X_{k} &= {H_{k}^{-1}}(Z_k-V_k)\\
\end{aligned}
$$
由于无法建模观测噪声 $V_k$ 的影响，所以
$$
X_{meak} = {H_{k}^{-1}}Z_k
$$


---

### 数据融合

至此，我们就知道了两个数据来源，一个是我们根据 公式 3 进行的预测值，一个是我们根据观测获得的测量值，这两个值都有误差，且误差都服从高斯分布，这就需要我们通过数据融合的思想，来获取更为精确的数据。

#### 后验估计

采用数据融合的思想，新的预测数据 $\hat{X}_{k}$ 由下式获得：
$$
\begin{aligned}
{\hat{X}_{k}}&={\hat{X}_k^-}+{G_k}(X_{meak}-{\hat{X}_k^-})\\
&={\hat{X}_k^-}+{G_k}({H_{k}^{-1}}Z_k-{\hat{X}_k^-})\\\\
&令 G_k = {K_k}{H_{k}}\\
&={\hat{X}_k^-}+{K_k}{H_{k}}({H_{k}^{-1}}Z_k-{\hat{X}_k^-})\\
&={\hat{X}_k^-}+{K_k}(Z_k-{H_{k}}{\hat{X}_k^-})\\
\end{aligned}
$$

$$
\tag{5}
{\hat{X}_{k}}={\hat{X}_k^-}+{K_k}(Z_k-{H_{k}}{\hat{X}_k^-})\\
$$

则**Kalman Gain** $K_k$ 的取值范围为：
$$
\begin{aligned}
&K_k\in[0,H{_k^{-1}}]\\
when:&K_k = 0,&\hat{X}_{k}&={\hat{X}_k^-}\\
when :&K_k = H{_k^{-1}},&\hat{X}_{k}&=X_{meak}
\end{aligned}
$$

这样就变成如下目标，寻找 $K_k$ 值，使得**后验误差** $X_k-\hat{X}_k$ 误差最小。
$$
e_k=X_k-\hat{X}_k
$$
#### 后验误差的协方差矩阵

$$
P_k=E(e,e^T)=
\begin{bmatrix}
\sigma_{e_1}^2 & \sigma_{e_1e_2} \\
\sigma_{e_2e_1} & \sigma_{e_2}^2 \\
\end{bmatrix}
$$
==希望误差最小，既$\hat{X}_k$ 越接近 $X_k$ ，既方差最小，既 $P_k$ 的 $tr(P_k)$ 最小。==
$$
\begin{aligned}
X_k-\hat{X}_k&= X_k-{\hat{X}_k^-}-{K_k}({Z_k}-{H_k}{\hat{X}_k^-})\\
&= X_k-{\hat{X}_k^-}-{K_k}{Z_k}+{K_k}{H_k}{\hat{X}_k^-}\\
&= X_k-{\hat{X}_k^-}-{K_k}({H_k}{X_k}+{V_k})+{K_k}{H_k}{\hat{X}_k^-}\\
&= X_k-{\hat{X}_k^-}-{K_k}{H_k}{X_k}-{K_k}{V_k}+{K_k}{H_k}{\hat{X}_k^-}\\
&= (X_k-{\hat{X}_k^-})-{K_k}{H_k}({X_k}-{\hat{X}_k^-})-{K_k}{V_k}\\
&= (I-{K_k}{H_k})({X_k}-{\hat{X}_k^-})-{K_k}{V_k}\\
&= (I-{K_k}{H_k}){e_k^-}-{K_k}{V_k}\\
\end{aligned}
$$
协方差矩阵：
$$
\begin{aligned}
P_k &= E(e,e^T)\\
&= E(({X_k}-{\hat{X}_k})({X_k}-{\hat{X}_k})^T)\\
&= E((I-{K_k}{H_k}){e_k^-}-{K_k}{V_k})((I-{K_k}{H_k}){e_k^-}-{K_k}{V_k})^T)\\
&= E((I-{K_k}{H_k}){e_k^-}-{K_k}{V_k})({{e_k^-}^T}(I-{K_k}{H_k})^T-{V_k^T}{K_k^T})^T)\\
&= E((I-{K_k}{H_k}){e_k^-}{{e_k^-}^T}(I-{K_k}{H_k})^T-(I-{K_k}{H_k}){e_k^-}{V_k^T}{K_k^T}
-{K_k}{V_k}{{e_k^-}^T}(I-{K_k}{H_k})^T+{K_k}{V_k}{V_k^T}{K_k^T})\\
&= E((I-{K_k}{H_k}){e_k^-}{{e_k^-}^T}(I-{K_k}{H_k})^T) - E((I-{K_k}{H_k}){e_k^-}{V_k^T}{K_k^T})
- E({K_k}{V_k}{{e_k^-}^T}(I-{K_k}{H_k})^T) + E({K_k}{V_k}{V_k^T}{K_k^T})\\
\end{aligned}
$$
其中：
$$
E((I-{K_k}{H_k}){e_k^-}{V_k^T}{K_k^T})
= (I-{K_k}{H_k})E({e_k^-}{V_k^T}){K_k^T}\\
E({K_k}{V_k}{{e_k^-}^T}(I-{K_k}{H_k})^T)
={K_k}E({V_k}{{e_k^-}^T})(I-{K_k}{H_k})^T
$$
在 $E({V_k}{{e_k^-}^T})$ 中 ${V_k}与{{e_k^-}^T}$  相互独立，故上面两个式子都为0.
$$
\begin{aligned}
P_k &= E(e,e^T)\\
&= E((I-{K_k}{H_k}){e_k^-}{{e_k^-}^T}(I-{K_k}{H_k})^T) + E({K_k}{V_k}{V_k^T}{K_k^T})\\
&= (I-{K_k}{H_k}){P_k^-}(I-{K_k}{H_k})^T + {K_k}{R_k}{K_k^T}\\
&= ({P_k^-}-{K_k}{H_k}{P_k^-})(I-{K_k}{H_k})^T + {K_k}{R_k}{K_k^T}\\
&= {P_k^-} - {K_k}{H_k}{P_k^-} - {P_k^-}{H_k^T}{K_k^T} + {K_k}{H_k}{P_k^-}{H_k^T}{K_k^T} + {K_k}{R_k}{K_k^T}\\

\end{aligned}
$$
#### 协方差矩阵的迹

$$
\begin{aligned}
& tr(P_k) &&\\
=& tr({P_k^-}) - tr({K_k}{H_k}{P_k^-}) - tr({P_k^-}{H_k^T}{K_k^T})+ tr({K_k}{H_k}{P_k^-}{H_k^T}{K_k^T}) + tr({K_k}{R_k}{K_k^T})\\
=& tr({P_k^-}) - 2tr({K_k}{H_k}{P_k^-}) + tr({K_k}{H_k}{P_k^-}{H_k^T}{K_k^T}) + tr({K_k}{R_k}{K_k^T})\\
\end{aligned}
\begin{aligned}
&因为：\\
&tr(A) = tr(A^T)\\
&{(({P_k^-}{H_k^T}){K_k^T})^T}={K_k}{H_k}{{P_k^-}^T}\\
&{P_k^-}是对称矩阵{P_k^-}={{P_k^-}^T}
\end{aligned}
$$
求 $K_k$ 的最小值，即 $tr(P_k)$ 对 $K_k$ 求导：
$$
\begin{aligned}
&\frac{d(tr(P_k))}{d(K_k)}=0 - 2({H_k}{{P_k^-}})^T+2{K_k}{H_k}{P_k^-}{H_k^T}+2{K_k}{R_k}&&&\\
\end{aligned}
\begin{aligned}
因为：&\\
&\frac{d(tr(AB))}{d(A)}=B^T\\
&\frac{d(tr(AB{A^T}))}{d(A)}=2AB\\
\end{aligned}\\
$$
令 $\frac{d(tr(P_k))}{d(K_k)}=0$ 
$$
\begin{aligned}
-2({H_k}{{P_k^-}})^T+2{K_k}{H_k}{P_k^-}{H_k^T}+2{K_k}{R_k}&=0&&&&&&&&&&\\
{K_k}{H_k}{P_k^-}{H_k^T}+{K_k}{R_k}&=({H_k}{{P_k^-}})^T\\
{K_k}({H_k}{P_k^-}{H_k^T}+{R_k})&={{P_k^-}^T}{H_k^T}\\
\end{aligned}
$$

#### 卡尔曼增益

$$
\tag{6}
{K_k}=\frac{{{P_k^-}^T}{H_k^T}}{{H_k}{P_k^-}{H_k^T}+{R_k}}\\
$$

同时：
$$
\begin{aligned}
P_k  &= {P_k^-} - {K_k}{H_k}{P_k^-} - {P_k^-}{H_k^T}{K_k^T} + {K_k}{H_k}{P_k^-}{H_k^T}{K_k^T} + {K_k}{R_k}{K_k^T}\\
&= {P_k^-} - {K_k}{H_k}{P_k^-} - {P_k^-}{H_k^T}{K_k^T} + {K_k}({H_k}{P_k^-}{H_k^T} + {K_k}{R_k}){K_k^T}\\
&= {P_k^-} - {K_k}{H_k}{P_k^-} - {P_k^-}{H_k^T}{K_k^T} + {P_k^-}{H_k^T}{K_k^T}\\
&= {P_k^-} - {K_k}{H_k}{P_k^-}\\
&= (I - {K_k}{H_k}){P_k^-}\\
\end{aligned}
$$

$$
\tag{7}
P_k = (I - {K_k}{H_k}){P_k^-}\\
$$



至此，公式推导结束。

### 卡尔曼的使用

1. 计算先验估计值
   $$
   {\mathbf{\hat{X}}_k^-} = {A_k}{\mathbf{\hat{X}}_{k-1}} + {B_k}{u_k}
   $$
   
2. 更新先验误差的协方差矩阵
   $$
   {P_k^-}= AP_{k-1}A^T+Q_{k-1}
   $$
   
3. 计算Kalman Gain
   $$
   {K_k}=\frac{{{P_k^-}^T}{H_k^T}}{{H_k}{P_k^-}{H_k^T}+{R_k}}\\
   $$
   
4. 计算后验估计
   $$
   {\hat{X}_{k}}={\hat{X}_k^-}+{K_k}(Z_k-{H_{k}}{\hat{X}_k^-})
   $$
   
5. 更新误差协方差矩阵

$$
P_k = (I - {K_k}{H_k}){P_k^-}
$$

### 引用

【1】[【卡尔曼滤波器】3_卡尔曼增益超详细数学推导 ～全网最完整_哔哩哔哩_bilibili](https://www.bilibili.com/video/BV1hC4y1b7K7/?spm_id_from=333.788&vd_source=87f4f47efa181a16a403b1e44adc1317)

【2】[通过简单直观的推导理解卡尔曼基础)Understanding the Basis of the Kalman Filter Via a Simple and Intuitive Derivation_tuszhangs的博客-CSDN博客](https://blog.csdn.net/weixin_39675633/article/details/103841638)

【3】[Understanding the Basis of the Kalman Filter Via a Simple and Intuitive Derivation (Lecture Note | IEEE Journals & Magazine | IEEE Xplore )](https://ieeexplore.ieee.org/document/6279585)

【4】[轻松理解卡尔曼滤波 - 知乎](https://zhuanlan.zhihu.com/p/444977764)

【5】[How a Kalman filter works, in pictures | Bzarg](https://www.bzarg.com/p/how-a-kalman-filter-works-in-pictures/)

【6】[kalman滤波理解一：理论框架_还差得远呢的博客-CSDN博客](https://blog.csdn.net/u011362822/article/details/95894836)

【7】[probability - Is the product of two Gaussian random variables also a Gaussian? - Mathematics Stack Exchange](https://math.stackexchange.com/questions/101062/is-the-product-of-two-gaussian-random-variables-also-a-gaussian)





