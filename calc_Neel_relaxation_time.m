function relaxtion_time =calc_Neel_relaxation_time(B,parameter)
%2023 1 13 布朗弛豫时间和粒径、磁场强度的关系（特征值计算）
%粒子体积
%磁核
Dc = parameter.Dc;

Vc =1/6*pi*Dc.^3;

%磁粒子饱和磁化强度(A/m为单位)
MS=parameter.MS;

%玻尔兹曼常数
k = parameter.k;

%温度
T = parameter.T;



%磁矩
m0 = MS*Vc;

%计算出的勒让德多项式中的时间项中的非磁场项
alpha_constant = m0/(k*T);

%各向异性常数
K = parameter.K;
beta = Vc/(k*T);

%计算出尼尔弛豫中的aplha_K
alpha_K = 2*K*beta;
%尼尔0弛豫时间
%注意，这里的弛豫时间并不是尼尔零场弛豫时间，它的零场弛豫时间有更复杂的形式
Nrt = Neel_relaxation_zero_field_2(parameter);

% if B==0
%     relaxtion_time = Nrt;
%     return
% end

% 定义N项后截断
N = 50;

alpha_assemble = alpha_constant*B;
% 这里参考ode15s的说明，定义了一个函数矩阵，该矩阵表示的线性方程组为
% y(1)'=矩阵第一项
% y(2)'=矩阵第二项
% 直到n
% 通过30项截断（即n=30）阻止无限向后递归，可以从结果看出n变大时值会越来越小
A = zeros(N,N);
n=1;
A(n,n)  = n*(n+1)/2*(-1+n*alpha_K/((2*n-1)*(2*n+1)) ...
    -(n+1)*alpha_K/((2*n+1)*(2*n+3)));
A(n,n+1) = -n*(n+1)*alpha_assemble/(2*(2*n+3));
A(n,n+2) = -n*(n+1)*(n+2)/(2*(2*n+3)*(2*n+5))*alpha_K;

n=2;
A(n,n)  = n*(n+1)/2*(-1+n*alpha_K/((2*n-1)*(2*n+1)) ...
    -(n+1)*alpha_K/((2*n+1)*(2*n+3)));
A(n,n+1) = -n*(n+1)*alpha_assemble/(2*(2*n+3));
A(n,n-1) = n*(n+1)*alpha_assemble/(2*(2*n-1));
A(n,n+2) = -n*(n+1)*(n+2)/(2*(2*n+3)*(2*n+5))*alpha_K;

n = N-1;
A(n,n)  = n*(n+1)/2*(-1+n*alpha_K/((2*n-1)*(2*n+1)) ...
    -(n+1)*alpha_K/((2*n+1)*(2*n+3)));
A(n,n+1) = -n*(n+1)*alpha_assemble/(2*(2*n+3));
A(n,n-1) = n*(n+1)*alpha_assemble/(2*(2*n-1));
A(n,n-2) = n*(n+1)*(n-1)/(2*(2*n-3)*(2*n-1))*alpha_K;

n = N;
A(n,n)  = n*(n+1)/2*(-1+n*alpha_K/((2*n-1)*(2*n+1)) ...
    -(n+1)*alpha_K/((2*n+1)*(2*n+3)));
A(n,n-1) = n*(n+1)*alpha_assemble/(2*(2*n-1));
A(n,n-2) = n*(n+1)*(n-1)/(2*(2*n-3)*(2*n-1))*alpha_K;

for n = 3:N-2
    A(n,n)  = n*(n+1)/2*(-1+n*alpha_K/((2*n-1)*(2*n+1)) ...
        -(n+1)*alpha_K/((2*n+1)*(2*n+3)));
    A(n,n+1) = -n*(n+1)*alpha_assemble/(2*(2*n+3));
    A(n,n-1) = n*(n+1)*alpha_assemble/(2*(2*n-1));
    A(n,n+2) = -n*(n+1)*(n+2)/(2*(2*n+3)*(2*n+5))*alpha_K;
    A(n,n-2) = n*(n+1)*(n-1)/(2*(2*n-3)*(2*n-1))*alpha_K;
end

[~,E] = eig(A);
eigenvalue = zeros(N,1);
for i=1:N
    eigenvalue(i) = E(i,i);
end
eigenvalue = real(eigenvalue);
max_eigenvalue = max(max(eigenvalue));

relaxtion_time = -Nrt/max_eigenvalue;



    function t = Neel_relaxation_zero_field_2(parameter)
    %零场尼尔弛豫时间，输入参数为粒子直径
    
    %磁粒子直径
    % D = D*1e-9;
    Dc = parameter.Dc;

    %计算磁粒子体积
    V=1/6*pi*Dc.^3;
    
    %玻尔兹曼常数
    k = parameter.k;
    %温度
    T = parameter.T;
    
    beta = V/(k*T);
    tt = beta;
    %阻尼常数
    alpha_dot = parameter.alpha_dot;
    
    
    %磁粒子饱和磁化强度(A/m为单位)
    MS=parameter.MS;
    
    %电子旋磁比
    gamma = parameter.gamma;
    
    t = beta.*(1+alpha_dot^2).*MS/(2.*gamma.*alpha_dot);
    

    end

end





