function relaxtion_time =calc_Brownian_relaxation_time(B,parameter)
%2023 1 13 布朗弛豫时间和粒径、磁场强度的关系（特征值计算）
%磁核
Dc = parameter.Dc;

%流体力学直径
Dh = parameter.Dh;

%粒子体积
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

%布朗0场弛豫时间
Brt = Brownian_relaxation_zero_field(parameter);
if B==0
    relaxtion_time = Brt;
    return
end

% 定义N项后截断
N = 100;

alpha_assemble = alpha_constant*B;
% 这里参考ode15s的说明，定义了一个函数矩阵，该矩阵表示的线性方程组为
% y(1)'=矩阵第一项
% y(2)'=矩阵第二项
% 直到n
% 通过30项截断（即n=30）阻止无限向后递归，可以从结果看出n变大时值会越来越小
A = zeros(N,N);
A(1,1) = -1;
A(1,2) = -0.2*alpha_assemble;
A(N,N) = -N*(N+1)/2;
A(N,N-1) = N*(N+1)*alpha_assemble/(2*(2*N-1));
for n = 2:N-1
    A(n,n)  = -n*(n+1)/2;
    A(n,n+1) = -n*(n+1)*alpha_assemble/(2*(2*n+3));
    A(n,n-1) = n*(n+1)*alpha_assemble/(2*(2*n-1));
end

[~,E] = eig(A);
eigenvalue = zeros(N,1);
for i=1:N
    eigenvalue(i) = E(i,i);
end
eigenvalue = real(eigenvalue);
max_eigenvalue = max(max(eigenvalue));

relaxtion_time = -Brt/max_eigenvalue;


    function t = Brownian_relaxation_zero_field(parameter)
        %202
        %零场布朗弛豫时间,用在后续参数调用中，返回值是弛豫时间
        %输入参数为粒子直径(流体力学直径)

        %液体粘度，单位Pa s
        % n = 1.0049*1e-3;
        n = parameter.n;

        %玻尔兹曼常数
        k = parameter.k;

        %温度
        T = parameter.T;

        %粒子体积
        Vh =1/6*pi*parameter.Dh.^3;

        %零场布朗弛豫时间
        t = 3*n*Vh/(k*T);
    end

end





