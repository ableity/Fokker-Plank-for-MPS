function [t,y] = Neel_relaxation_response(time,H,parameter)
% 2023/1/13 尼尔弛豫
% 2023/1/9修改
% 参考论文Dependence of Brownian and Néel relaxation times on magnetic field strength
% 
% 2023/1/2修改
% 计算方波下的尼尔弛豫响应,该程序只能计算从0到ampulitude的突变响应
% 输入的ampulitude为方波的幅值
% 输出的t为时间，y(:,2)为响应的磁场幅值
% 这里计算的y，其实是粒子磁矩概论分布的勒让德变换系数（y第二个维度，另一个维度是时间）
% 它的第2项（论文中一般是从0下标的第一个）刚好是磁场幅值
% 具体参考下面的论文
% Dependence of Brownian and Néel relaxation times on magnetic field strength

% 定义粒径，在尼尔弛豫中应该是指磁粒子的磁核粒square_wave_Neel_relaxation_response.m径
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
%各向异性常数
K = parameter.K;
beta = Vc/(k*T);
%计算出的勒让德多项式中的时间项中的非磁场项
alpha_constant = m0/(k*T);

%计算出尼尔弛豫中的aplha_K
alpha_K = 2*K*beta;

%尼尔0场弛豫时间
Nrt = Neel_relaxation_zero_field_2(parameter);


% 仿真的时间范围，
tspan = 0:1e-7:time(end);
% 初始值，初始值定义为第一项为0.5，其它项均为0
% 论文中讲了第一项a0为0.5是由概论密度函数的性质得出的（概论积分后总为1）
% 其他项均为0，个人理解为是，在初始条件下（外加磁场为0），磁粒子的角度应该是均匀
% 分布的，这导致了其它项均为0
% 定义N项后截断
N = 50;
% 初值
y0 = zeros(1,N);

% 求解线性方程组
[t,y] = ode15s(@(t,y) odefcn(t,y,Nrt,alpha_constant,alpha_K,time,H,N), tspan, y0);


    function dydt = odefcn(t,y,Nrt,alpha_constant,alpha_K,time,H,N)
        % 定义N项后截断
        alpha_assemble = alpha_constant*alpha(time,H,t);
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


        b = zeros(N,1);
        b(1) = 0.5*alpha_assemble;
        b(2) = 0.5*alpha_K;
        dydt = (A*y+b)/Nrt;
        

    end

    function out = alpha(time,H,t)
        %定义磁场的函数
        %20230212 将外磁场的实测值输入
        %根据实测值H和time序列，生成磁场
        %该方法使用插值，可以生成高于实测数据的采样率数据

        out = interp1(time,H,t);
        out(isnan(out))=0;
%         if t<= 10e-5
%             out = 0;
%         else
%             out = 6e-3;
%         end
    end
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
    
    %阻尼常数
    alpha_dot = parameter.alpha_dot;
    
    
    %磁粒子饱和磁化强度(A/m为单位)
    MS=parameter.MS;
    
    %电子旋磁比
    gamma = parameter.gamma;
    
    t = beta.*(1+alpha_dot^2).*MS/(2.*gamma.*alpha_dot);
    
    
    end

end