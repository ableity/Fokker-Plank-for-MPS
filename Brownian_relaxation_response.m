function [t,y] = Brownian_relaxation_response(time,H,parameter)
% 调整采样率，采样时间等参数，使其符合实际
% 2023/1/9修改
% 参考论文Dependence of Brownian and Néel relaxation times on magnetic field strength
% 
% 2023/1/2修改
% 计算方波下的布朗弛豫响应,该程序只能计算从0到ampulitude的突变响应
% 输入的ampulitude为方波的幅值
% 输出的t为时间，y(:,2)为响应的磁场幅值
% 这里计算的y，其实是粒子磁矩概论分布的勒让德变换系数（y第二个维度，另一个维度是时间）
% 它的第2项（论文中一般是从0下标的第一个）刚好是磁场幅值
% 具体参考下面的论文
% Dependence of Brownian and Néel relaxation times on magnetic field strength

% 定义粒径，在布朗弛豫中应该是指磁粒子的总粒径，即磁核外的部分也要算上
% 20 nm粒径
% Dc = 17e-9;
Dc = parameter.Dc;

%流体力学直径
% Dh = 31e-9;
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

% 仿真的时间范围，
tspan = 0:1e-7:time(end);

% 初始值，初始值定义为第一项为0.5，其它项均为0
% 论文中讲了第一项a0为0.5是由概论密度函数的性质得出的（概论积分后总为1）
% 其他项均为0，个人理解为是，在初始条件下（外加磁场为0），磁粒子的角度应该是均匀
% 分布的，这导致了其它项均为0
% 定义N项后截断
N = 100;
% 初值
y0 = zeros(1,N);

% 求解线性方程组
[t,y] = ode15s(@(t,y) odefcn(t,y,Brt,alpha_constant,time,H,N), tspan, y0);


    function dydt = odefcn(t,y,Brt,alpha_constant,time,H,N)
        % 定义N项后截断
        alpha_assemble = alpha_constant*alpha(time,H,t);
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
        b = zeros(N,1);
        b(1) = 0.5*alpha_assemble;
        
        dydt = (A*y+b)/Brt;
        

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