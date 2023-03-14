clc
clear
close all
%% 仿真synomag-D
t = 0:1e-6:1/(25e3)*10;
H = 10*1e-3*sin(2*pi*25e3*t);
parameter = parameter_of_simulation();

[t,ub,un] = simdata_synomag(H,t,parameter);


figure
subplot(2,1,1)
plot(t,ub,LineWidth=2)
title("synomag-D布朗弛豫响应",fontsize = 30)
subplot(2,1,2)
plot(t,un,LineWidth=2)
title("synomag-D尼尔弛豫响应",fontsize = 30)
clc
clear


%% 仿真布朗弛豫和尼尔弛豫信号
t = 0:1e-6:1/(25e3)*10;
H = 10*1e-3*sin(2*pi*25e3*t);
parameter = parameter_of_simulation();

[t_b,u_brownian] = simdata_Brownian(H,t,parameter);
[t_n,u_neel] = simdata_Neel(H,t,parameter);

figure
subplot(2,1,1)
plot(t_n,u_brownian,LineWidth=2)
title("布朗弛豫响应",fontsize = 30)
subplot(2,1,2)
plot(t_n,u_neel,LineWidth=2)
title("尼尔弛豫响应",fontsize = 30)
clc
clear




%% 弛豫时间随磁场的变化曲线

% 定义参数，默认参数见函数parameter_of_simulation
% 根据需要可以自行对参数修改，修改方式见下语句
% parameter_of_simulation(Dh=20e-9,Dc=20e-9,n=1.01e-3)
parameter = parameter_of_simulation();

%存储布朗和尼尔弛豫时间的向量
B = zeros(1,100);
N = zeros(1,100);
%外磁场从0到99mT
for H = 0:99
    B(H+1) = calc_Brownian_relaxation_time(H*1e-3,parameter);
    N(H+1) = calc_Neel_relaxation_time(H*1e-3,parameter);
end

%绘图
figure
plot(0:99,log(B)/log(10),LineWidth=2)
hold on
plot(0:99,log(N)/log(10),LineWidth=2)

ylabel("log(t)",fontsize=20)
xlabel("磁场,mT",fontsize=20)
legend("布朗弛豫","尼尔弛豫",fontsize=20)
title("弛豫时间随磁场变换图",fontsize=30)


clc
clear

%% 弛豫时间随黏度的变化
n = 0.9:0.1:5;
n=n*1e-3;

relaxtime = zeros(10,42);
for H=1:10
    for i =1:42
        parameter = parameter_of_simulation(n=n(i));
        relaxtime(H,i)=calc_Brownian_relaxation_time(H*1e-3,parameter);
    end
end
relaxtime = relaxtime*1e6;

figure
plot(n,relaxtime',LineWidth=3)
xlabel("黏度 Pa.s",fontsize=20)
ylabel("弛豫时间 µs",fontsize=20)
legend(num2str((1:10)')+"mT",fontsize=20)
title("不同磁场下布朗弛豫时间和黏度的关系",fontsize=25)