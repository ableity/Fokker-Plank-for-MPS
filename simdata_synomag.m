function [t,ub,un]=simdata_synomag(H,time,parameter)

% 20230306李蕾
% 根据磁粒子的粒径分布
Dc = 25:33;
weight = normpdf(Dc,29,4);
weight = weight/sum(sum(weight));
K = (0.150*Dc+9.5)*1e3;

ub = 0;
un = 0;
for i = 1:9
    parameter.Dc = Dc(i)*1e-9;
    parameter.Dh = (Dc(i)+16.9)*1e-9;
    parameter.K = K(i);
    [~,u_B] = simdata_Brownian(H,time,parameter);
    [t,u_N] = simdata_Neel(H,time,parameter);


    ub = ub + u_B*weight(i);
    un = un + u_N*weight(i);
end





end