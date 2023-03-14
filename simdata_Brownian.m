function [t,u_sim]=simdata_Brownian(H,time,parameter)

% 20230301李蕾
% 将参数全部移入该脚本
% 保存于parameter

% 20230215 李蕾
% 根据输入的磁场进行仿真


%[t,M2] = square_wave_Brownian_relaxation_response_2(time,H*1e-3,parameter);
[t,M2] = Brownian_relaxation_response(time,H,parameter);
u_sim = diff(M2(:,1));
t = t(2:end);
u_sim = u_sim(1:10:end);
t = t(1:10:end);



end

