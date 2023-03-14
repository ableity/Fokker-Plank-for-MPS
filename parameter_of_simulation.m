function parameter = parameter_of_simulation(options)

% 2023 03 06
% 李蕾
% 参数集合，如果是使用的synomag-D的脚本，则粒径和K等参数会在脚本中充值
    arguments
        options.Dc double = 20*1e-9
        options.Dh double = 20*1e-9
        options.MS double = 474000
        options.k double = 1.380649e-23
        options.T double = 300
        options.n double = 1.0049e-3
        options.K double = 20000
        options.alpha_dot double = 0.1
        options.gamma double = 1.75*1e11
    end



%磁核粒径
parameter.Dc = options.Dc;
%水合粒径
parameter.Dh = options.Dh;
%磁粒子饱和磁化强度(A/m为单位)
parameter.MS = options.MS;
%玻尔兹曼常数
parameter.k = options.k;
%温度
parameter.T = options.T;
%液体粘度，单位Pa s
parameter.n = options.n;
% parameter.n = 1.0049e-3;

% 尼尔弛豫特有参数
%各向异性常数
parameter.K = options.K;

%阻尼常数
parameter.alpha_dot = options.alpha_dot;

%电子旋磁比
parameter.gamma = options.gamma;

end