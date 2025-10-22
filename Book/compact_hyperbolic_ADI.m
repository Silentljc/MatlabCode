%紧双曲ADI

clear;
clc;
tic;
M = [10, 20, 40];  % 空间网格数量
N = [100, 400, 1600];  % 时间网格数量
error_inf = zeros(1, length(M));
Norm = zeros(1, length(M)-1);

for p = 1:length(M)
    h = 1 / M(p);  % 空间步长
    tau = 1 / N(p);  % 时间步长
    x = 0:h:1;
    y = 0:h:1;
    t = 0:tau:1;
    
    Numerical = zeros(M(p)+1, M(p)+1, N(p)+1);  % u
    numerical = zeros(M(p)+1, M(p)-1);  % u*

    r = tau^2 / (h^2);
    alpha = 1/12 - r/2;
    beta = 10/12 + r;
    
    % 构造三对角矩阵系数
    a = alpha * ones(1, M(p)-2);
    b = beta * ones(1, M(p)-1);
    c = alpha * ones(1, M(p)-2);

    % 设置初值
    for i = 1:M(p)+1
        for j = 1:M(p)+1
            Numerical(i,j,1) = exp(1/2*(x(i)+y(j)));  % 初值
        end
    end
    
    % 设置边界条件
    for j = 1:M(p)+1
        for k = 1:N(p)+1
            Numerical(1,j,k) = exp(1/2*y(j)-t(k));  % u(0,y,t)
        end
    end
    
    for j = 1:M(p)+1
        for k = 1:N(p)+1
            Numerical(M(p)+1,j,k) = exp(1/2*(1+y(j))-t(k));  % u(1,y,t)
        end
    end
    
    for i = 1:M(p)+1
        for k = 1:N(p)+1
            Numerical(i,1,k) = exp(1/2*x(i)-t(k));  % u(x,0,t)
        end
    end
    
    for i = 1:M(p)+1
        for k = 1:N(p)+1
            Numerical(i,M(p)+1,k) = exp(1/2*(1+x(i))-t(k));  % u(x,1,t)
        end
    end
    
    % 定义函数句柄
    f_func = @(x,y,t) 1/2*exp(1/2*(x+y)-t);
    fun = @(x,y,t) exp(1/2*(x+y)-t);
    psi = @(x,y) -exp(1/2*(x+y));
    rho = @(x,y) -exp(1/2*(x+y));
    
    % 计算精确解
    Accurate = zeros(M(p)+1, M(p)+1, N(p)+1);
    for i = 1:M(p)+1
        for j = 1:M(p)+1
            for k = 1:N(p)+1
                Accurate(i,j,k) = fun(x(i), y(j), t(k));
            end
        end
    end

    % 计算Numerical（x,y,1）=u(i,j,1)
    for j = 1:M(p)-1  % 固定j
        numerical(1,j) = alpha * Numerical(1,j,2) + ...
            beta * Numerical(1,j+1,2) + ...
            alpha * Numerical(1,j+2,2);  % u*0j1
        numerical(M(p)+1,j) = alpha * Numerical(M(p)+1,j,2) + ...
            beta * Numerical(M(p)+1,j+1,2) + ...
            alpha * Numerical(M(p)+1,j+2,2);  % u*mj1
        
        % 循环生成右端列向量
        numerical_right_vector = zeros(M(p)-1, 1);
        for i = 1:M(p)-1
            numerical_right_vector(i) = (1/144) * ( ...
                Numerical(i,j,1) + tau*psi(x(i),y(j)) - tau^3*rho(x(i),y(j))/3 + ...
                tau^2*f_func(x(i),y(j),t(2))/2 + ...
                10*(Numerical(i+1,j,1) + tau*psi(x(i+1),y(j)) - tau^3*rho(x(i+1),y(j))/3 + ...
                tau^2*f_func(x(i+1),y(j),t(2))/2) + ...
                Numerical(i+2,j,1) + tau*psi(x(i+2),y(j)) - tau^3*rho(x(i+2),y(j))/3 + ...
                tau^2*f_func(x(i+2),y(j),t(2))/2) + ...
                (1/144) * 10 * ( ...
                Numerical(i,j+1,1) + tau*psi(x(i),y(j+1)) - tau^3*rho(x(i),y(j+1))/3 + ...
                tau^2*f_func(x(i),y(j+1),t(2))/2 + ...
                10*(Numerical(i+1,j+1,1) + tau*psi(x(i+1),y(j+1)) - tau^3*rho(x(i+1),y(j+1))/3 + ...
                tau^2*f_func(x(i+1),y(j+1),t(2))/2) + ...
                Numerical(i+2,j+1,1) + tau*psi(x(i+2),y(j+1)) - tau^3*rho(x(i+2),y(j+1))/3 + ...
                tau^2*f_func(x(i+2),y(j+1),t(2))/2) + ...
                (1/144) * ( ...
                Numerical(i,j+2,1) + tau*psi(x(i),y(j+2)) - tau^3*rho(x(i),y(j+2))/3 + ...
                tau^2*f_func(x(i),y(j+2),t(2))/2 + ...
                10*(Numerical(i+1,j+2,1) + tau*psi(x(i+1),y(j+2)) - tau^3*rho(x(i+1),y(j+2))/3 + ...
                tau^2*f_func(x(i+1),y(j+2),t(2))/2) + ...
                Numerical(i+2,j+2,1) + tau*psi(x(i+2),y(j+2)) - tau^3*rho(x(i+2),y(j+2))/3 + ...
                tau^2*f_func(x(i+2),y(j+2),t(2))/2);
        end
        
        % 添加边界贡献
        numerical_right_vector(1) = numerical_right_vector(1) - alpha * numerical(1,j);
        numerical_right_vector(M(p)-1) = numerical_right_vector(M(p)-1) - alpha * numerical(M(p)+1,j);
        
        numerical(2:M(p),j) = chase(a, b, c, numerical_right_vector);
    end
    
    for i = 1:M(p)-1  % 固定i
        Numerical_right_vector = zeros(M(p)-1, 1);
        for j = 1:M(p)-1
            Numerical_right_vector(j) = numerical(i+1,j);
        end
        
        Numerical_right_vector(1) = Numerical_right_vector(1) - alpha * Numerical(i+1,1,2);
        Numerical_right_vector(M(p)-1) = Numerical_right_vector(M(p)-1) - alpha * Numerical(i+1,M(p)+1,2);
        
        Numerical(i+1,2:M(p),2) = chase(a, b, c, Numerical_right_vector);
    end

    % 核心部分
    for k = 2:N(p)
        for j = 1:M(p)-1  % 固定j
            numerical(1,j) = alpha * Numerical(1,j,k+1) + ...
                beta * Numerical(1,j+1,k+1) + ...
                alpha * Numerical(1,j+2,k+1);  % u*0j
            numerical(M(p)+1,j) = alpha * Numerical(M(p)+1,j,k+1) + ...
                beta * Numerical(M(p)+1,j+1,k+1) + ...
                alpha * Numerical(M(p)+1,j+2,k+1);  % u*mj
            
            % 循环生成右端列向量
            numerical_right_vector = zeros(M(p)-1, 1);
            for i = 1:M(p)-1
                numerical_right_vector(i) = (r/24) * ( ...
                    Numerical(i,j,k-1) - 2*Numerical(i+1,j,k-1) + Numerical(i+2,j,k-1) + ...
                    10*(Numerical(i,j+1,k-1) - 2*Numerical(i+1,j+1,k-1) + Numerical(i+2,j+1,k-1)) + ...
                    Numerical(i,j+2,k-1) - 2*Numerical(i+1,j+2,k-1) + Numerical(i+2,j+2,k-1)) + ...
                    (r/24) * ( ...
                    Numerical(i,j,k-1) - 2*Numerical(i,j+1,k-1) + Numerical(i,j+2,k-1) + ...
                    10*(Numerical(i+1,j,k-1) - 2*Numerical(i+1,j+1,k-1) + Numerical(i+1,j+2,k-1)) + ...
                    Numerical(i+2,j,k-1) - 2*Numerical(i+2,j+1,k-1) + Numerical(i+2,j+2,k-1)) - ...
                    (1/144) * ( ...
                    Numerical(i,j,k-1) + 10*Numerical(i+1,j,k-1) + Numerical(i+2,j,k-1) + ...
                    10*(Numerical(i,j+1,k-1) + 10*Numerical(i+1,j+1,k-1) + Numerical(i+2,j+1,k-1)) + ...
                    Numerical(i,j+2,k-1) + 10*Numerical(i+1,j+2,k-1) + Numerical(i+2,j+2,k-1)) + ...
                    (2/144) * ( ...
                    Numerical(i,j,k) + 10*Numerical(i+1,j,k) + Numerical(i+2,j,k) + ...
                    10*(Numerical(i,j+1,k) + 10*Numerical(i+1,j+1,k) + Numerical(i+2,j+1,k)) + ...
                    Numerical(i,j+2,k) + 10*Numerical(i+1,j+2,k) + Numerical(i+2,j+2,k)) + ...
                    (tau^2/144) * ( ...
                    f_func(x(i),y(j),t(k)) + 10*f_func(x(i+1),y(j),t(k)) + f_func(x(i+2),y(j),t(k)) + ...
                    10*(f_func(x(i),y(j+1),t(k)) + 10*f_func(x(i+1),y(j+1),t(k)) + f_func(x(i+2),y(j+1),t(k))) + ...
                    f_func(x(i),y(j+2),t(k)) + 10*f_func(x(i+1),y(j+2),t(k)) + f_func(x(i+2),y(j+2),t(k)));
            end
            
            % 添加边界贡献
            numerical_right_vector(1) = numerical_right_vector(1) - alpha * numerical(1,j);
            numerical_right_vector(M(p)-1) = numerical_right_vector(M(p)-1) - alpha * numerical(M(p)+1,j);
            
            numerical(2:M(p),j) = chase(a, b, c, numerical_right_vector);
        end
        
        for i = 1:M(p)-1  % 固定i
            Numerical_right_vector = zeros(M(p)-1, 1);
            for j = 1:M(p)-1
                Numerical_right_vector(j) = numerical(i+1,j);
            end
            
            Numerical_right_vector(1) = Numerical_right_vector(1) - alpha * Numerical(i+1,1,k+1);
            Numerical_right_vector(M(p)-1) = Numerical_right_vector(M(p)-1) - alpha * Numerical(i+1,M(p)+1,k+1);
            
            Numerical(i+1,2:M(p),k+1) = chase(a, b, c, Numerical_right_vector);
        end
    end

    error = abs(Numerical(:,:,end) - Accurate(:,:,end));
    error_inf(p) = max(error(:));
    
    % 绘图
    figure(p);
    [X, Y] = meshgrid(y, x);
    
    subplot(1,3,1);
    surf(X, Y, Accurate(:,:,end));
    xlabel('x'); ylabel('y'); zlabel('Accurate');
    title('精确解');
    grid on;
    
    subplot(1,3,2);
    surf(X, Y, Numerical(:,:,end));
    xlabel('x'); ylabel('y'); zlabel('Numerical');
    title('数值解');
    grid on;
    
    subplot(1,3,3);
    surf(X, Y, error);
    xlabel('x'); ylabel('y'); zlabel('error');
    title('误差');
    grid on;
end

% 计算收敛阶
for k = 2:length(M)
    X = error_inf(k-1) / error_inf(k);
    Norm(k-1) = X;
end

figure(length(N)+1);
plot(1:length(N)-1, Norm, '-b^');
xlabel('序号'); ylabel('误差阶数');
title('紧致ADI格式误差阶');
grid on;
toc


