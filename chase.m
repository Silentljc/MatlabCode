function x = chase(a, b, c, d)
% 追赶法求解三对角线性方程组 Ax = d
% 输入:
%   a - 下对角线元素向量 (长度 n-1)
%   b - 主对角线元素向量 (长度 n)  
%   c - 上对角线元素向量 (长度 n-1)
%   d - 右端向量 (长度 n)
% 输出:
%   x - 解向量

    n = length(b);
    % 初始化
    alpha = zeros(n, 1);
    beta = zeros(n, 1);
    x = zeros(n, 1);
    
    % 追过程
    alpha(1) = b(1);
    beta(1) = d(1) / alpha(1);
    for i = 2:n
        alpha(i) = b(i) - a(i-1) * c(i-1) / alpha(i-1);
        beta(i) = (d(i) - a(i-1) * beta(i-1)) / alpha(i);
    end
    
    % 赶过程
    x(n) = beta(n);
    for i = n-1:-1:1
        x(i) = beta(i) - c(i) * x(i+1) / alpha(i);
    end
end