% five_point_solver.m
% 用五点差分法求解 - (u_xx + u_yy) = (pi^2-1) e^x sin(pi y)
% 区间 x in (0,2), y in (0,1) ，Dirichlet 边界：
% u(0,y)=sin(pi y), u(2,y)=e^2 sin(pi y), u(x,0)=0, u(x,1)=0
% 精确解 u = e^x sin(pi y)

clear; clc; close all

% 网格细化参数（可改）
M = 40;   % x 方向分 M 等份 => 内部点在 i=1..M-1 (总网格 M+1 个节点)
N = 20;   % y 方向分 N 等份 => 内部点在 j=1..N-1

Lx = 2; Ly = 1;
h = Lx / M;
k = Ly / N;

% 网格节点坐标（包含边界）
x = linspace(0,Lx,M+1);
y = linspace(0,Ly,N+1);

% 内部点数
nx = M-1;
ny = N-1;
numUnknown = nx * ny;

% 右端 f(x,y) 在内部点
[Xint, Yint] = meshgrid(x(2:end-1), y(2:end-1)); % 注意 meshgrid 的顺序: 列为 x，行为 y
% 我们将按行主序（y方向外层）把点展平
f_interior = (pi^2 - 1) * exp(Xint) .* sin(pi * Yint);

% 稀疏矩阵系数
% 对角 a = 2/h^2 + 2/k^2，邻点 b = -1/h^2 (左右), c = -1/k^2 (上下)
Ax = 1/h^2;
Ay = 1/k^2;
diag_center = 2*Ax + 2*Ay;

% 预分配稀疏矩阵三对角位置（五点模板）
% 我们使用 i index for x (1..nx), j index for y (1..ny)
% 全局编号 p = (j-1)*nx + i
i_idx = []; j_idx = []; val = [];

% helper to push entry
push = @(ii,jj,v) deal(i_idx(end+1+0), j_idx(end+1+0), val(end+1+0)); % dummy to avoid warnings
% we'll not use push; directly append

i_idx = []; j_idx = []; val = [];

for j = 1:ny
    for i = 1:nx
        p = (j-1)*nx + i;
        % center
        i_idx(end+1) = p; j_idx(end+1) = p; val(end+1) = diag_center;
        % left neighbor (i-1,j)
        if i>1
            q = p - 1;
            i_idx(end+1) = p; j_idx(end+1) = q; val(end+1) = -Ax;
        end
        % right neighbor (i+1,j)
        if i<nx
            q = p + 1;
            i_idx(end+1) = p; j_idx(end+1) = q; val(end+1) = -Ax;
        end
        % bottom neighbor (i,j-1)
        if j>1
            q = p - nx;
            i_idx(end+1) = p; j_idx(end+1) = q; val(end+1) = -Ay;
        end
        % top neighbor (i,j+1)
        if j<ny
            q = p + nx;
            i_idx(end+1) = p; j_idx(end+1) = q; val(end+1) = -Ay;
        end
    end
end

A = sparse(i_idx, j_idx, val, numUnknown, numUnknown);

% 组装 RHS（含边界贡献）
F = reshape(f_interior', numUnknown, 1); % 注意 reshape 的方向，保持与编号一致
% 将边界值移到 RHS：左、右、下、上边界
% 边界函数
g_left  = @(y) sin(pi*y);         % x=0
g_right = @(y) exp(2) .* sin(pi*y); % x=2
g_bottom = @(x) zeros(size(x));   % y=0
g_top    = @(x) zeros(size(x));   % y=1

% 对每个内部点添加边界贡献
for j = 1:ny
    for i = 1:nx
        p = (j-1)*nx + i;
        xi = x(i+1);    % 实际 x 坐标（i 从 1 对应 x(2)）
        yj = y(j+1);
        % 若左邻为边界 (i==1) ： u_{i-1,j} = g_left(yj)
        if i == 1
            F(p) = F(p) + Ax * g_left(yj);
        end
        % 若右邻为边界 (i==nx) ： u_{i+1,j} = g_right(yj)
        if i == nx
            F(p) = F(p) + Ax * g_right(yj);
        end
        % 若下邻为边界 (j==1) ： u_{i,j-1} = g_bottom(xi)
        if j == 1
            F(p) = F(p) + Ay * g_bottom(xi);
        end
        % 若上邻为边界 (j==ny) ： u_{i,j+1} = g_top(xi)
        if j == ny
            F(p) = F(p) + Ay * g_top(xi);
        end
    end
end

% 求解线性系统
Uvec = A \ F;

% 把解展开为网格矩阵（含边界）
U = zeros(ny+2, nx+2); % 行对应 y (0..N), 列对应 x (0..M)
% 填入内部解（注意 reshape 与转置）
Uint = reshape(Uvec, nx, ny)'; % first index row (y), second index col (x)
U(2:end-1, 2:end-1) = Uint;

% 填入边界（Dirichlet）
% x=0 列
U(:,1) = sin(pi * y');           % y is row vector
% x=2 列
U(:,end) = exp(2) * sin(pi * y') ;
% y=0 行
U(1,:) = 0;
% y=1 行
U(end,:) = 0;

% 精确解在网格上的取值
[Xg, Yg] = meshgrid(x, y);
U_exact = exp(Xg) .* sin(pi * Yg);

% 误差计算
err = abs(U - U_exact);
max_err = max(err(:));
% 离散 L2 误差（近似）
L2_err = sqrt( sum(err(:).^2) * h * k );

fprintf('网格: M=%d, N=%d, 内部未知数=%d\n', M, N, numUnknown);
fprintf('||error||_inf = %g\n', max_err);
fprintf('approx L2 error = %g\n', L2_err);

% 绘图：数值解、解析解、误差
figure;
surf(Xg, Yg, U); xlabel('x'); ylabel('y'); zlabel('U_num'); title('数值解 U');
shading interp; view(45,30);

figure;
surf(Xg, Yg, U_exact); xlabel('x'); ylabel('y'); zlabel('U_exact'); title('精确解 U');
shading interp; view(45,30);

figure;
surf(Xg, Yg, err); xlabel('x'); ylabel('y'); zlabel('abs error'); title('绝对误差 |U-U_{exact}|');
shading interp; view(45,30);

% 可选：打印部分网格误差
disp('最大误差的位置 (近似)：');
[indy,indx] = find(err == max_err);
disp([x(indx), y(indy), max_err]);
