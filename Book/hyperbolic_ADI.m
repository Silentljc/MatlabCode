% 双曲ADI

tic;
clear
clc
M=[10,20,40,80,160];% 空间网格数量
N=M;% 时间网格数量
for p=1:length(M)
    h=1/M(p);% 这里定义空间步长等距
    tau=1/N(p); % 时间步长
    x=0:h:1;
    y=0:h:1;
    t=0:tau:1;
    Numerical(M(p)+1,M(p)+1,N(p)+1)=0;%u
    numerical(M(p)+1,M(p)-1)=0;%u*

    r = tau^2/(h^2);
    alpha = -r/2;
    beta = 1+r;
    % 求解u*ij和uij过程中构造三对角矩阵
    % a 表示下对角线元素
    % b 表示主对角线元素
    % c 表示上对角线元素
    a=alpha*ones(M(p)-2,1);
    b=beta*ones(M(p)-1,1);
    c=alpha*ones(M(p)-2,1);

    for i=1:M(p)+1
        for j=1:M(p)+1
            Numerical(i,j,1)=exp(1/2*(x(i)+y(j)));% 初值Numerical（x,y,0）=u(i,j,0)
        end
    end
    for j=1:M(p)+1
        for k=1:N(p)+1
            Numerical(1,j,k)=exp(1/2*y(j)-t(k));%  边值Numerical（0,y,t）=u(0,j,k)
        end
    end
    for j=1:M(p)+1
        for k=1:N(p)+1
            Numerical(M(p)+1,j,k)=exp(1/2*(1+y(j))-t(k));%  边值Numerical（1,y,t）=u(m,j,k)
        end
    end
    for i=1:M(p)+1
        for k=1:N(p)+1
            Numerical(i,1,k)=exp(1/2*x(i)-t(k));%  边值Numerical（x,0,t）=u(i,0,k)
        end
    end
    for i=1:M(p)+1
        for k=1:N(p)+1
            Numerical(i,M(p)+1,k)=exp(1/2*(1+x(i))-t(k));%  边值Numerical（x,1,t）=u(i,m,k)
        end
    end
    f = @(x,y,t) 1/2*exp(1/2*(x+y)-t);
    fun = @(x,y,t) exp(1/2*(x+y)-t);
    psi = @(x,y) -exp(1/2*(x+y));
    rou = @(x,y) -exp(1/2*(x+y));
    % 计算精确解
    for i=1:M(p)+1
        for j=1:M(p)+1
            for k=1:N(p)+1  % 修正：应该是N(p)+1而不是M(p)+1
                Accurate(i,j,k)=fun(x(i),y(j),t(k));
            end
        end
    end

    % 计算Numerical（x,y,1）=u(i,j,1)
    for j=1:M(p)-1% 固定j
        numerical(1,j)=alpha*Numerical(1,j,2)+...
            beta*Numerical(1,j+1,2)+...
            alpha*Numerical(1,j+2,2);%  u*0j1
        numerical(M(p)+1,j)=alpha*Numerical(M(p)+1,j,2)+...
            beta*Numerical(M(p)+1,j+1,2)+...
            alpha*Numerical(M(p)+1,j+2,2);%  u*mj1
        % 循环生成1式右端列向量
        for i=1:M(p)-1
            numerical_right_vector(i,1)=Numerical(i+1,j+1,1)+tau*psi(x(i+1),y(j+1))-tau^3*rou(x(i+1),y(j+1))/3+...
                tau^2*f(x(i+1),y(j+1),t(2))/2;
        end
        % 添加边界贡献
        numerical_right_vector(1,1)=numerical_right_vector(1,1)-alpha*numerical(1,j);
        numerical_right_vector(M(p)-1,1)=numerical_right_vector(M(p)-1,1)-alpha*numerical(M(p)+1,j);
        numerical(2:M(p),j)=chase(a,b,c,numerical_right_vector);
    end
    for i=1:M(p)-1 % 固定i
        for j=1:M(p)-1
            Numerical_right_vector(j,1)=numerical(i+1,j);
        end
        Numerical_right_vector(1,1)=Numerical_right_vector(1,1)-alpha*Numerical(i+1,1,2);
        Numerical_right_vector(M(p)-1,1)=Numerical_right_vector(M(p)-1,1)-alpha*Numerical(i+1,M(p)+1,2);
        Numerical(i+1,2:M(p),2)=chase(a,b,c,Numerical_right_vector);
    end

    % 核心部分
    for k=2:N(p)
        for j=1:M(p)-1% 固定j
            numerical(1,j)=alpha*Numerical(1,j,k+1)+...
                beta*Numerical(1,j+1,k+1)+...
                alpha*Numerical(1,j+2,k+1);%  u*0j
            numerical(M(p)+1,j)=alpha*Numerical(M(p)+1,j,k+1)+...
                beta*Numerical(M(p)+1,j+1,k+1)+...
                alpha*Numerical(M(p)+1,j+2,k+1);%  u*mj
            % 循环生成1式右端列向量
            for i=1:M(p)-1
                numerical_right_vector(i,1)=2*Numerical(i+1,j+1,k)+f(x(i+1),y(j+1),t(k))*tau^2+...
                    (r/2)*(Numerical(i,j+1,k-1)+Numerical(i+2,j+1,k-1))+...
                    (r/2)*(Numerical(i+1,j,k-1)+Numerical(i+1,j+2,k-1))-(1+r+r)*Numerical(i+1,j+1,k-1);
            end
            % 添加边界贡献
            numerical_right_vector(1,1)=numerical_right_vector(1,1)-alpha*numerical(1,j);
            numerical_right_vector(M(p)-1,1)=numerical_right_vector(M(p)-1,1)-alpha*numerical(M(p)+1,j);
            numerical(2:M(p),j)=chase(a,b,c,numerical_right_vector);
        end
        for i=1:M(p)-1 % 固定i
            for j=1:M(p)-1
                Numerical_right_vector(j,1)=numerical(i+1,j);
            end
            Numerical_right_vector(1,1)=Numerical_right_vector(1,1)-alpha*Numerical(i+1,1,k+1);
            Numerical_right_vector(M(p)-1,1)=Numerical_right_vector(M(p)-1,1)-alpha*Numerical(i+1,M(p)+1,k+1);
            Numerical(i+1,2:M(p),k+1)=chase(a,b,c,Numerical_right_vector);
        end
    end

    error=abs(Numerical(:,:,:)-Accurate(:,:,:));  
    error_inf(p)=max(error(:));
    figure(p)
    [X,Y]=meshgrid(y,x);
    subplot(1,3,1)
    surf(X,Y,Accurate(:,:,end));  % 修正：使用end而不是M(p)
    xlabel('x');ylabel('y');zlabel('Accurate');  % 修正：标签改为Accurate
    title('精确解图像');
    grid on;
    subplot(1,3,2)
    surf(X,Y,Numerical(:,:,end));  % 修正：使用end而不是M(p)
    xlabel('x');ylabel('y');zlabel('Numerical');
    title('数值解图像');
    grid on;
    subplot(1,3,3)
    surf(X,Y,error(:,:,end));
    xlabel('x');ylabel('y');zlabel('error');
    title('误差图像');
    grid on;
end
for k=2:length(M)
    X=error_inf(k-1)/error_inf(k);
    Norm(k-1)=X;  
end
figure(length(N)+1)
plot(1:length(N)-1,Norm,'-b^');
xlabel('序号');ylabel('误差阶数');
title('ADI格式误差阶');
grid on;
toc;