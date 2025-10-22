tic;
clear
clc
M=[5,10,20,40,80];% 空间网格数量
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
    alpha = 1/12-r/2;
    beta = 10/12+r;
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
%     f=inline('1/2*exp(1/2*(x+y)-t)','x','y','t');% f(x,y,t)
%     fun=inline('exp(1/2*(x+y)-t)','x','y','t');% 精确解
%     psi=inline('-exp(1/2*(x+y))','x','y');% psi(x,y)
%     rou=inline('-exp(1/2*(x+y))','x','y');% f(x,y,t)
    f = @(x,y,t) 1/2*exp(1/2*(x+y)-t);
    fun = @(x,y,t) exp(1/2*(x+y)-t);
    psi = @(x,y) -exp(1/2*(x+y));
    rou = @(x,y) -exp(1/2*(x+y));
    % 计算精确解
    for i=1:M(p)+1
        for j=1:M(p)+1
            for k=1:M(p)+1
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
                numerical_right_vector(i,1)=(1/144)*(Numerical(i,j,1)+tau*psi(x(i),y(j))-tau^3*rou(x(i),y(j))/3+...
                    tau^2*f(x(i),y(j),t(2))/2+10*(Numerical(i+1,j,1)+tau*psi(x(i+1),y(j))-tau^3*rou(x(i+1),y(j))/3+...
                    tau^2*f(x(i+1),y(j),t(2))/2)+Numerical(i+2,j,1)+tau*psi(x(i+2),y(j))-tau^3*rou(x(i+2),y(j))/3+...
                    tau^2*f(x(i+2),y(j),t(2))/2)+...
                    (1/144)*10*(Numerical(i,j+1,1)+tau*psi(x(i),y(j+1))-tau^3*rou(x(i),y(j+1))/3+...
                    tau^2*f(x(i),y(j+1),t(2))/2+10*(Numerical(i+1,j+1,1)+tau*psi(x(i+1),y(j+1))-tau^3*rou(x(i+1),y(j+1))/3+...
                    tau^2*f(x(i+1),y(j+1),t(2))/2)+Numerical(i+2,j+1,1)+tau*psi(x(i+2),y(j+1))-tau^3*rou(x(i+2),y(j+1))/3+...
                    tau^2*f(x(i+2),y(j+1),t(2))/2)+...
                    (1/144)*(Numerical(i,j+2,1)+tau*psi(x(i),y(j+2))-tau^3*rou(x(i),y(j+2))/3+...
                    tau^2*f(x(i),y(j+2),t(2))/2+10*(Numerical(i+1,j+2,1)+tau*psi(x(i+1),y(j+2))-tau^3*rou(x(i+1),y(j+2))/3+...
                    tau^2*f(x(i+1),y(j+2),t(2))/2)+Numerical(i+2,j+2,1)+tau*psi(x(i+2),y(j+2))-tau^3*rou(x(i+2),y(j+2))/3+...
                    tau^2*f(x(i+2),y(j+2),t(2))/2);
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
                numerical_right_vector(i,1)=(r/24)*(Numerical(i,j,k-1)-2*Numerical(i+1,j,k-1)+Numerical(i+2,j,k-1)+...
                    10*(Numerical(i,j+1,k-1)-2*Numerical(i+1,j+1,k-1)+Numerical(i+2,j+1,k-1))+...
                    Numerical(i,j+2,k-1)-2*Numerical(i+1,j+2,k-1)+Numerical(i+2,j+2,k-1))+...
                    (r/24)*(Numerical(i,j,k-1)-2*Numerical(i,j+1,k-1)+Numerical(i,j+2,k-1)+...
                    10*(Numerical(i+1,j,k-1)-2*Numerical(i+1,j+1,k-1)+Numerical(i+1,j+2,k-1))+...
                    Numerical(i+2,j,k-1)-2*Numerical(i+2,j+1,k-1)+Numerical(i+2,j+2,k-1))-...
                    (1/144)*(Numerical(i,j,k-1)+10*Numerical(i+1,j,k-1)+Numerical(i+2,j,k-1)+...
                    10*(Numerical(i,j+1,k-1)+10*Numerical(i+1,j+1,k-1)+Numerical(i+2,j+1,k-1))+...
                    Numerical(i,j+2,k-1)+10*Numerical(i+1,j+2,k-1)+Numerical(i+2,j+2,k-1))+...
                    (2/144)*(Numerical(i,j,k)+10*Numerical(i+1,j,k)+Numerical(i+2,j,k)+...
                    10*(Numerical(i,j+1,k)+10*Numerical(i+1,j+1,k)+Numerical(i+2,j+1,k))+...
                    Numerical(i,j+2,k)+10*Numerical(i+1,j+2,k)+Numerical(i+2,j+2,k))+...
                    (tau^2/144)*(f(x(i),y(j),t(k))+10*f(x(i+1),y(j),t(k))+f(x(i+2),y(j),t(k))+...
                    10*(f(x(i),y(j+1),t(k))+10*f(x(i+1),y(j+1),t(k))+f(x(i+2),y(j+1),t(k)))+...
                    f(x(i),y(j+2),t(k))+10*f(x(i+1),y(j+2),t(k))+f(x(i+2),y(j+2),t(k)));
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

    error=abs(Numerical(:,:,M(p)+1)-Accurate(:,:,M(p)+1));
    error_inf(p)=max(max(error));
    figure(p)
    [X,Y]=meshgrid(y,x);
    subplot(1,3,1)
    surf(X,Y,Accurate(:,:,M(p)));
    xlabel('x');ylabel('y');zlabel('Numerical');
    title('the image of Accurate rusult');
    grid on;
    subplot(1,3,2)
    surf(X,Y,Numerical(:,:,M(p)));
    xlabel('x');ylabel('y');zlabel('Numerical');
    title('the image of Numerical');
    grid on;
    subplot(1,3,3)
    surf(X,Y,error);
    xlabel('x');ylabel('y');zlabel('error');
    title('the image of error Numerical');
    grid on;
end
for k=2:length(M)
    X=error_inf(k-1)/error_inf(k);
    Norm(k-1)=X;  
end
figure(length(N)+1)
plot(1:length(N)-1,Norm,'-b^');
xlabel('序号');ylabel('误差阶数');
title('紧ADI格式误差阶');
grid on;
toc;

