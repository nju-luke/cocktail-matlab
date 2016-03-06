clear
clc
 
%观察信号矩阵 X 为 K*N 矩阵
%源信号矩阵 S 为 K*N 矩阵
%混合矩阵 A 为 K*K 矩阵
 
K=2;
%生成 K 个源信号, 存储到 S
f1=30;f2=15;f3=120;
Fs=3000;
N=1000;
t=0:1/Fs:(N-1)/Fs;
 
figure(1);
%S(1,:)=sin(2*pi*f3*t).*[1+1.5*sin(2*pi*f2*t)]; 
S(1,:)=0.3*cos(2*pi*f1*t+0.5*sin(2*pi*f2*t))+0.2*sin(2*pi*f3*t);
subplot(2,1,1);
plot(t,S(1,:),'k');
set(gca,'fontname','Times New Roman','fontsize',9);
title('信号1');
xlabel('time/ms','fontname','Times New Roman','fontsize',9);
ylabel('amplitude','fontname','Times New Roman','fontsize',9);
hold on;
 
S(2,:)=sin(2*pi*f1*t)+sin(2*pi*f2*t); 
subplot(2,1,2);
plot(t,S(2,:),'k');
set(gca,'fontname','Times New Roman','fontsize',9);
title('信号2');xlabel('time/ms','fontname','Times New Roman','fontsize',9);
ylabel('amplitude','fontname','Times New Roman','fontsize',9);
hold on;
 
%生成 K 个混合信号,存储到 X, A 中存储混合矩阵
A=rand(2);
%A=rand(1,2);
X=A*S;
figure(2);
subplot(2,1,1);
plot(t,X(1,:),'k');
set(gca,'fontname','Times New Roman','fontsize',9);
title('观测信号1','fontname','宋体','fontsize',9);xlabel('time/ms','fontname','Times New Roman','fontsize',9);
ylabel('amplitude','fontname','Times New Roman','fontsize',9);
subplot(2,1,2);
plot(t,X(2,:),'k');
set(gca,'fontname','Times New Roman','fontsize',9);
title('观测信号2','fontname','宋体','fontsize',9);xlabel('time/ms','fontname','Times New Roman','fontsize',9);
ylabel('amplitude','fontname','Times New Roman','fontsize',9);

%D=corrcoef(X1,c7);
% X 数据中心化
m=mean(X,2);
for i=1:N
    X(:,i)=X(:,i)-m;
end
 
% X 数据白化
covMat=cov(X');
[E,D]=eig(covMat);
V=E*D^(-0.5)*E';
X=V*X;
figure(7);xlabel('time/ms');
ylabel('amplitude/um');
plot(t,X);
title('白化后的信号');

%%
[row, col]=size(X);
a=rand(row);
%FastIca算法
%求完后 W 是 混合矩阵 A 的近似
W=a;
maxTimex=1000;
for p=1:row
    W(:,p)=W(:,p)/norm(W(:,p));
    times=1;
    while times<maxTimex
        pred=W(:,p);        
        method=2;                                                           %选择方法：1.定点   2. 最大熵信息
        switch method
            case {1}                
                W(:,p)=1/N*X*((W(:,p)'*X).^3)'-3*W(:,p);                    %定点法
            case {2}
               for j=1:row
                    temp(j,1)=mean(X(j,:).*tanh(W(:,p)'*X));                %最大熵信息
               end
                W(:,p)=temp-mean(1-(tanh(W(:,p)'*X)).^2)*W(:,p);
        end

        %采用收缩策略，类似于正交化的过程
        sum=zeros(K,1);
        for i=1:p-1
            sum=sum+W(:,p)'*W(:,i)*W(:,i);
        end
        W(:,p)=W(:,p)-sum;
        W(:,p)=W(:,p)/norm(W(:,p));
        
        %收敛条件
        if abs(dot(W(:,p),W(:,p)))>1-1e-22
            break
        end
        
        times=times+1;
        if times==maxTimex
               'Not converged!'
        end        
    end
end
 
%推算源信号
out=W'*X;
 
figure;t=1:1000;
for k=1:2  
    subplot(2,1,k);
    plot(t,out(k,:),'k');
    set(gca,'fontname','Times New Roman','fontsize',9);
    xlabel('time/ms','fontname','Times New Roman','fontsize',9);
ylabel('amplitude/um','fontname','Times New Roman','fontsize',9);
    hold on;
  
end;
