clear
clc
 
%�۲��źž��� X Ϊ K*N ����
%Դ�źž��� S Ϊ K*N ����
%��Ͼ��� A Ϊ K*K ����
 
K=2;
%���� K ��Դ�ź�, �洢�� S
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
title('�ź�1');
xlabel('time/ms','fontname','Times New Roman','fontsize',9);
ylabel('amplitude','fontname','Times New Roman','fontsize',9);
hold on;
 
S(2,:)=sin(2*pi*f1*t)+sin(2*pi*f2*t); 
subplot(2,1,2);
plot(t,S(2,:),'k');
set(gca,'fontname','Times New Roman','fontsize',9);
title('�ź�2');xlabel('time/ms','fontname','Times New Roman','fontsize',9);
ylabel('amplitude','fontname','Times New Roman','fontsize',9);
hold on;
 
%���� K ������ź�,�洢�� X, A �д洢��Ͼ���
A=rand(2);
%A=rand(1,2);
X=A*S;
figure(2);
subplot(2,1,1);
plot(t,X(1,:),'k');
set(gca,'fontname','Times New Roman','fontsize',9);
title('�۲��ź�1','fontname','����','fontsize',9);xlabel('time/ms','fontname','Times New Roman','fontsize',9);
ylabel('amplitude','fontname','Times New Roman','fontsize',9);
subplot(2,1,2);
plot(t,X(2,:),'k');
set(gca,'fontname','Times New Roman','fontsize',9);
title('�۲��ź�2','fontname','����','fontsize',9);xlabel('time/ms','fontname','Times New Roman','fontsize',9);
ylabel('amplitude','fontname','Times New Roman','fontsize',9);

%D=corrcoef(X1,c7);
% X �������Ļ�
m=mean(X,2);
for i=1:N
    X(:,i)=X(:,i)-m;
end
 
% X ���ݰ׻�
covMat=cov(X');
[E,D]=eig(covMat);
V=E*D^(-0.5)*E';
X=V*X;
figure(7);xlabel('time/ms');
ylabel('amplitude/um');
plot(t,X);
title('�׻�����ź�');

%%
[row, col]=size(X);
a=rand(row);
%FastIca�㷨
%����� W �� ��Ͼ��� A �Ľ���
W=a;
maxTimex=1000;
for p=1:row
    W(:,p)=W(:,p)/norm(W(:,p));
    times=1;
    while times<maxTimex
        pred=W(:,p);        
        method=2;                                                           %ѡ�񷽷���1.����   2. �������Ϣ
        switch method
            case {1}                
                W(:,p)=1/N*X*((W(:,p)'*X).^3)'-3*W(:,p);                    %���㷨
            case {2}
               for j=1:row
                    temp(j,1)=mean(X(j,:).*tanh(W(:,p)'*X));                %�������Ϣ
               end
                W(:,p)=temp-mean(1-(tanh(W(:,p)'*X)).^2)*W(:,p);
        end

        %�����������ԣ��������������Ĺ���
        sum=zeros(K,1);
        for i=1:p-1
            sum=sum+W(:,p)'*W(:,i)*W(:,i);
        end
        W(:,p)=W(:,p)-sum;
        W(:,p)=W(:,p)/norm(W(:,p));
        
        %��������
        if abs(dot(W(:,p),W(:,p)))>1-1e-22
            break
        end
        
        times=times+1;
        if times==maxTimex
               'Not converged!'
        end        
    end
end
 
%����Դ�ź�
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
