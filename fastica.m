function out=fastica(X,method)
% �������X����Ϊ����������Ϊ������
% ѡ�񷽷�method��1.����   2. �������Ϣ
[row, col]=size(X);

mean_X=mean(X,2);
X=X-repmat(mean_X,1,col);
[V,D]=eig(X*X');
X=V*D^(-1/2)*V'*X;

%FastIca�㷨
%����� W �� ��Ͼ��� A �Ľ���
W=rand(row);
maxTimex=1000;
for p=1:row
    W(:,p)=W(:,p)/norm(W(:,p));
    times=1;
    while times<maxTimex
        pred=W(:,p);        
        switch method
            case {1}                
                W(:,p)=1/col*X*((W(:,p)'*X).^3)'-3*W(:,p);                    %���㷨
            case {2}
               for j=1:row
                    temp(j,1)=mean(X(j,:).*tanh(W(:,p)'*X));                %�������Ϣ
               end
                W(:,p)=temp-mean(1-(tanh(W(:,p)'*X)).^2)*W(:,p);
        end

        %�����������ԣ��������������Ĺ���
        sum=zeros(row,1);
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
 