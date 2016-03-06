function out=fastica(X,method)
% 输入矩阵X，行为特征数，列为采样数
% 选择方法method：1.定点   2. 最大熵信息
[row, col]=size(X);

mean_X=mean(X,2);
X=X-repmat(mean_X,1,col);
[V,D]=eig(X*X');
X=V*D^(-1/2)*V'*X;

%FastIca算法
%求完后 W 是 混合矩阵 A 的近似
W=rand(row);
maxTimex=1000;
for p=1:row
    W(:,p)=W(:,p)/norm(W(:,p));
    times=1;
    while times<maxTimex
        pred=W(:,p);        
        switch method
            case {1}                
                W(:,p)=1/col*X*((W(:,p)'*X).^3)'-3*W(:,p);                    %定点法
            case {2}
%               for j=1:row
%                    temp(j,1)=mean(X(j,:).*tanh(W(:,j)'*X));                %最大熵信息
%               end
%                W(:,p)=temp-mean(1-(tanh(W(:,p)'*X)).^2)*W(:,p);
                 W(:,p)=mean(X.*tanh(W'*X),2)-mean(1-(tanh(W(:,p)'*X)).^2)*W(:,p);
        end

        %采用收缩策略，类似于正交化的过程
        sum_1=zeros(row,1);
        for i=1:p-1
            sum_1=sum_1+W(:,p)'*W(:,i)*W(:,i);
        end
        W(:,p)=W(:,p)-sum_1;
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
 
