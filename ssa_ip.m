function y=ssa_ip(x,L,p)%y为重构信号，在插补缺失数据时用到的ssa
%x为时间序列，L为窗口长度
%第一步，计算轨迹矩阵
    N=length(x);
    K=N-L+1;
    X=zeros(L,K);
    for i=1:L
        for j=1:K
            X(i,j)=x(i+j-1);
        end
    end
%第二步，奇异值分解
    [U,S,V]=svd(X);%矩阵X经svd分解后所得，sigma的平方是XX'和X'X的特征值，其对应的特征向量分别是U和V的列向量
    sigma=diag(S);
    sigma=-sigma;
    [sigma2,t]=sort(sigma);%排序
    sigma2=-sigma2;
    %d=length(sigma2);

    Z=zeros(L,K);
    for i=1:p
        Z=Z+sigma2(i)*U(:,t(i))*V(:,t(i))';
    end
%第四步，重构
    y=zeros(N,1);
	L2=min(L,K);
	K2=max(L,K);
    
    for k=1:N
        if k>=1 && k<=L2
           for q=1:k
              y(k)=y(k)+Z(q,k-q+1); 
           end
           y(k)=y(k)/k;
        end
        if k>=L2 && k<=K2
           for q=1:L2
              y(k)=y(k)+Z(q,k-q+1); 
           end
           y(k)=y(k)/L2;
        end
        if k>=K2 && k<=N
           for q=(k-K2+1):(N-K2+1)
              y(k)=y(k)+Z(q,k-q+1);
           end
           y(k)=y(k)/(N-k+1);
        end
    end
end