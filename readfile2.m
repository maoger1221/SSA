function [data1]=readfile2(data)
    k=0;sum=0;n=1;
    for i=2:length(data)
        if data(i,4)==data(i-1,4)
            sum=sum+data(i,9);
            k=k+1;
        else
            data1(n)=sum/k;
            n=n+1;
            k=1;
            sum=data(i,9);
        end
    end
end