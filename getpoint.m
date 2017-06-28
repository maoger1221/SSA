function y=getpoint(data,start,num)%选插值点,先选取start点和start点+1，再向左向右各选七个，共16个点，num为缺失点个数
    y=[data(start);data(start+1)];
    for i=1:7
        y=[data(start-i*(num+1));y;data(start+1+i*(num+1))];
    end
end