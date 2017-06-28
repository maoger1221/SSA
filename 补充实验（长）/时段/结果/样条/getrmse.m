function rmse=getrmse(data1,data2)
    sum=0;
    for i=1:length(data1)
        sum=sum+(data1(i)-data2(i))*(data1(i)-data2(i));
    end
    rmse=sqrt(sum/length(data1));

end