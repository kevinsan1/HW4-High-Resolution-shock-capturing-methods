function df = limiter(nc,beta,f)

df(1)=0;
df(nc)=0;
for i=2:nc-1,
    df1=f(i+1)-f(i);
    df2=f(i)-f(i-1);
    df(i)=limit(df1,df2,beta);
end

