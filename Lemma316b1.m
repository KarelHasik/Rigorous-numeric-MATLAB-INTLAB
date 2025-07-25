tic %here we consider the sub−interval $[−0.25, −0.13]$
for j=1:300
m1(j)=infsup(-0.13-(j-1)*0.0004,-0.13-(j-1)*0.0004);
m2(j)=infsup(-0.13-j*0.0004,-0.13-j*0.0004);
a2(j)=-m1(j);
A(j)=-(37/24)*m2(j)/(1+m2(j))-1+(1+m2(j))*log(1+m2(j))/m2(j);
N(j)= floor(sup((A(j)-a2(j)))*10^3)+1;
h=intval('0.001');
mY(j)= hull(m2(j),m1(j));
for k=1:N(j)
W(j,k)= sup(dL(hull(a2(j)+(k-1)*h,a2(j)+k*h),mY(j)));
end
end
max(max(W))
toc
%ans =0.815039474318645, Elapsed time is 34.254236 seconds.