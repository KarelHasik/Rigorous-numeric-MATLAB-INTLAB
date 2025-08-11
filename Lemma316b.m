X(1)=intval('-0.25');
m(1)=inf(X(1));
b=intval('0.009401');
Z(2)=LR(b,X(1));
m(2)=inf(Z(2));
X(2)=intval(m(2));
j=2;
while (m(j-1)<m(j)) && (m(j) <=-0.009297)
    j=j+1;
    Z(j)=LR(b,X(j-1));
    m(j)=inf(Z(j));
    X(j)=intval(m(j)); 
end
j
m(j-1)
m(j)
%j =5; ans =-0.009298870400434, ans =-0.009296711442578