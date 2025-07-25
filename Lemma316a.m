X(1)=intval('-0.61');
m(1)=inf(X(1));
Z(2)=LR(-X(1),X(1));
m(2)=inf(Z(2));
X(2)=intval(m(2));
j=2;
while (m(j-1)<m(j)) && (m(j) <-0.009)
    j=j+1;
    Z(j)=LR(-X(j-1),X(j-1));
    m(j)=inf(Z(j));
    X(j)=intval(m(j)); 
end
j
m(j-1)
m(j)
