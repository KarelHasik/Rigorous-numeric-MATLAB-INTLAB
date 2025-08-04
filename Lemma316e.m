X=intval('-0.61');
longprecision(32);
j=1;

while (inf(X)<=-0.24837) && (inf(C(X)) >sup(X))
    % Save rounding mode
rndold = getround;

setround(-1); a_inf = inf(C(X));   % round down
setround(1);  a_sup = inf(C(X));   % round up

setround(rndold);
X = infsup(a_inf, a_sup);             % now contains the true infimum

    j=j+1; 
end

disp(['j = ', num2str(j), ', inf(X) = ', num2str(inf(X),10), ', inf(C(X)) = ', num2str(inf(C(X)),10)]) 

function y=C(m)
Am= -(intval(37)./intval(24))*m/(1+m)-1+(1+m)*log(1+m)/m; 
y=LR(Am,m);
end

%j = 92, inf(X) = -0.2483696639, inf(C(X)) = -0.2483692024