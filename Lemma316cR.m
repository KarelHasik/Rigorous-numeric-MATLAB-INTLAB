X = intval('-0.009');
j = 1;
rndold = getround;  % Save rounding mode
while (sup(X) >= -1/9) && (sup(C(X)) < inf(X))
    setround(1);                         % upward rounding
    X = intval(sup(C(X)));
    setround(rndold);                    % restore rounding mode
    j = j + 1; 
end

disp('Final result:')
disp(['sup(X) = ', num2str(sup(X))])
disp(['j = ', num2str(j)])
disp(['sup(C(X)) = ', num2str(sup(C(X)))])


function y = C(m)
  a = -intval(37)/intval(24);
  one = intval(1);
  Am = a*m/(one + m) - one + (one + m)*log(one + m)/m;
  y = LR(Am, m);  
end
