%Conjugate Gradient Method implementation for minimizing a function.  
%The tolerance is value of residual. Function should be inputed as matrix.
%by Huiting Su.

function Xk = CGM(A,b,Xk,tol)
    n=length(b);
    rk=A*Xk-b; %rk is residual in k^(th) iteration
    pk=-rk;

    for i=1:n
        alphak = -rk.'* pk / (pk'*A*pk);
        rk1 = rk + alphak*A*pk;      %rk1 is r_(k+1)
        Xk = Xk + alphak*pk;
        if sqrt(rk1)<tol
              break;
        end
        pk = -rk1 + (rk1.'*A*pk) /( pk.'*A*pk)*pk;
        rk=rk1;
    end
end