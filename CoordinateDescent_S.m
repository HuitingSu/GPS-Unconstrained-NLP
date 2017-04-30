%Coordinate Descent method implementation for minimizing a function.  
%The tolerance is the distance between two successive points. Input
%function is "MultiV_f.m"
%by Huiting Su.

function CoordinateDescent_S(X0,tol)  %x0 should be a row vector
    tic
    Xk=X0;
    n=length(X0);
    X = sym('x',[1 n]); 
    fX=MultiV_f(X);
    
    while 1
        Xkj=Xk;
        for j=1:n
            
            gXj=diff(fX,X(j));
            gk=eval(subs(gXj,X,Xkj));
            dk=-1*gk;  %dk is 1-dimension
    
            %line search
            t=0.1; c1=0.1; c2=0.9; beta=0.1;
            fk=MultiV_f(Xkj);
 
            while 1
                Xk1j=Xkj;
                Xk1j(j)=Xkj(j) + t.*transpose(dk);
                fk1=MultiV_f(Xk1j);  %update function value and gradient at point x_(k+1)
                gk1=eval(subs(gXj,X,Xk1j));  %calculate gradient at point x_(k+1)
             
                if ( fk1 <= fk + c1*t*dk*gk ) && ( dk*gk1 >= c2*dk*gk )  %wolfe condition
                    break
                end
                %break
                t=beta*t;
            end
            Xkj=Xk1j;
        end
        Xk1=Xk1j;
        f=0;
        for i=1:n 
            if  abs(Xk1(i)-Xk(i))>tol
                f=f+1;
            end
        end
        if  f==0
            break
        end
        Xk=Xk1;
        
    end
    display(Xk1);
    toc
end

