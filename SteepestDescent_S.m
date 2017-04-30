%Steepest Descent Method implementation for minimizing a function.  
%The tolerance is the distance between two successive points. Input
%function is "MultiV_f.m"
%When tolerence is low, like 0.0001, and the input function is complex, the running time can be long. 
%by Huiting Su.

function SteepestDescent_S(X0,tol)  %x0 should be a row vector
    tic
    Xk=X0;
    n=length(X0);
    X = sym('x',[1 n]); 
    fX=MultiV_f(X);
    gX=gradient(fX);
    while 1
        gk=eval(subs(gX,X,Xk));    %calculate gradient at point x_k
        if abs(gk)<tol
            break
        end
        dk=-1.*gk;
        display(Xk);
        %line search
        t=0.05; c1=0.1; c2=0.9; beta=0.1;
        fk=MultiV_f(Xk);
        while 1
            Xk1=Xk + t.*transpose(dk);
            fk1=MultiV_f(Xk1);  %update function value and gradient at point x_(k+1)
            gk1=eval(subs(gX,X,Xk1));  %calculate gradient at point x_(k+1)
            
            if ( fk1 < fk + c1*t*transpose(dk)*gk ) || ( transpose(dk)*gk1 > c2*transpose(dk)*gk )  %wolfe condition
                break
            end
            t=beta*t;
        end
        
        Xk=Xk1;
    end
    display(Xk1);
    toc
end

