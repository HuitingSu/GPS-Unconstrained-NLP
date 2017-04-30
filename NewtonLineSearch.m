%Newton's Line Search Algorithm implementation for minimizing a function.  
%The tolerance is abs of gradient. Input function is "MultiV_f.m".
%by Huiting Su.

function NewtonLineSearch(X0,tol)  %x0 should be a row vector
    tic
    Xk=X0;
    n=length(X0);
    t0=1000000;
    while 1
        gk=g_Xk(Xk);   %calculate gradient at point x_k
%         if sqrt(gk.'*gk)<tol
%             break
%         end
 
%check that the derivative along all dimension are below tol 
        
        if abs(MultiV_f(Xk))<tol
            break
        end
        
        Hk=H_Xk(Xk);
        dk=-1.*inv(Hk)*gk;
        display(Xk);
        %line search
        
        t=t0; c1=0.1; c2=0.9; beta=0.1;k=0;
        fk=MultiV_f(Xk);
        while 1
            Xk1=Xk + t.*dk;
            fk1=MultiV_f(Xk1);  %update function value and gradient at point x_(k+1)
            gk1=g_Xk(Xk1);  %calculate gradient at point x_(k+1)
            
            if ( fk1 < fk + c1*t*transpose(dk)*gk ) && ( transpose(dk)*gk1 > c2*transpose(dk)*gk )  %wolfe condition
                break
            end
            k=k+1;
            if k>=15
                t=t0*100;
            end
            t=beta*t;
        end
        Xk=Xk1;
    end
    display(Xk1);
    fmin=MultiV_f(Xk1);
    display(fmin);
    toc
end

