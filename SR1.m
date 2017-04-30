%quasi-Newton method SR1 implementation for solving equations
%The tolerance is abs of gradient. Input function is "MultiV_f.m".
%by Huiting Su.

function SR1(X0,tol)   %x0 should be a row vector
    tic    
    Xk=X0;
    n=length(X0);
    Hk=eye(n);
    while 1
        gk=g_Xk(Xk);     %calculate gradient at point x_k
%         if sqrt(gk.'*gk)<tol
%             break
%         end
    
    %check that the derivative along all dimension are below tol 
        f=0;
        for i=1:n 
            if  abs(gk(i))>tol
                f=f+1;
            end
        end
        if  f==0
            break
        end
        
        dk=-1.*Hk*gk;
        display(Xk);
        %line search
        t=0.1; c1=0.1; c2=0.9; beta=0.1;
        fk=MultiV_f(Xk);
        while 1
            Xk1=Xk + t.*transpose(dk);
            fk1=MultiV_f(Xk1);  %update function value and gradient at point x_(k+1)
            gk1=g_Xk(Xk1);  %calculate gradient at point x_(k+1)
            
            if ( fk1 < fk + c1*t*transpose(dk)*gk ) || ( transpose(dk)*gk1 > c2*transpose(dk)*gk )  %wolfe condition
                break
            end
            t=beta*t;
        end
        sk=(Xk1-Xk).';  %column
        yk=gk1-gk;  %column
        Hk=Hk+(sk-Hk*yk)*(sk-Hk*yk).'/((sk-Hk*yk).'*yk);
        
        Xk=Xk1;
    end
    display(Xk1);
    fmin=MultiV_f(Xk1);
    display(fmin);
    toc
end

