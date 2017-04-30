%quasi-Newton method SR1 implementation for solving equations
%The tolerance is abs of gradient. Input function is "MultiV_f.m".
%by Huiting Su.

function SR1_S(X0,tol)   %x0 should be a row vector
    tic    
    Xk=X0;
    n=length(X0);
    X = sym('x',[1 n]); 
    fX=MultiV_f(X);
    gX=gradient(fX);
    Hk=eye(n);
    while 1
        gk=eval(subs(gX,X,Xk));    %calculate gradient at point x_k
        if abs(gk)<tol
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
            gk1=eval(subs(gX,X,Xk1));  %calculate gradient at point x_(k+1)
            
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
    toc
end

