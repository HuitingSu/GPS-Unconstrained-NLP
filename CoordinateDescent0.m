%Coordinate Descent method implementation for minimizing a function.  
%The tolerance is derivative on each coordinate. Input function is "MultiV_f.m"
%by Huiting Su.

function CoordinateDescent0(X0,tol)  %x0 should be a row vector
    tic
    Xk=X0;
    n=length(X0);
   
    while 1
        Xkj=Xk;
        counter=0;
        for j=1:n
            
            %patial derivative of X(j) at Xkj
 
            h= 0.00001; %sqrt(eps)*Xk(i);
            f_Sample=zeros(5,1);
            for i=1:5
                X_Sample=Xkj;
                X_Sample(j)=Xkj(j)+(i-3)*h;
                f_Sample(i)=MultiV_f(X_Sample);
            end
            gkj= (1*f_Sample(1)-8*f_Sample(2)+0*f_Sample(3)+8*f_Sample(4)-1*f_Sample(5))/(12*h);
            
            %if the derivative is already smaller than tolerance, then just
            %skip this dimension
            if abs(gkj)<tol
                continue
            end
            counter=counter+1;
            dk=-1*gkj;  %dk is 1-dimension
    
            %line search
            t=1000; c1=0.1; c2=0.9; beta=0.1;  %need to choose a relatively large t
            fk=MultiV_f(Xkj);
 
            while 1
                Xk1j=Xkj;
                Xk1j(j)=Xkj(j) + t.*transpose(dk);
                fk1=MultiV_f(Xk1j);  %update function value and gradient at point x_(k+1)
                
                %calculate patial derivative at point x_(k+1)
                for i=1:5
                    X_Sample=Xk1j;
                    X_Sample(j)=Xk1j(j)+(i-3)*h;
                    f_Sample(i)=MultiV_f(X_Sample);
                end
                gk1j= (1*f_Sample(1)-8*f_Sample(2)+0*f_Sample(3)+8*f_Sample(4)-1*f_Sample(5))/(12*h);

                if ( fk1 <= fk + c1*t*transpose(dk)*gkj ) && ( transpose(dk)*gk1j >= c2*transpose(dk)*gkj )  %wolfe condition
                    break
                end
                t=beta*t;
            end
            Xkj=Xk1j;
        end
        if counter==0
            break
        end
        Xk1=Xk1j;
%         f=0;
%         for i=1:n 
%             if  abs(Xk1(i)-Xk(i))>tol
%                 f=f+1;
%             end
%         end
%         if  f==0
%             break
%         end
       
%         gk1=g_Xk(Xk1);
%         if sqrt(gk1.'*gk1)<tol
%             break
%         end

        %check that the derivative along all dimension are below tol 
%         f=0;
%         for i=1:n 
%             if  abs(gk1(i))>tol
%                 f=f+1;
%             end
%         end
%         if  f==0
%             break
%         end


        Xk=Xk1;
        
    end
    display(Xk1);
    fmin=MultiV_f(Xk1);
    display(fmin);
    toc
end

