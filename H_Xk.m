function [Hk]=H_Xk(Xk)
    n=length(Xk);
    Hk=zeros(n,n);
    f_Sample=zeros(5,1);
    h= 0.00001; %sqrt(eps);
    
    for i=1:n
        X_Sample=Xk;
        X_Sample(i)=Xk(i)-h;
        f_Sample(1)=MultiV_f(X_Sample);
        f_Sample(2)=MultiV_f(Xk);
        X_Sample=Xk;
        X_Sample(i)=Xk(i)+h;
        f_Sample(3)=MultiV_f(X_Sample);
        
        Hk(i,i)=(f_Sample(1)-2*f_Sample(2)+f_Sample(3))/h^2;
    end
    
    for i=1:n      
        for j=1:n
            if i~=j
                X_Sample=Xk;
                X_Sample(i)=Xk(i)-h;
                X_Sample(j)=Xk(j)-h;
                f_Sample(1)=MultiV_f(X_Sample);
                
                X_Sample=Xk;
                X_Sample(i)=Xk(i)-h;
                X_Sample(j)=Xk(j)+h;
                f_Sample(2)=MultiV_f(X_Sample);
                
                X_Sample=Xk;
                X_Sample(i)=Xk(i)+h;
                X_Sample(j)=Xk(j)-h;
                f_Sample(3)=MultiV_f(X_Sample);
                
                X_Sample=Xk;
                X_Sample(i)=Xk(i)+h;
                X_Sample(j)=Xk(j)+h;
                f_Sample(4)=MultiV_f(X_Sample);

                Hk(i,j)=(f_Sample(1) - f_Sample(2) - f_Sample(3) + f_Sample(4))/(4*h^2);
            end
        end
    end
    
end