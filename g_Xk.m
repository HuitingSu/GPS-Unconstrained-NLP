function [gk]=g_Xk(Xk)
    n=length(Xk);
    gk=zeros(n,1);
    f_Sample=zeros(5,1);
    for i=1:n
         h= 0.00001; %sqrt(eps)*Xk(i);
        for j=1:5
            X_Sample=Xk;
            X_Sample(i)=Xk(i)+(j-3)*h;
            f_Sample(j)=MultiV_f(X_Sample);
        end
            gk(i)= (1*f_Sample(1)-8*f_Sample(2)+0*f_Sample(3)+8*f_Sample(4)-1*f_Sample(5))/(12*h);
    end
    
end