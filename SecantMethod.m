%Secant Method implementation for finding root of equation by Huiting Su.
function SecantMethod(x0,x1,tol)
    tic        
    while abs(input_func(x1))>tol
        f_x0=input_func(x0);
        f_x1=input_func(x1);
        t=x1-f_x1*(x1-x0)/(f_x1-f_x0);
        x0=x1;
        x1=t;
    end    
    display(x1);
    toc
end
