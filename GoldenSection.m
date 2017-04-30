%Golden Section Method implementation for minimizing a strictly unimodal function.  --by Huiting Su.
%The tolerance is the distance between two successive points. 
function GoldenSection(a,b,tol)   %(a,b) is the initial searching interval
    tic
    GoldenRatio =(sqrt(5)+1)/2;  %Define Golden Ratio
    c=b-(b-a)/GoldenRatio;  %Find first c&d points
    d=a+(b-a)/GoldenRatio; 
    while abs(c-d)>tol   
        if input_func(c)<input_func(d)
            b=d;
        else
            a=c;
        end
        c=b-(b-a)/GoldenRatio;
        d=a+(b-a)/GoldenRatio;
    end
    display((b+a)/2);
    toc
end

