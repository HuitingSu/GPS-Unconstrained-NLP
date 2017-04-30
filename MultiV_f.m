%This is a multi-variate function. Function variables are x(1), x(2),
%x(3), ..., x(n). There is no limitation for n.
function y=MultiV_f(x)   %x y z b
%y=100*(x(2) - x(1)^2)^2 + (1-x(1))^2;
%y=x(1)^2 + (x(2)-5)^2 + (x(3)-4)^4; %+ (x(4)-1)^2+exp(x(4));   %+ 8*x(1)*x(4) + 12* x(2)*x(3);
P=zeros(6,3);
    degrees=[ 50.7	 -118.4
        7.3	 -121.1
        36.5  -82.5
        51.5  -54.7
        16.5  -163.3
        31.3  -19.2];   % longitude & latitude of the satellites  
    G=zeros(1,3);    % Cartesian coordinate of the receiver
    degreesG=[40.424031	-86.910387];  % longitude & latitude of the receiver
    r=26578000;  %orbit radius
    Er=6378000;   %earth radius
    
    % Transform the degree into Cartesian coordinate
    for i=1:6
        P(i,1)=r*cos(degrees(i,1)/360*2*pi)*cos(degrees(i,2)/360*2*pi);
        P(i,2)=r*cos(degrees(i,1)/360*2*pi)*sin(degrees(i,2)/360*2*pi);
        P(i,3)=r*sin(degrees(i,1)/360*2*pi);
    end
        G(1)=Er*cos(degreesG(1)/360*2*pi)*cos(degreesG(2)/360*2*pi);  %[261.7;-4848.3; 4135.7]
        G(2)=Er*cos(degreesG(1)/360*2*pi)*sin(degreesG(2)/360*2*pi);
        G(3)=Er*sin(degreesG(1)/360*2*pi);
    G=[261.7,-4848.3, 4135.7].*1000;    
    % Generate distance from satellites to receiver     
    d=zeros(6,1);
    for i=1:6   % Calculate accurate distance
        d(i) = norm(G-P(i,:));       %sqrt((P(i,1)-G(1))^2+(P(i,2)-G(2))^2+(P(i,3)-G(3))^2);
    end
    b_starc=1.6025; %10*rand()-5;
    Tsc=zeros(6,1);
    for i=1:6
       Tsc(i)=(0-b_starc)-d(i);
    end
    % Calculate total error
    e=0;
    for i=1:6
       e=e+abs( norm(x(1:3)-P(i,:).') - abs(x(4)+Tsc(i)) );
    end
    
    y=e;
end
