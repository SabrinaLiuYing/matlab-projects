
%Thus my personal version of the problem is :
% Using Secant method to implement
% Using the function V(t)=2*t^2+sin(3*t)*exp(-1*t/8)

% the given informations
r1 = 2;
l1 = 3;
l2 = 3;
r2 = 4;
N = 100;
tol = 10^(-4);
a = 0;
b = N;

%volume at given point
VA = 4/3*pi*r1^3;
VB = VA + 1/3*pi*r2^2*l1;
VC = VB + pi*r2^2*l2;
VD = VC + 2/3*pi*r2^3;

%Times, and print times:
TA = slover(@(t) (V(t)-VA),a,b,tol);
TB = slover(@(t) (V(t)-VB),a,b,tol);
TC = slover(@(t) (V(t)-VC),a,b,tol);
TD = slover(@(t) (V(t)-VD),a,b,tol);

%print out each times
TA
TB
TC
TD

%simulate
simulate(r1,l1,l2,r2,N,tol);


%function of simulate
function simulate(r1,l1,l2,r2,N,tol)
    %given information
    H = 2*r1+l1+l2+r2;
    Subint = H/N;
    %volume formula
    Vsp = @(h) (pi/3*(3*r1*h^2-h^3));
    Vcone = @(h) (1/3*pi*(h*r2/l1)^2*h);
    Vcyl = @(h) (pi*r2^2*h);
    Vhs = @(h) (pi/3*(3*r2^2*h-h^3));
    %volume at each point
    VA = Vsp(2*r1);
    VB = VA+Vcone(l1);
    VC = VB+Vcyl(l2);
    %height at each point
    Hsp = 2*r1;
    Hcone = Hsp+l1;
    Hcyl = Hcone+l2;
    %initial T and range
    T = zeros(N);
    a=0;
    b=N;
    
    %hh is the height at hi, tt is t(hi)
    hh= 0;
    tt = 1;
    %loop
    while hh < H
        % when height is less than point A
        if hh <= Hsp
            T(tt)= slover(@(t) (V(t)-Vsp(hh)),a,b,tol);
        % when height is greater than point A but less than point B
        elseif hh <= Hcone
            T(tt)= slover(@(t) (V(t)-Vcone(hh-2*r1)-VA),a,b,tol);
        % when height is greater than point B but less than point C
        elseif hh <= Hcyl
            T(tt)= slover(@(t) (V(t)-Vcyl(hh-2*r1-l1)-VB),a,b,tol);
        % when height is greater than point C but less than point D, (H)
        else
            T(tt)= slover(@(t) (V(t)-Vhs(hh-2*r1-l1-l2)-VC),a,b,tol);
        end
        hh = hh+Subint;
        tt = tt+1;
    end
    %plot the graph
    plot(linspace(0,H,N),T)
    title("Dependency t(h) The time when water reaches level h")
    xlabel("Height(h)")
    ylabel("time(t)")

end




% function of V(t)
function V = V(t)
    V = 2*t^2+sin(3*t)*exp(-1*t/8);
end

% secant method
function x = slover(f, a, b, tol)
    x_0=a;
    x_1=b;
    while abs(x_1-x_0) >= tol
        if abs(f(x_1)-f(x_0)) < tol
            break;
        end
        x_temp=x_1-f(x_1)*((x_1-x_0)/(f(x_1)-f(x_0)));
        x_0 = x_1;
        x_1 = x_temp;
    end
    x=x_1;
end














