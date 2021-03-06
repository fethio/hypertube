alpha= 0.5;        %
mu=0.07;           %
k1=0.296;
k2=6;
phi=pi/6
ll=1;               %long. stretch
A=10;               %inner radii
x=15;               %upper integral boundary
y_0=10;             %initial guess of a
p_i=0.01;           %inflation(inner) pressure

BC = @(a) trr(a,alpha,mu,k1,k2,phi,ll,A,x,p_i);
fplot(BC,[10,12])
value_a = fsolve(BC,y_0)
