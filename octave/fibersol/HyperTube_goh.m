%HyperTube: Derivation of an Analytical Solution to pressurized thick, soft and fiber
%embedded cylinder by ofarukbkaya
syms R a A k1 k2 alpha mu phi c %R:   Reference radial coordinate
                   %A:   Reference inner radii
                   %a:   Current inner radii
                   %k1:  GOH model fiber material property #1
                   %k2:  GOH model fiber material property #2
%phi = pi/6;        %phi: Fiber angle, measured from R-Theta plane
ll = 1;            %ll: lambda_l:   Longitudinal stretch (in Z direction)
f = sqrt(a^2 + (R^2 - A^2)/ll);  %f: r: Current radii
lr = diff(f,R,1)    %lr: lambda_r:   Radial stretch
lt = f/R;          %lt: lambda_theta:  Circumferential Stretch
I4 = ll^2*sin(phi)^2 + lt^2*cos(phi)^2;
                   %4th invariant of deformation tensor C
assume(R,'positive')
assumeAlso(A,'positive')
assumeAlso(a,'positive')
assumeAlso(R-A >= 0)
assumeAlso(a>A)

%the row below (p(R)) is obtained by extracting hydrostatic pressure term p(R)
%from the integration of Cauchy Equation of Motion
syms  p_i                                   %inner pressure
%p_i = 0.01 %MPa

%fplot(ttt_fp(a),[8.5,11.5])

%pause

p = c + (mu/(2*ll))*(log(lt^2*ll)-(a^2*ll-A^2)/(f^2*ll))
    %+ (alpha/(1-alpha))*H_1(A,k1,k2,phi,a,x,ll);
%p_A = subs(p,[mu, A, k1, k2, R],[0.07, 10, 0.296, 65, 10])


trr = (1-alpha)*(-p+mu*lr^2);               %Radial Cauchy Stress
trr_A = subs(trr,[R],[A]) == -p_i;            %BC @ R=A
c_=solve(trr_A,c)
%c_ = subs(solve(trr_A,c),[alpha,mu,A,k1,k2,R,p_i,phi],[0,0.07,10,0.296,65,10,0.01,pi/6])
%syms B                                      %Reference outer radii
%trr_B = @(a) subs(trr,[alpha,mu,A,k1,k2,R,p_i,phi,c],[0,0.07,10,0.296,65,15,0.01,pi/4,c_])
%fplot(trr_B(a),[8.5,9])

%z=subs(trr_B,[a],[11]);
%a_sol = solve(trr_B,a,'Real',true)
%trr_A_ = p_A == (p_i/(alpha - 1)) + - (100*mu)/(- A^2 + a^2 + 100)
%trr_A = subs(trr_A,[p_A],[rhs(trr_A_)])
%p = subs(p,[],[])
