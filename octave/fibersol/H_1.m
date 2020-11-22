function value = H_1(a,k1,k2,phi,ll,A,x)

  f = @(R) sqrt(a^2 + (R^2 - A^2)/ll);  %f: r: Current radii
  I4 = @(R) ll^2*sin(phi)^2 + ((sqrt(a^2 + (R^2 - A^2)/ll))/R)^2*cos(phi)^2;

H_1_= @(R) -cos(phi)^2*(4*k1*(sqrt(a^2 + (R^2 - A^2)/ll))*...
    (ll^2*sin(phi)^2 + ((sqrt(a^2 + (R^2 - A^2)/ll))/R)^2*cos(phi)^2-1)*...
    exp(k2*(ll^2*sin(phi)^2 + ((sqrt(a^2 + (R^2 - A^2)/ll))/R)^2*cos(phi)^2-1)^2)/R^2);
# H_1_(x);
%H_1__= @(R) H_1_(R,y_0);

value = Romberg(H_1_,A,x);
