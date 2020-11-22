function sig_rr = trr(a,alpha,mu,k1,k2,phi,ll,A,x,p_i)

c = -(p_i + ((mu*(log(a^2/A^2) + (A^2 - a^2)/a^2))/2 -...
	     (A^2*mu)/a^2)*(alpha - 1))/(alpha - 1);

int_m = (mu/(2*ll))*(log((sqrt(a^2 + (x^2 - A^2)/ll)/x)^2*ll)-...
    (a^2*ll-A^2)/(sqrt(a^2 + (x^2 - A^2)/ll)^2*ll));

p =  c + int_m + alpha/(1-alpha)*H_1(a,k1,k2,phi,ll,A,x);

sig_rr =  (1-alpha)*(-p + mu*(x/(ll*(a^2 - (A^2 - x^2)/ll)^(1/2)))^2);
