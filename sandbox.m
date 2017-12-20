Nx_prime = 1;
Ny_prime = 1;


if ( ~isprime(Nx_prime) ||  ~isprime(Ny_prime) )
    error('Needs two prime numbers!');
end

ba = zeros(Nx_prime,Ny_prime);

% K is associated with Nx_prime and M is associated with Ny_prime.

% a simple method to implement the equations is to evaluate mod(x^2,r) for
% all x from 1 to r. The resulting values give the locations (I) in Cr
% that contains +1. All other terms in Cr are -1.
Cr = zeros(1,Nx_prime)-1;
cr_idx = unique(sort(mod((1:Nx_prime).^2,Nx_prime)))+1;
Cr(cr_idx) = 1;

Cs = zeros(1,Ny_prime)-1;
cs_idx = unique(sort(mod((1:Ny_prime).^2,Ny_prime)))+1;
Cs(cs_idx) = 1;

for ix = 1:Nx_prime
    for jy = 1:Ny_prime
        if ix == 1
            ba(ix,jy) = 0;
        elseif ( ix ~= 1 && jy == 1 )
            ba(ix,jy) = 1;
        elseif ( Cr(ix)*Cs(jy) == 1 )
            ba(ix,jy) = 1;
        else
            ba(ix,jy) = 0;
        end
    end
end
