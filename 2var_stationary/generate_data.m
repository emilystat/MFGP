pix = 50;
n_simu = 100;

Nn=50;
x = linspace(0,20,Nn);
y = linspace(0,20,Nn);
[lx,ly]=meshgrid(x,y);
locs = [reshape(lx,pix^2,1), reshape(ly,pix^2,1)];
lat = locs(:,1);
lon = locs(:,2);
dtemp = pdist2(locs, locs);

% Note: diag(dtemp)=0, causing diag(C11)=NAN, need to adjust
dtemp_adj = dtemp;
for i=1:pix^2
    dtemp_adj(i,i) = 0.00001;   
end


p = 2; alpha = 1; 
nu = zeros(p,p);
beta = zeros(p,p);
sigma = zeros(p,p);

for j=1:2
    for k = 1:2
        nu(j,k) = 0.5 + 0.5*(6-j-k)/(2*p-2);
        beta(j,k) = 0.8^(abs(j-k));
    end    
end

for j=1:2
    for k = 1:2
        sigma(j,k) = j * k * gamma(nu(j,k))/gamma(nu(j,k)+1) * sqrt( gamma(nu(j,j)+1) * gamma(nu(k,k)+1) / gamma(nu(j,j))/gamma(nu(k,k)) ) * beta(j,k);      
    end
end

Bessel_11 = besselk(nu(1,1), dtemp_adj .* alpha); 
Bessel_22 = besselk(nu(2,2), dtemp_adj .* alpha);
Bessel_12 = besselk(nu(1,2), dtemp_adj .* alpha);

C11 = sigma(1,1).* ( (dtemp_adj.*alpha).^nu(1,1) ) ./ ( 2^(nu(1,1)-1)*gamma(nu(1,1)) ) .* Bessel_11;
C22 = sigma(2,2).* ( (dtemp_adj.*alpha).^nu(2,2) ) ./ ( 2^(nu(2,2)-1)*gamma(nu(2,2)) ) .* Bessel_22;
C12 = sigma(1,2).* ( (dtemp_adj.*alpha).^nu(1,2) ) ./ ( 2^(nu(1,2)-1)*gamma(nu(1,2)) ) .* Bessel_12;


% generate Y1,Y2,Y3 together
Sigma = [C11, C12; C12', C22];

rng(1002)
Y = mvnrnd(zeros(p*pix^2,1), Sigma, n_simu )';

Sigma_eps = blkdiag( 0.05 * eye(pix^2), 0.2 * eye(pix^2) );
eps = mvnrnd(zeros(p*pix^2,1), Sigma_eps, n_simu)';

Y_simu = Y + eps;

hist(Y(1:2500,1))
hist(Y(2501:5000,1))
hist(eps(1:2500,1))
hist(eps(2501:5000,1))


% plot data 

figure(195);
for indii=1:100
    disp(indii)
    title(['sample ' num2str(indii)]);
    subplot(1,2,1)
    imagesc(reshape(Y_simu( 1:pix^2 ,indii),pix,pix));
    colorbar;
    
    subplot(1,2,2)
    imagesc(reshape(Y_simu( (pix^2+1):(2*pix^2) ,indii),pix,pix)); 
    title(['sample ' num2str(indii)]);
    colorbar;
     
     
    pause(1);
end



save Y_simu.mat Y_simu
