pix = 50;
n_simu = 150;

Nn=50;
x = linspace(0,20,Nn);
y = linspace(0,20,Nn);
[lx,ly]=meshgrid(x,y);
locs = [reshape(lx,pix^2,1), reshape(ly,pix^2,1)];
lat = locs(:,1);
lon = locs(:,2);
dtemp = pdist2(locs, locs);

sig2_11 = 1; sig2_2_1 = 0.5; 
kappa_11 = 1; kappa_2_1 = 3; 
nu_11 = 1.5; nu_2_1 = 0.5;

% eta0: unit length for each dimension
eta0 = (20-0)/50;

% Note: diag(dtemp)=0, causing diag(C11)=diag(C_2_1)=0, need to adjust
dtemp_adj = dtemp;
for i=1:pix^2
    dtemp_adj(i,i) = 0.00001;   
end

% C11(s,u)
Bessel_11 = besselk(nu_11,kappa_11 .* dtemp_adj); 
C11 = sig2_11/( 2^(nu_11 - 1)*gamma(nu_11) ) .* ( kappa_11 .* dtemp_adj ).^nu_11 .* Bessel_11;

% C2|1(s,u)
Bessel_2_1 = besselk(nu_2_1,kappa_2_1 .* dtemp_adj);
C_2_1 = sig2_2_1/( 2^(nu_2_1 - 1)*gamma(nu_2_1) ) .* ( kappa_2_1 .* dtemp_adj ).^nu_2_1 .* Bessel_2_1;


Delta = [0 0]; r = 50; A = 0.1;
dtemp2 = pdist2(locs, locs+Delta);
dtemp_s_to_0 = pdist2(locs, zeros(2500,2));

% B matrix, in which each element is b(s,v)
B = zeros(pix^2, pix^2);
[locx, locy] = find(dtemp2 < r/20 | dtemp2 == r/20);
store = A .* ( 1.-( (dtemp_s_to_0 - dtemp2 )/r ).^2 ).^2;

for i = 1:size(locx)
    disp(i)
    B(locx(i), locy(i)) = store(locx(i), locy(i)) ;
end

% generate Y1
rng(1001)
Y1 = mvnrnd(zeros(pix^2,1), C11, n_simu )';
Y2 = B * Y1 .* (eta0^2);


mean(var(Y1))*0.05
mean(var(Y2))*0.05

rng(2002)
eps1 = mvnrnd(zeros(pix^2,1), 0.05*eye(pix^2), n_simu)';
eps2 = mvnrnd(zeros(pix^2,1), 0.04*eye(pix^2), n_simu)';

Y_simu = [Y2; Y1] + [eps2; eps1];

% plot data 
figure(195);
for indii=1:100
    disp(indii)
    title(['sample ' num2str(indii)]);
    subplot(1,2,1)
    imagesc(reshape(Y_simu( 1:pix^2 ,indii),pix,pix));
    colorbar;
    subplot(1,2,2)
    imagesc(reshape(Y_simu( (pix^2+1):end ,indii),pix,pix)); 
    colorbar;
    title(['sample ' num2str(indii)]);
    
    pause(1);
end

save Y_simu.mat Y_simu

hist(Y_simu(1:2500,1))
hist(Y_simu(2501:end,1))
mean(Y_simu(1:2500,1))
mean(Y_simu(2501:end,1))
