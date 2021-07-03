load('Y_simu.mat');

pix = 50;
n_simu=100;

Y_simu(2501:5000, :) = Y_simu(2501:5000, :)./2;
Y_simu(5001:7500, :) = Y_simu(5001:7500, :)./4;

% plot data

% figure(195);
% for indii=1:100
%     disp(indii)
%     title(['sample ' num2str(indii)]);
%     subplot(1,3,1)
%     imagesc(reshape(Y_simu( 1:pix^2 ,indii),pix,pix));
%     colorbar;
%     
%     subplot(1,3,2)
%     imagesc(reshape(Y_simu( (pix^2+1):(2*pix^2) ,indii),pix,pix)); 
%     title(['sample ' num2str(indii)]);
%     colorbar;
%      
%     subplot(1,3,3)
%     imagesc(reshape(Y_simu( (2*pix^2+1):end ,indii),pix,pix)); 
%     title(['sample ' num2str(indii)]);
%     colorbar;
%      
%     pause(0.5);
% end

% 
% hist(Y_simu(1:pix^2,1))
% hist(Y_simu( (1+pix^2):end,1))


Nn=50;
x = linspace(0,20,Nn);
y = linspace(0,20,Nn);
[lx,ly]=meshgrid(x,y);
locs = [reshape(lx,pix^2,1), reshape(ly,pix^2,1)];
lat = locs(:,1);
lon = locs(:,2);

missing_loc1 = find( locs(:,1) > 3 & locs(:,1) < 13 & locs(:,2) > 3.5 & locs(:,2) < 6.5 );
missing_loc2 = find( locs(:,1) > 7 & locs(:,1) < 17 & locs(:,2) > 13.5 & locs(:,2) < 16.5 );


missing_blk = [ missing_loc1 ; missing_loc2 ];
remain = setdiff( 1:Nn^2 , missing_blk );
missing_points = randsample( remain, 200 )';
missing_locs_x1 = missing_points;
missing_locs_x2 = [missing_loc1; missing_points];
missing_locs_x3 = [missing_loc1; missing_loc2; missing_points];

data_o = Y_simu;
missing_all = [missing_locs_x1; missing_locs_x2 + pix^2; missing_locs_x3 + 2 * pix^2];
data_o(missing_all, :) = nan;


%figure(196);
%subplot(1,3,1)
%imagesc(reshape(data_o( 1:pix^2 ,1),pix,pix),'AlphaData',~isnan(reshape(data_o( 1:pix^2 ,1),pix,pix)));
%colorbar;
%subplot(1,3,2)
%imagesc(reshape(data_o( (1+pix^2):(2*pix^2) ,1),pix,pix),'AlphaData',~isnan(reshape(data_o( (1+pix^2):(2*pix^2) ,1),pix,pix)));
%colorbar;
%subplot(1,3,3)
%imagesc(reshape(data_o( (1+2*pix^2):end ,1),pix,pix),'AlphaData',~isnan(reshape(data_o( (1+2*pix^2):end ,1),pix,pix)));
%colorbar;

 
non_missing_x1 = setdiff( 1:Nn^2, missing_locs_x1 );
non_missing_x2 = setdiff( 1:Nn^2, missing_locs_x2 );
non_missing_x3 = setdiff( 1:Nn^2, missing_locs_x3 );
non_missing_all = [ (non_missing_x1)'; (pix^2 + non_missing_x2)'; (2*pix^2 + non_missing_x3)' ];

n1 = size(non_missing_x1, 2);
n2 = size(non_missing_x2, 2);
n3 = size(non_missing_x3, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax=0.5:1:49.5;
ay=0.5:1:49.5;
[allx,ally]=meshgrid(ax,ay);
locs=[reshape(allx,pix^2,1), reshape(ally,pix^2,1)];
n = pix^2;

% find positions of missing values of x1, x2.
[missing_x1, missing_y1] = find(isnan(reshape(data_o( 1:pix^2 ,1),pix,pix))); % find all non missing locations
missing_x1 = missing_x1 - 0.5;
missing_y1 = missing_y1 - 0.5;

[missing_x2, missing_y2] = find(isnan(reshape(data_o( (pix^2+1):(2*pix^2) ,1 ),pix,pix))); % find all non missing locations
missing_x2 = missing_x2 - 0.5;
missing_y2 = missing_y2 - 0.5;

[missing_x3, missing_y3] = find(isnan(reshape(data_o( (2*pix^2+1):end ,1 ),pix,pix))); % find all non missing locations
missing_x3 = missing_x3 - 0.5;
missing_y3 = missing_y3 - 0.5;

% define initial basis-function
a=5;
b=5;
crit = 1.5;
cx = linspace(-pix/a/2, pix+pix/a/2, a+2);
cy = linspace(-pix/b/2, pix+pix/b/2, b+2);
cx=cx(2:a+1);
cy=cy(2:b+1);
[cx,cy] = meshgrid(cx, cy);

basis1 = [reshape(cx, a*b, 1), reshape(cy, a*b, 1)];
d = pdist(basis1);
radius1 = crit*min(d);
	
r1 = a*b;
S1_x1 = sparse(n, r1);
record_delete = [];
for j = 1:r1    
    d0 = ((missing_x1-basis1(j,1)).^2+(missing_y1-basis1(j,2)).^2).^(1/2);
    if mean( d0<radius1 ) >= 0.6  %% if the basis is too close to any of the missing values, then delete this basis
        record_delete = [record_delete; j];
    end      
    d=((locs(:,1)-basis1(j,1)).^2+(locs(:,2)-basis1(j,2)).^2).^(1/2);  %d is the distance between the points in locs and each points in basis1
	S1_x1(:, j) = (1-(d./radius1).^2).^2.*(d<=radius1);             
end                                                                     

S1_x1( :, record_delete ) = [];                                                                 

%figure(1001);
%for j=1:r1
%     imagesc(reshape(S1_x1(:,j),pix,pix));colorbar;
%     pause(.1);
%end

a=a*2;
b=b*2;
cx = linspace(-pix/a/2, pix+pix/a/2, a+2);
cy = linspace(-pix/b/2, pix+pix/b/2, b+2);
cx=cx(2:a+1);
cy=cy(2:b+1);
[cx,cy] = meshgrid(cx, cy);

basis1 = [reshape(cx, a*b, 1), reshape(cy, a*b, 1)];
d = pdist(basis1);
radius1 = crit*min(d);
	
r1 = a*b;
S2_x1 = sparse(n, r1);
record_delete = [];
for j = 1:r1    
    d0 = ((missing_x1-basis1(j,1)).^2+(missing_y1-basis1(j,2)).^2).^(1/2);
    if mean( d0<radius1 ) >= 0.6  %% if the basis is too close to any of the missing values, then delete this basis
        record_delete = [record_delete; j];
    end      
    d=((locs(:,1)-basis1(j,1)).^2+(locs(:,2)-basis1(j,2)).^2).^(1/2);  %d is the distance between the points in locs and each points in basis1
	S2_x1(:, j) = (1-(d./radius1).^2).^2.*(d<=radius1);              
end                                                                     

S2_x1( :, record_delete ) = [];


%figure(1001);
%for j=1:r1
%     imagesc(reshape(S2_x1(:,j),pix,pix));colorbar;
%     pause(.1);
%end

a=a*2;
b=b*2;
cx = linspace(-pix/a/2, pix+pix/a/2, a+2);
cy = linspace(-pix/b/2, pix+pix/b/2, b+2);
cx=cx(2:a+1);
cy=cy(2:b+1);
[cx,cy] = meshgrid(cx, cy);

basis1 = [reshape(cx, a*b, 1), reshape(cy, a*b, 1)];
d = pdist(basis1);
radius1 = crit*min(d);
	
r1 = a*b;
S3_x1 = sparse(n, r1);
record_delete = [];
for j = 1:r1    
    d0 = ((missing_x1-basis1(j,1)).^2+(missing_y1-basis1(j,2)).^2).^(1/2);
    if mean( d0<radius1 ) >= 0.6  %% if the basis is too close to any of the missing values, then delete this basis
        record_delete = [record_delete; j];
    end      
    d=((locs(:,1)-basis1(j,1)).^2+(locs(:,2)-basis1(j,2)).^2).^(1/2);  %d is the distance between the points in locs and each points in basis1
	S3_x1(:, j) = (1-(d./radius1).^2).^2.*(d<=radius1);              
end                                                                    

S3_x1( :, record_delete ) = [];


%  figure(1001);
%  for j=1:r1
%      imagesc(reshape(S3_x1(:,j),pix,pix));colorbar;
%      pause(.001);
%  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define initial basis-function for X2
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=5;
b=5;
crit = 1.5;
cx = linspace(-pix/a/2, pix+pix/a/2, a+2);
cy = linspace(-pix/b/2, pix+pix/b/2, b+2);
cx=cx(2:a+1);
cy=cy(2:b+1);
[cx,cy] = meshgrid(cx, cy);

basis1 = [reshape(cx, a*b, 1), reshape(cy, a*b, 1)];
d = pdist(basis1);
radius1 = crit*min(d);
	
r1 = a*b;
S1_x2 = sparse(n, r1);
record_delete = [];
for j = 1:r1    
    d0 = ((missing_x2-basis1(j,1)).^2+(missing_y2-basis1(j,2)).^2).^(1/2);
    if mean( d0<radius1 ) >= 0.6  %% if the basis is too close to any of the missing values, then delete this basis
        record_delete = [record_delete; j];
    end      
    d=((locs(:,1)-basis1(j,1)).^2+(locs(:,2)-basis1(j,2)).^2).^(1/2);  
	S1_x2(:, j) = (1-(d./radius1).^2).^2.*(d<=radius1);              
end                                                                     

S1_x2( :, record_delete ) = [];                                                                  

% figure(1001);
%  for j=1:r1
%      imagesc(reshape(S1_x2(:,j),pix,pix));colorbar;
%      pause(.001);
%  end


a=a*2;
b=b*2;
cx = linspace(-pix/a/2, pix+pix/a/2, a+2);
cy = linspace(-pix/b/2, pix+pix/b/2, b+2);
cx=cx(2:a+1);
cy=cy(2:b+1);
[cx,cy] = meshgrid(cx, cy);

basis1 = [reshape(cx, a*b, 1), reshape(cy, a*b, 1)];
d = pdist(basis1);
radius1 = crit*min(d);
	
r1 = a*b;
S2_x2 = sparse(n, r1);
record_delete = [];
for j = 1:r1    
    d0 = ((missing_x2-basis1(j,1)).^2+(missing_y2-basis1(j,2)).^2).^(1/2);
    if mean( d0<radius1 ) >= 0.6  %% if the basis is too close to any of the missing values, then delete this basis
        record_delete = [record_delete; j];
    end      
    d=((locs(:,1)-basis1(j,1)).^2+(locs(:,2)-basis1(j,2)).^2).^(1/2);  
	S2_x2(:, j) = (1-(d./radius1).^2).^2.*(d<=radius1);              
end                                                                    

S2_x2( :, record_delete ) = [];

%figure(1001);
%  for j=1:r1
%      imagesc(reshape(S2_x2(:,j),pix,pix));colorbar;
%      pause(.001);
%  end


a=a*2;
b=b*2;
cx = linspace(-pix/a/2, pix+pix/a/2, a+2);
cy = linspace(-pix/b/2, pix+pix/b/2, b+2);
cx=cx(2:a+1);
cy=cy(2:b+1);
[cx,cy] = meshgrid(cx, cy);

basis1 = [reshape(cx, a*b, 1), reshape(cy, a*b, 1)];
d = pdist(basis1);
radius1 = crit*min(d);
	
r1 = a*b;
S3_x2 = sparse(n, r1);
record_delete = [];
for j = 1:r1    
    d0 = ((missing_x2-basis1(j,1)).^2+(missing_y2-basis1(j,2)).^2).^(1/2);
    if mean( d0<radius1 ) >= 0.6  %% if the basis is too close to any of the missing values, then delete this basis
        record_delete = [record_delete; j];
    end      
    d=((locs(:,1)-basis1(j,1)).^2+(locs(:,2)-basis1(j,2)).^2).^(1/2);  
	S3_x2(:, j) = (1-(d./radius1).^2).^2.*(d<=radius1);              
end                                                                     

S3_x2( :, record_delete ) = [];

%figure(1001);
%  for j=1:r1
%      imagesc(reshape(S3_x2(:,j),pix,pix));colorbar;
%      pause(.001);
%  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define initial basis-function for X3
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=5;
b=5;
crit = 1.5;
cx = linspace(-pix/a/2, pix+pix/a/2, a+2);
cy = linspace(-pix/b/2, pix+pix/b/2, b+2);
cx=cx(2:a+1);
cy=cy(2:b+1);
[cx,cy] = meshgrid(cx, cy);

basis1 = [reshape(cx, a*b, 1), reshape(cy, a*b, 1)];
d = pdist(basis1);
radius1 = crit*min(d);
	
r1 = a*b;
S1_x3 = sparse(n, r1);
record_delete = [];
for j = 1:r1    
    d0 = ((missing_x3-basis1(j,1)).^2+(missing_y3-basis1(j,2)).^2).^(1/2);
    if mean( d0<radius1 ) >= 0.6  %% if the basis is too close to any of the missing values, then delete this basis
        record_delete = [record_delete; j];
    end      
    d=((locs(:,1)-basis1(j,1)).^2+(locs(:,2)-basis1(j,2)).^2).^(1/2);  
	S1_x3(:, j) = (1-(d./radius1).^2).^2.*(d<=radius1);              
end                                                                     

S1_x3( :, record_delete ) = [];                                                                  

% figure(1001);
%  for j=1:r1
%      imagesc(reshape(S1_x2(:,j),pix,pix));colorbar;
%      pause(.001);
%  end


a=a*2;
b=b*2;
cx = linspace(-pix/a/2, pix+pix/a/2, a+2);
cy = linspace(-pix/b/2, pix+pix/b/2, b+2);
cx=cx(2:a+1);
cy=cy(2:b+1);
[cx,cy] = meshgrid(cx, cy);

basis1 = [reshape(cx, a*b, 1), reshape(cy, a*b, 1)];
d = pdist(basis1);
radius1 = crit*min(d);
	
r1 = a*b;
S2_x3 = sparse(n, r1);
record_delete = [];
for j = 1:r1    
    d0 = ((missing_x3-basis1(j,1)).^2+(missing_y3-basis1(j,2)).^2).^(1/2);
    if mean( d0<radius1 ) >= 0.6  %% if the basis is too close to any of the missing values, then delete this basis
        record_delete = [record_delete; j];
    end      
    d=((locs(:,1)-basis1(j,1)).^2+(locs(:,2)-basis1(j,2)).^2).^(1/2);  
	S2_x3(:, j) = (1-(d./radius1).^2).^2.*(d<=radius1);              
end                                                                    

S2_x3( :, record_delete ) = [];

%figure(1001);
%  for j=1:r1
%      imagesc(reshape(S2_x2(:,j),pix,pix));colorbar;
%      pause(.001);
%  end


a=a*2;
b=b*2;
cx = linspace(-pix/a/2, pix+pix/a/2, a+2);
cy = linspace(-pix/b/2, pix+pix/b/2, b+2);
cx=cx(2:a+1);
cy=cy(2:b+1);
[cx,cy] = meshgrid(cx, cy);

basis1 = [reshape(cx, a*b, 1), reshape(cy, a*b, 1)];
d = pdist(basis1);
radius1 = crit*min(d);
	
r1 = a*b;
S3_x3 = sparse(n, r1);
record_delete = [];
for j = 1:r1    
    d0 = ((missing_x3-basis1(j,1)).^2+(missing_y3-basis1(j,2)).^2).^(1/2);
    if mean( d0<radius1 ) >= 0.6  %% if the basis is too close to any of the missing values, then delete this basis
        record_delete = [record_delete; j];
    end      
    d=((locs(:,1)-basis1(j,1)).^2+(locs(:,2)-basis1(j,2)).^2).^(1/2);  
	S3_x3(:, j) = (1-(d./radius1).^2).^2.*(d<=radius1);              
end                                                                     

S3_x3( :, record_delete ) = [];

%figure(1001);
%  for j=1:r1
%      imagesc(reshape(S3_x2(:,j),pix,pix));colorbar;
%      pause(.001);
%  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create block matrices for Basis (S). 
S1_blk = [S1_x1 S2_x1];
S1_blk_o = S1_blk(non_missing_x1,:);

S2_blk = [S1_x2 S2_x2];
S2_blk_o = S2_blk(non_missing_x2,:);

S3_blk = [S1_x3 S2_x3];
S3_blk_o = S3_blk(non_missing_x3,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = blkdiag( S1_blk, S2_blk, S3_blk);
S_o = blkdiag( S1_blk_o, S2_blk_o, S3_blk_o);
[n1, r1] = size(S1_blk_o);
[n2, r2] = size(S2_blk_o);
[n3, r3] = size(S3_blk_o);
n = [n1 n2 n3];
r = [r1 r2 r3];


% fit FGP                                                  
dtemp = pdist2(locs, locs);
[indi, indj] = find(dtemp>0 & dtemp<=1); 
H = sparse(indi, indj, ones(length(indi),1), size(locs,1), size(locs,1));             

I = speye(pix^2);                                                
B1 = I( non_missing_x1,: );      
B2 = I( non_missing_x2,: );
B3 = I( non_missing_x3,: );
Ao = blkdiag(B1, B2, B3);

%%% Create block matrices for Basis (S) for missing predicted observations. 
S1_blk_pred = S1_blk(missing_locs_x1,:);
S2_blk_pred = S2_blk(missing_locs_x2,:);
S3_blk_pred = S3_blk(missing_locs_x3,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = blkdiag( S1_blk, S2_blk, S3_blk);
S_pred = blkdiag( S1_blk_pred, S2_blk_pred, S3_blk_pred);
[n1_pred, r1_pred] = size(S1_blk_pred);
[n2_pred, r2_pred] = size(S2_blk_pred);
[n3_pred, r3_pred] = size(S3_blk_pred);
n_pred = [n1_pred n2_pred n3_pred];
r_pred = [r1_pred r2_pred r3_pred];

%%% Create block indentity matrix for missing predicted Xi.
A1_blk_pred = I(missing_locs_x1,:);
A2_blk_pred = I(missing_locs_x2,:);
A3_blk_pred = I(missing_locs_x3,:);
A_pred = blkdiag(A1_blk_pred, A2_blk_pred, A3_blk_pred);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
for i=1:10
    disp(i)
    Y_o = Y_simu(non_missing_all, i);  
    
    x1_trend = mean( Y_simu(1:pix^2, i) );
    x2_trend = mean( Y_simu((pix^2+1):(2*pix^2), i) ); 
    x3_trend = mean( Y_simu((2*pix^2+1):(3*pix^2), i) );  
    
    x123_trend = [ones(size(non_missing_x1,2),1).*x1_trend; ones(size(non_missing_x2,2),1).*x2_trend; ones(size(non_missing_x3,2),1).*x3_trend ];
    
    resid = Y_o - x123_trend;
    
    
    %%% Use variogram to estimate sig2eps
    d = variogram( locs(non_missing_x1,:), resid(1:length(non_missing_x1) ) );
    a1 = d.distance(1);
    a2 = d.distance(2);
    b1 = d.val(1);
    b2 = d.val(2);
    sig2eps1_new = b1 - (b2-b1)/(a2-a1)*a1;

    if sig2eps1_new < 0
       sig2eps1_new = 1e-9;
    end

    d = variogram( locs(non_missing_x2,:), resid( (size(non_missing_x1,2)+1):( size(non_missing_x1,2)+size(non_missing_x2,2) ) ) );
    a1 = d.distance(1);
    a2 = d.distance(2);
    b1 = d.val(1);
    b2 = d.val(2);
    sig2eps2_new = b1 - (b2-b1)/(a2-a1)*a1;

    if sig2eps2_new < 0
       sig2eps2_new = 1e-9;
    end
    
    
    d = variogram( locs(non_missing_x3,:), resid( ( size(non_missing_x1,2)+size(non_missing_x2,2)+1 ):end ) );
    a1 = d.distance(1);
    a2 = d.distance(2);
    b1 = d.val(1);
    b2 = d.val(2);
    sig2eps3_new = b1 - (b2-b1)/(a2-a1)*a1;

    if sig2eps3_new < 0
       sig2eps3_new = 1e-9;
    end
    
    

    sig2eps_new = [sig2eps1_new sig2eps2_new sig2eps3_new];

    aprop1=3.5; aprop2=2.5; aprop3=3.5;  
    aprop = [aprop1; aprop2; aprop3];
    Z1=aprop1*resid(1:n1 );
    Z2=aprop2*resid( (n1+1):(n1+n2) );
    Z3=aprop3*resid( (n1+n2+1): end );
    
    Z=[Z1; Z2; Z3];       
    sig2eps(1) = aprop1.^2.*sig2eps_new(1);
    sig2eps(2) = aprop2.^2.*sig2eps_new(2);
    sig2eps(3) = aprop3.^2.*sig2eps_new(3);
    
    %%%%%%% EM estimation 
    tic
    [K, gamma, tau2, alpha, iter,epsilon,time]=EM_MCAR_3Vars(S_o, Ao, 3, r, n, Z,H,sig2eps);
    toc

    Q1 = (speye(size(H,1)) - gamma(1)*H)/tau2(1);
    Q2 = (speye(size(H,1)) - gamma(2)*H)/tau2(2);
    Q3 = (speye(size(H,1)) - gamma(3)*H)/tau2(3);
    P1 = alpha(1) * speye(size(H,1)) + alpha(2) * H;
    P2 = alpha(3) * speye(size(H,1)) + alpha(4) * H;
    P3 = alpha(5) * speye(size(H,1)) + alpha(6) * H;


    Sigma_CAR_inv = [ Q1+P1'*Q2*P1+P2'*Q3*P2 , P2'*Q3*P3-P1'*Q2, -P2'*Q3;
                      (P2'*Q3*P3-P1'*Q2)' ,      Q2+P3'*Q3*P3,   -P3'*Q3;
                      -(P2'*Q3)' ,              -(P3'*Q3)',          Q3 ];

   
    x123_trend_pred = [ones(size(missing_locs_x1,1),1).*x1_trend; ones(size(missing_locs_x2,1),1).*x2_trend; ones(size(missing_locs_x3,1),1).*x3_trend ];
      
    % Do the prediction
    [Y_pred] = Predict_missings_v3(S_o, S_pred, Ao, A_pred, Z, sig2eps, n, r, n_pred, x123_trend_pred, Sigma_CAR_inv, K, aprop);
       
    % Compute CRPS
    Y1_pred(:,i) = Y_pred(1:size(missing_locs_x1,1));
    Y2_pred(:,i) = Y_pred((size(missing_locs_x1,1)+1):(size(missing_locs_x1,1)+size(missing_locs_x2,1) ) ).*2;
    Y3_pred(:,i) = Y_pred((size(missing_locs_x1,1)+size(missing_locs_x2,1)+1):end).*4;
    
    Y1_true = Y_simu(missing_locs_x1,i);
    Y2_true = Y_simu(pix^2 + missing_locs_x2,i).*2;
    Y3_true = Y_simu(2*pix^2 + missing_locs_x3,i).*4;

    [CRPS(i,1)] = crps( [Y1_pred(:,i) Y1_pred(:,i)], Y1_true ); 
    [CRPS(i,2)] = crps( [Y2_pred(:,i) Y2_pred(:,i)], Y2_true ); 
    [CRPS(i,3)] = crps( [Y3_pred(:,i) Y3_pred(:,i)], Y3_true ); 

    MSPE(i,1) =   mean( ( Y1_true - Y1_pred(:,i) ).^2 );
    MSPE(i,2) =   mean( ( Y2_true - Y2_pred(:,i) ).^2 );
    MSPE(i,3) =   mean( ( Y3_true - Y3_pred(:,i) ).^2 );
      
    [CRPS_block1(i,1)] = nan;
    [CRPS_block1(i,2)] = crps( [ Y2_pred(1:length(missing_loc1), i ) Y2_pred(1:length(missing_loc1), i ) ], Y2_true(1:length(missing_loc1) ) );
    [CRPS_block1(i,3)] = crps( [ Y3_pred(1:length(missing_loc1), i ) Y3_pred(1:length(missing_loc1), i ) ], Y3_true(1:length(missing_loc1) ) );

    [CRPS_block2(i,1)] = nan;
    [CRPS_block2(i,2)] = nan;
    [CRPS_block2(i,3)] = crps( [ Y3_pred( (length(missing_loc1)+1):(length(missing_loc1)+length(missing_loc2) ), i )  Y3_pred( (length(missing_loc1)+1):(length(missing_loc1)+length(missing_loc2) ), i ) ], Y3_true( (length(missing_loc1)+1):(length(missing_loc1)+length(missing_loc2) ) ) );

    [CRPS_points(i,1)] = crps( [ Y1_pred( :, i )  Y1_pred( :, i ) ], Y1_true );
    [CRPS_points(i,2)] = crps( [ Y2_pred( (length(missing_loc1)+1):end, i )  Y2_pred( (length(missing_loc1)+1):end, i ) ], Y2_true( (length(missing_loc1)+1):end ) );
    [CRPS_points(i,3)] = crps( [ Y3_pred( (length(missing_loc1)+length(missing_loc2)+1):end, i )  Y3_pred( (length(missing_loc1)+length(missing_loc2)+1):end, i ) ], Y3_true( (length(missing_loc1)+length(missing_loc2)+1):end ) );

    
    MSPE_block1(i,1) = nan;
    MSPE_block1(i,2) = mean( ( Y2_true(1:length(missing_loc1) ) - Y2_pred(1:length(missing_loc1),i) ).^2 );
    MSPE_block1(i,3) = mean( ( Y3_true(1:length(missing_loc1) ) - Y3_pred(1:length(missing_loc1),i) ).^2 );

    MSPE_block2(i,1) = nan;
    MSPE_block2(i,2) = nan;
    MSPE_block2(i,3) = mean( ( Y3_true( (length(missing_loc1)+1):(length(missing_loc1)+length(missing_loc2) ) ) - Y3_pred( (length(missing_loc1)+1):(length(missing_loc1)+length(missing_loc2) ),i ) ).^2 );

    MSPE_points(i,1) = mean( ( Y1_true - Y1_pred(:, i ) ).^2 );
    MSPE_points(i,2) = mean( ( Y2_true((length(missing_loc1)+1):end ) - Y2_pred( (length(missing_loc1)+1):end, i ) ).^2 );
    MSPE_points(i,3) = mean( ( Y3_true((length(missing_loc1)+length(missing_loc2)+1):end ) - Y3_pred( (length(missing_loc1)+length(missing_loc2)+1):end, i ) ).^2 );

    disp(MSPE)


end  % end for loop
toc

mean(CRPS)
mean(MSPE)

save Y1_pred.mat Y1_pred
save Y2_pred.mat Y2_pred
save Y3_pred.mat Y3_pred
save CRPS.mat CRPS
save MSPE.mat MSPE
save CRPS_block1.mat CRPS_block1
save CRPS_block2.mat CRPS_block2
save CRPS_points.mat CRPS_points
save MSPE_block1.mat MSPE_block1
save MSPE_block2.mat MSPE_block2
save MSPE_points.mat MSPE_points