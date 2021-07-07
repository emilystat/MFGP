load('Y_simu.mat');

pix = 50;
n_simu=100;

Y_simu((1+pix^2):end,:) = Y_simu((1+pix^2):end,:)./2;

% plot data 

% figure(195);
% for indii=1:100
%     disp(indii)
%     title(['sample ' num2str(indii)]);
%     subplot(1,2,1)
%     imagesc(reshape(Y_simu( 1:pix^2 ,indii),pix,pix));
%     subplot(1,2,2)
%     imagesc(reshape(Y_simu( (pix^2+1):end ,indii),pix,pix)); 
%     title(['sample ' num2str(indii)]);
%     
%     pause(2);
% end

Nn=50;
x = linspace(0,20,Nn);
y = linspace(0,20,Nn);
[lx,ly]=meshgrid(x,y);
locs = [reshape(lx,pix^2,1), reshape(ly,pix^2,1)];
lat = locs(:,1);
lon = locs(:,2);

missing_loc1 = find( locs(:,1) > 2.5 & locs(:,1) < 4.3 & locs(:,2) > 2.5 & locs(:,2) < 10 );
missing_loc2 = find( locs(:,1) > 7 & locs(:,1) < 17 & locs(:,2) > 13 & locs(:,2) < 17 );

missing_blk = [ missing_loc1 ; missing_loc2 ];
remain = setdiff( 1:Nn^2 , missing_blk );
missing_points = randsample( remain, 200 )';
missing_locs_x1 = [missing_loc1; missing_points] ;
missing_locs_x2 = [missing_loc1; missing_loc2; missing_points];

data_o = Y_simu;
missing_all = [missing_locs_x1; missing_locs_x2 + pix^2];
data_o(missing_all, :) = nan;

imagesc(reshape(data_o( 1:pix^2 ,1),pix,pix),'AlphaData',~isnan(reshape(data_o( 1:pix^2 ,1),pix,pix)));
imagesc(reshape(data_o( (1+pix^2):(2*pix^2) ,1),pix,pix),'AlphaData',~isnan(reshape(data_o( (1+pix^2):(2*pix^2) ,1),pix,pix)));


non_missing_x1 = setdiff( 1:Nn^2, missing_locs_x1 );
non_missing_x2 = setdiff( 1:Nn^2, missing_locs_x2 );
non_missing_all = [ (non_missing_x1)'; (pix^2 + non_missing_x2)' ];

n1 = size(non_missing_x1, 2);
n2 = size(non_missing_x2, 2);

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

[missing_x2, missing_y2] = find(isnan(reshape(data_o( (pix^2+1):end ,1 ),pix,pix))); % find all non missing locations
missing_x2 = missing_x2 - 0.5;
missing_y2 = missing_y2 - 0.5;


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
    if mean( d0<radius1) >= 0.6  %% if the basis is too close to any of the missing values, then delete this basis
        record_delete = [record_delete; j];
    end      
    d=((locs(:,1)-basis1(j,1)).^2+(locs(:,2)-basis1(j,2)).^2).^(1/2);  % d is the distance between the points in locs and each points in basis1
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
    if mean( d0<radius1) >= 0.6  %% if the basis is too close to any of the missing values, then delete this basis
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
    if mean( d0<radius1) >= 0.6  %% if the basis is too close to any of the missing values, then delete this basis
        record_delete = [record_delete; j];
    end      
    d=((locs(:,1)-basis1(j,1)).^2+(locs(:,2)-basis1(j,2)).^2).^(1/2);  % d is the distance between the points in locs and each points in basis1
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
    if mean( d0<radius1) >= 0.6  %% if the basis is too close to any of the missing values, then delete this basis
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
    if mean( d0<radius1) >= 0.6  %% if the basis is too close to any of the missing values, then delete this basis
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
    if mean( d0<radius1) >= 0.6  %% if the basis is too close to any of the missing values, then delete this basis
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
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create block matrices for Basis (S). 
S1_blk = [S1_x1 S2_x1];
S1_blk_o = S1_blk(non_missing_x1,:);

S2_blk = [S1_x2 S2_x2];
S2_blk_o = S2_blk(non_missing_x2,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = blkdiag( S1_blk, S2_blk);
S_o = blkdiag( S1_blk_o, S2_blk_o);
[n1, r1] = size(S1_blk_o);
[n2, r2] = size(S2_blk_o);
n = [n1 n2];
r = [r1 r2];


% fit FGP                                                  
dtemp = pdist2(locs, locs);
[indi, indj] = find(dtemp>0 & dtemp<=1); 
H = sparse(indi, indj, ones(length(indi),1), size(locs,1), size(locs,1));              

I = speye(pix^2);                                                
B1 = I( non_missing_x1,: );      
B2 = I( non_missing_x2,: );
Ao = blkdiag(B1, B2);

%%% Create block matrices for Basis (S) for missing predicted observations. 
S1_blk_pred = S1_blk(missing_locs_x1,:);
S2_blk_pred = S2_blk(missing_locs_x2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = blkdiag( S1_blk, S2_blk);
S_pred = blkdiag( S1_blk_pred, S2_blk_pred);
[n1_pred, r1_pred] = size(S1_blk_pred);
[n2_pred, r2_pred] = size(S2_blk_pred);
n_pred = [n1_pred n2_pred];
r_pred = [r1_pred r2_pred];

%%% Create block indentity matrix for missing predicted Xi.
A1_blk_pred = I(missing_locs_x1,:);
A2_blk_pred = I(missing_locs_x2,:);
A_pred = blkdiag(A1_blk_pred, A2_blk_pred);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:n_simu   
    disp(i)
    Y_o = Y_simu(non_missing_all, i);  
    
    x1_trend = mean( Y_simu(1:pix^2, i) );
    x2_trend = mean( Y_simu((pix^2+1):(2*pix^2), i) );    
    x12_trend = [ones(size(non_missing_x1,2),1).*x1_trend; ones(size(non_missing_x2,2),1).*x2_trend ];
    
    resid = Y_o - x12_trend;
    
    
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

    d = variogram( locs(non_missing_x2,:), resid( (size(non_missing_x1,2)+1):end ) );
    a1 = d.distance(1);
    a2 = d.distance(2);
    b1 = d.val(1);
    b2 = d.val(2);
    sig2eps2_new = b1 - (b2-b1)/(a2-a1)*a1;

    if sig2eps2_new < 0
       sig2eps2_new = 1e-9;
    end

    sig2eps_new = [sig2eps1_new sig2eps2_new];

    aprop1=2; aprop2=2;   
    
    Z1=aprop1*resid(1:n1 );
    Z2=aprop2*resid( (n1+1): end );
    Z=[Z1; Z2];       
    sig2eps(1) = aprop1^2*sig2eps_new(1);
    sig2eps(2) = aprop2^2*sig2eps_new(2);
    %%%%%%% EM estimation 
    tic
    [K, gamma, tau2, alpha, iter,epsilon,time]=EM_MCAR_2Vars(S_o, Ao, 2, r, n, Z,H,sig2eps);
    toc

    Q1 = (speye(size(H,1)) - gamma(1)*H)/tau2(1);
    Q2 = (speye(size(H,1)) - gamma(2)*H)/tau2(2);
    P = alpha(1) * speye(size(H,1)) + alpha(2) * H;
    
    Sigma_CAR_inv = [Q1+P'*Q2*P, -P'*Q2; -Q2*P, Q2];
   
    x12_trend_pred = [ones(size(missing_locs_x1,1),1).*x1_trend; ones(size(missing_locs_x2,1),1).*x2_trend ];
      
    % Do the prediction
    aprop = [aprop1 aprop2];
    [Y_pred] = Predict_missings_v2(S_o, S_pred, Ao, A_pred, Z, sig2eps, n, r, n_pred, x12_trend_pred, Sigma_CAR_inv, K, aprop);
    
       
    % Compute CRPS
    Y1_pred(:,i) = Y_pred(1:size(missing_locs_x1,1));
    Y2_pred(:,i) = Y_pred((size(missing_locs_x1,1)+1):end).*2;
    Y1_true = Y_simu(missing_locs_x1,i);
    Y2_true = Y_simu(pix^2 + missing_locs_x2,i).*2;

    [CRPS(i,1)] = crps( [Y1_pred(:,i) Y1_pred(:,i)], Y1_true ); 
    [CRPS(i,2)] = crps( [Y2_pred(:,i) Y2_pred(:,i)], Y2_true ); 

    MSPE(i,1) =   mean( ( Y1_true - Y1_pred(:,i) ).^2 );
    MSPE(i,2) =   mean( ( Y2_true - Y2_pred(:,i) ).^2 );    
    
    [CRPS_block1(i,1)] = crps( [ Y1_pred(1:length(missing_loc1), i ) Y1_pred(1:length(missing_loc1), i ) ], Y1_true(1:length(missing_loc1) ) );
    [CRPS_block1(i,2)] = crps( [ Y2_pred(1:length(missing_loc1), i ) Y2_pred(1:length(missing_loc1), i ) ], Y2_true(1:length(missing_loc1) ) );

    [CRPS_block2(i,1)] = nan;
    [CRPS_block2(i,2)] = crps( [ Y2_pred( (length(missing_loc1)+1):(length(missing_loc1)+length(missing_loc2) ), i )  Y2_pred( (length(missing_loc1)+1):(length(missing_loc1)+length(missing_loc2) ), i ) ], Y2_true( (length(missing_loc1)+1):(length(missing_loc1)+length(missing_loc2) ) ) );

    [CRPS_points(i,1)] = crps( [ Y1_pred( (length(missing_loc1)+1):end, i )  Y1_pred( (length(missing_loc1)+1):end, i ) ], Y1_true( (length(missing_loc1)+1):end ) );
    [CRPS_points(i,2)] = crps( [ Y2_pred( (length(missing_loc1)+length(missing_loc2)+1):end, i )  Y2_pred( (length(missing_loc1)+length(missing_loc2)+1):end, i ) ], Y2_true( (length(missing_loc1)+length(missing_loc2)+1):end ) );

    
    MSPE_block1(i,1) = mean( ( Y1_true(1:length(missing_loc1) ) - Y1_pred(1:length(missing_loc1),i) ).^2 );
    MSPE_block1(i,2) = mean( ( Y2_true(1:length(missing_loc1) ) - Y2_pred(1:length(missing_loc1),i) ).^2 );

    MSPE_block2(i,1) = nan;
    MSPE_block2(i,2) = mean( ( Y2_true( (length(missing_loc1)+1):(length(missing_loc1)+length(missing_loc2) ) ) - Y2_pred( (length(missing_loc1)+1):(length(missing_loc1)+length(missing_loc2) ),i ) ).^2 );

    MSPE_points(i,1) = mean( ( Y1_true((length(missing_loc1)+1):end ) - Y1_pred( (length(missing_loc1)+1):end, i ) ).^2 );
    MSPE_points(i,2) = mean( ( Y2_true((length(missing_loc1)+length(missing_loc2)+1):end ) - Y2_pred( (length(missing_loc1)+length(missing_loc2)+1):end, i ) ).^2 );



    disp(CRPS)
    disp(MSPE)
end  % end for loop
    

mean(CRPS)
mean(MSPE)

save Y1_pred.mat Y1_pred
save Y2_pred.mat Y2_pred
save CRPS.mat CRPS
save MSPE.mat MSPE
save CRPS_block1.mat CRPS_block1
save CRPS_block2.mat CRPS_block2
save CRPS_points.mat CRPS_points
save MSPE_block1.mat MSPE_block1
save MSPE_block2.mat MSPE_block2
save MSPE_points.mat MSPE_points
