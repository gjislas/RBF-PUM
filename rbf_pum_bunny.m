% Radial kernel
ep=510;  % Shape parameter
phi=@(r2) 1./sqrt(1+ep^2*r2);         % IMQ
dphi=@(r2) -ep^2./sqrt(1+ep^2*r2).^3; % Derivative of IMQ over r

% Weight function for global Shepard partition of unity weighting
wf = @(e,r) r.^4.*(5*spones(r)-4*r);

pu_n = 50; % # of subdomains
wep = 35; % weight parameter
cellradius = 1/wep; % subdomain radius
sub_n = 700; % # of nodes chosen on each subdomain
pts_n = 10000; % global # of nodes

load('Data3D_Bunny3.mat'); % coarser data set 
pu_x = dsites; 
pu_xx = pu_x(:,1); pu_yy = pu_x(:,2); pu_zz = pu_x(:,3);

%basis evaluation matrix for subdomain center candidates
dx = repmat(pu_xx,1,length(pu_xx))-repmat(pu_xx,1,length(pu_xx))';
dy = repmat(pu_yy,1,length(pu_xx))-repmat(pu_yy,1,length(pu_xx))';
dz = repmat(pu_zz,1,length(pu_xx))-repmat(pu_zz,1,length(pu_xx))';
r2 = (dx.^2+dy.^2+dz.^2);
AA = phi(r2);

%column pivoting QR
%RBF centers selection
[QA,RA,EA] = qr(AA,'vector');
[~,~,Ep] = qr(QA(:,1:(pu_n)).','vector');

%subdomain center selection
pu_xx = pu_xx(sort(Ep(1:pu_n))); pu_yy = pu_yy(sort(Ep(1:pu_n))); pu_zz = pu_zz(sort(Ep(1:pu_n)));
cellctrs = [pu_xx pu_yy pu_zz];

% Candidate points
load('Data3D_Bunny0.mat'); %finer data set
lcpoints=[dsites(:,1),dsites(:,2),dsites(:,3)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%local candidate set selection

index = [];
evaltree = KDTreeSearcher(lcpoints);

for j = 1:pu_n
     % Find evaluation points in cell j
   [idx,~] = rangesearch(evaltree,...
        cellctrs(j,:),cellradius);
    eidx = sort(idx{1});
    if (~isempty(eidx))
        xp = lcpoints(eidx,1); yp = lcpoints(eidx,2); 
        zp = lcpoints(eidx,3);
        
        xc = xp; yc = yp; zc = zp;
        x = [xp yp zp];
        
    %basis evaluation matrix for local candidate points
    dx = repmat(xc,1,length(xc))-repmat(xc,1,length(xc))';
    dy = repmat(yc,1,length(xc))-repmat(yc,1,length(xc))';
    dz = repmat(zc,1,length(xc))-repmat(zc,1,length(xc))';
    r2 = (dx.^2+dy.^2+dz.^2);
    AA = phi(r2);

    %ensures # nodes chosen is not less than # of nodes available 
    sub_n_temp = sub_n;
    if sub_n > length(xp)
      sub_n_temp = length(xp);
    else
      sub_n_temp = sub_n;
    end
    
    %column pivoting QR factorizaiton
    %RBF centers selection
    [QA,RA,EA] = qr(AA.','vector');
    [~,~,Ep] = qr(QA(:,1:sub_n_temp).','vector');

    %local node selection
    index = [index eidx(Ep(1:sub_n_temp))];
    end
end

%removes duplicates
cpoints = lcpoints(unique(index),:);
normals = normals(unique(index),:);

% candidate points
xx = cpoints(:,1); yy = cpoints(:,2); zz = cpoints(:,3);

%basis evaluation matrix for candiate points
dx = repmat(xx,1,length(xx))-repmat(xx,1,length(xx))';
dy = repmat(yy,1,length(xx))-repmat(yy,1,length(xx))';
dz = repmat(zz,1,length(xx))-repmat(zz,1,length(xx))';
r2 = (dx.^2+dy.^2+dz.^2);
AA = phi(r2);

%Row normalization
for k = 1:size(AA,1)
  AA(k,:) = AA(k,:)./norm(AA(k,:));
end

%QR factorization with small data set
%[QA,RA,EA] = qr(AA.','vector');
%[~,~,Ep] = qr(QA(:,1:rA).','vector');

%Randomized QR factorization for large data sets
%RBF centers selection
wBinary(AA.','A.bin','double'); %writes basis evaluation matrix to binary 
[m,n] = size(AA.'); 
system(sprintf('./intel_HQRRP_double %d %d A.bin idx.bin',m,n)) %QR in C
Ep = rBinary('idx.bin',1,m,'I'); %reads binary basis evaluation matrix

%Node selection process
wBinary(AA(:,Ep(1:n)).','A.bin','double'); %writes to binary
[m,n] = size(AA(:,Ep(1:n)).'); 
system(sprintf('./intel_HQRRP_double %d %d A.bin idx.bin',m,n)) %QR in C
Ep = rBinary('idx.bin',1,m,'I'); %reads binary basis evaluation matrix
%Node selection
xx = xx(sort(Ep(1:pts_n))); yy = yy(sort(Ep(1:pts_n))); zz = zz(sort(Ep(1:pts_n)));
normals = normals(sort(Ep(1:pts_n)),:);

%Plot of interpolation nodes
figure(1)
plot3(xx,yy,zz,'.')
title('nodes')
saveas(gcf,'nodes.fig')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epoints = [xx yy zz];
evaltree = KDTreeSearcher(epoints);

%% Laplacian computation process

%distance between subdomain centers and points
DM_eval = DistanceMatrixCSRBF(epoints,cellctrs,wep);
%weights for Shepards method
SEM = wf(ep,DM_eval);
SEM = spdiags(1./(SEM*ones(pu_n,1)),0,pts_n,pts_n)*SEM;


LL = zeros(length(epoints(:,1)),length(epoints(:,1)),3);
%%
for j = 1:pu_n
    %Indentify points in each subdomain
    [idx,~] = rangesearch(evaltree,cellctrs(j,:),cellradius);
    
    eidx = sort(idx{1});
    if (~isempty(eidx))
        xp = epoints(eidx,1); yp = epoints(eidx,2); 
        zp = epoints(eidx,3);
        
        x = [xp yp zp];
        
        N=size(xp,1);
        NR=normals(eidx,:);
        
        %compute the surface laplacian
        xij=repmat(x(:,1),[1 N]); xij=xij-xij.'; nxi=repmat(NR(:,1),[1 N]);
        yij=repmat(x(:,2),[1 N]); yij=yij-yij.'; nyi=repmat(NR(:,2),[1 N]);
        zij=repmat(x(:,3),[1 N]); zij=zij-zij.'; nzi=repmat(NR(:,3),[1 N]);
        r2=xij.^2 + yij.^2 + zij.^2; A=dphi(r2);
        DPx=((1-nxi.^2).*xij - nxi.*nyi.*yij - nxi.*nzi.*zij).*A;
        DPy=(-nxi.*nyi.*xij + (1-nyi.^2).*yij - nyi.*nzi.*zij).*A;
        DPz=(-nxi.*nzi.*xij - nyi.*nzi.*yij + (1-nzi.^2).*zij).*A;
        A=phi(r2);
        DPx=(DPx/A); DPy=(DPy/A); DPz=(DPz/A);
        Lap=DPx*DPx+DPy*DPy+DPz*DPz;  %surface Laplacian

      %calls weights for shepards method
      W = eye(length(eidx));
      for k = 1:length(eidx)
          W(k,k) = SEM(eidx(k),j);
      end
    
      %Calls approrite values from data vector 
      %M*V  = V(eidx{j})
      M = zeros(length(eidx),length(epoints(:,1)));
      for k = 1:length(eidx)
          M(k,eidx(k)) = 1;
      end
      
      %Builds vector of appropriate size
      MM = zeros(length(epoints(:,1)),length(eidx));
      for k = 1:length(eidx)
          N(eidx(k),k) = 1;
      end

      LL(:,:,j) = MM*W*Lap*M;
    end
end 
Lap = sum(LL,3);

%Plot of the eigenvalues of the Laplacian matrix
figure(2)
[V,D] = eig(L);
plot(diag(D),'.')
title('eigs')
saveas(gcf,'eigs.fig')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Turing spot pattern using SBDF2
N = length(epoints);
xx = epoints(:,1); yy = epoints(:,2); zz = epoints(:,3);
del=0.00005; d=0.516; tau1=0.02; tau2=0.2; alp=0.899; bet=-0.91; gam=-alp;
ufun=@(u,v) alp*u.*(1-tau1*v.^2)+v.*(1-tau2*u);
vfun=@(u,v) bet*v.*(1+alp*tau1/bet*u.*v)+u.*(gam+tau2*v);

% Initial condition
stream=RandStream('mrg32k3a','seed',7122005); u=rand(stream,N,2)-0.5; 
v=u(:,2); u=u(:,1);

%surface triangulation
[tria]=MyCrustOpen(x);

%Implicit systems for SBDF2
dt=0.1; %time step 
tfinal=300; %final time

%LU decomposition of differentiation matrix
Du=1.5*eye(N)-del*d*dt*(Lap); [Lu,Uu,pu]=lu(Du,'vector');
Dv=1.5*eye(N)-del*dt*(Lap);   [Lv,Uv,pv]=lu(Dv,'vector');
optsu.UT=true; optsl.LT=true;

% One step of backward Euler to bootstrap SBDF2
rhsu=u+dt*ufun(u,v); rhsv=v+dt*vfun(u,v);
u0=u; u=(eye(N)-del*d*dt*Lap)\rhsu;
v0=v; v=(eye(N)-del*dt*Lap)\rhsv;

%time stepping
for j=1:tfinal/dt
   rhsu=2*u-0.5*u0+dt*(2*ufun(u,v)-ufun(u0,v0)); % SBDF2
   rhsv=2*v-0.5*v0+dt*(2*vfun(u,v)-vfun(u0,v0));
   u0=u; u=linsolve(Uu,linsolve(Lu,rhsu(pu),optsl),optsu);
   v0=v; v=linsolve(Uv,linsolve(Lv,rhsv(pv),optsl),optsu);
end

figure(3)
trisurf(tria,xx,yy,zz,u)%plot della superficie
title(sprintf('t=%d.fig',tfinal))
axis equal
axis vis3d
shading interp
colorbar
saveas(gcf,'bunny_spots.fig')
