%Script to test externally calculated Green's function matrix by comparing it to "inverse-r" Green's function (already implemented in hmmvp).
%Interactions are given by: 1/sqrt(R2+delta).^order, where R2 is distance squared.
%Before running this script, compile C++ code and matlab hmmvp code.

hmmvp_dir='/home/camcat/code/hmmvp';
hmmvp_bin=[hmmvp_dir '/bin/hmmvpbuild_omp'];
hmat_dir='/home/camcat/output/hmmvp-test';

path(path,[hmmvp_dir '/matlab']);

clearvars X G r r2

%Coordinates array:
N=50;
X(1,:)=1:1:N;
X(2:3,:)=0;

order=2;
delta=0.3;

%Parameters common to all cases:
cc.eta=3;
cc.err_method='brem-fro';
cc.tol = 1e-6; %[-6, -9, -12, -16]={very loose, loose, moderate, tight}.
cc.nthreads=3;
cc.command='compress';
cc.allow_overwrite=1;

%Set up structure for "inverse-r"
c0=cc;
c0.write_hmat_filename=[hmat_dir '/inv'];
c0.greens_fn='inverse-r';
c0.X=X;
c0.order=order;
c0.delta=delta;
c0.kvf=[c0.write_hmat_filename '.kvf'];
kvf('Write', c0.kvf, c0, 1);

%Set up Green's function matrix. I add few elements since I don't want it to be square, in case I switched indices in my implementation.
X2=[X X(:,1:5)+X(:,end)];
N2=length(X2)
%Here I pad *receivers*
for i=1:N2
  for j=1:N
     r2(i,j)=sum((X2(:,i)-X(:,j)).^2);
  end
end
r=sqrt(r2+delta);
G=r.^-order;

c1=cc;
c1.write_hmat_filename=[hmat_dir '/ext'];
c1.greens_fn='external';
c1.G=G;
c1.Nr=N2;
c1.Ns=N;
c1.xr=X2;
c1.xs=X;
c1.kvf=[c1.write_hmat_filename '.kvf'];
kvf('Write', c1.kvf, c1, 1);

%Here I switch sources, receivers:
c2=c1;
c2.Nr=c1.Ns;
c2.Ns=c1.Nr;
c2.xr=c1.xs;
c2.xs=c1.xr;
c2.G=c1.G';
c2.write_hmat_filename=[hmat_dir '/ext2'];
c2.kvf=[c2.write_hmat_filename '.kvf'];
kvf('Write', c2.kvf, c2, 1);

%Run hmmvp:
system([hmmvp_bin ' ' c0.kvf]);
system([hmmvp_bin ' ' c1.kvf]);
system([hmmvp_bin ' ' c2.kvf]);

%Read H-matrices and calculate a MVP for a sinusoidal input:
x=sin(X(1,:)/5)';
x2=zeros(N2,1); %Pad sources for c2.
x2(1:N)=x;

[id0 nnz] = hmmvp('init',c0.write_hmat_filename);
y0 = hmmvp('mvp',id0,x);

[id1 nnz] = hmmvp('init',c1.write_hmat_filename);
y1 = hmmvp('mvp',id1,x);

[id2 nnz] = hmmvp('init',c2.write_hmat_filename);
y2 = hmmvp('mvp',id2,x2);

plot(c0.X(1,:),y0); hold on
plot(c1.xr(1,:),y1,'--');
plot(c2.xr(1,:),y2,'o');
title('MVP')
legend({'inverse-r',['ext (' num2str(c1.Nr) 'x' num2str(c1.Ns) ')'],...
                    ['ext (' num2str(c2.Nr) 'x' num2str(c2.Ns) ')']});

