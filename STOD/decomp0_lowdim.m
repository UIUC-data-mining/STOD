% Learn the number of components and perform eigen decomposition
% use randomized linear algebra for dimensionality reduction 
% Chi Wang
% chiw@microsoft.com
function [issmalldata, iseigsuccess, isasym, K,U0,D0] = ...
    decomp0_lowdim(dwmat, options)
% K - number of topics
% ALPHA0 - summation of alpha_1, ... , alpha_K
% options.K - the range of possible K
% options.N - number of outer iterations, initialization
% options.n - number of inner iterations, fixed point iteration

% folder - the place to store the intermediate data

% D - number of documents
% W - number of words

issmalldata = false;
iseigsuccess = false;
isasym = false;
K=[];
U0=[];
D0=[];

doclen = sum(dwmat,2);    % the length of each document
ind = doclen<3;
dwmat(ind,:)=[];   % filter documents with length<3
doclen(ind,:)=[];
[D,W] = size(dwmat);

if D<10 || nnz(dwmat) < 50
   issmalldata = true;
   return;
end

upper = max(options.K);
R = sparse(randi([0 1], W, upper*3));
wdmat = dwmat'; 
doclen2 = 1./doclen./(doclen-1);   % the l(l-1) of each document
diagdoclen2 = sparse(1:D,1:D,doclen2);
% diag mat of 1/l(l-1) of each doc
sqdoclen2 = sparse(1:D, 1:D, sqrt(doclen2)); 
% punorm2 = sqdoclen2*dwmat;    % the p(w|d)/(l-1) of each doc
% % p(w1,w2)*D, M2=(ALPHA0+1)*(M2 in paper) = 
% (wdmat*punorm2-diag(sum(punorm2)))/D*(ALPHA0+1)
diagc = spdiags(wdmat*doclen2,0,W,W);
Ew1w2R = (wdmat*diagdoclen2*(dwmat*R)-diagc*R)/D;
[U0,~] = svds(Ew1w2R,upper);

Ew1w2U = sqdoclen2*(dwmat*U0);
U = sqrt(diagc)*U0;
Ew1w2 = (Ew1w2U'*Ew1w2U-U'*U)/D;
% Ew1w2=(punorm2' * punorm2-spdiags(wdmat * doclen2, 0, W, W))/D;
if any(any(Ew1w2 ~= Ew1w2',1),2) == 1
    disp('Ew1w2 asym');
    isasym = true;
    return;
end
% nnzEw1w2=nnz(Ew1w2);
% disp(['nonzero # co-occured word pairs = ' num2str(nnzEw1w2)])
[U1,D0]=eig(Ew1w2); 
U1 = U1(:,upper:-1:1);
D0 = diag(D0);
D0 = D0(upper:-1:1);
while true  
  if D0(upper)>0
     iseigsuccess = true;
     break;
  else
     upper = upper - 1;
     if upper == 1
        break;
     end
     disp(['decreasing upper bound to ' int2str(upper)]);
  end
end

% if not learning K, return directly, otherwise select K
if length(options.K)==1
   K = upper;
   return;
end

% if flag has never been 0, return
if ~iseigsuccess
   return;
end

K = selectK(D0(1:upper), options);
U0=U0*U1(:,1:K);
D0=diag(D0(1:K));
