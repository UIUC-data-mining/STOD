% Learn the number of components and perform eigen decomposition
% Chi Wang
% chiw@microsoft.com
function [issmalldata, iseigsuccess, isasym, K,U0,D0] = decomp0(...
    dwmat, options)
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

wdmat = dwmat'; 
doclen2 = 1./doclen./(doclen-1);   % the l(l-1) of each document
% diag mat of 1/l(l-1) of each doc
sqdoclen2 = sparse(1:D, 1:D, sqrt(doclen2)); 
punorm2 = sqdoclen2*dwmat;    % the p(w|d)/(l-1) of each doc

% % p(w1,w2)*D, M2=(ALPHA0+1)*(M2 in paper) = 
% (wdmat*punorm2-diag(sum(punorm2)))/D*(ALPHA0+1)
Ew1w2=(punorm2' * punorm2-spdiags(wdmat * doclen2, 0, W, W))/D;
if any(any(Ew1w2 ~= Ew1w2',1),2) == 1
    disp('Ew1w2 asym');
    isasym = true;
    return;
end
nnzEw1w2=nnz(Ew1w2);
%disp(['nonzero # co-occured word pairs = ' num2str(nnzEw1w2)])
upper = max(options.K);
while true  
  [U0,D0,flag]=eigs(Ew1w2,upper,'la'); 
  if flag == 0
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
if iseigsuccess == false
   return;
end

K = selectK(diag(D0), options);
U0=U0(:,1:K);
D0=D0(1:K,1:K);
