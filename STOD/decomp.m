% Scalable Tensor Orthogonal Decomposition for LDA
% Chi Wang
% chiw@microsoft.com
function [isnegeig, wtmat,ALPHA] = decomp(dwmat,k,ALPHA0, U0,D0,options)
% dwmat - sparse matrix of document-word matrix
% k - number of topics
% ALPHA0 - summation of alpha_1, ... , alpha_T
% options.N
% options.n
% options.lr

% Output:
% isnegeig: return true if there is a negative eigen value, 
% so this alpha cannot be selected
% wtmat: word-topic matrix

% D - number of documents
% W - number of words

% initialize return value;
isnegeig = false;
wtmat = [];
ALPHA = [];

doclen = sum(dwmat,2);    % the length of each document
ind = doclen<3;
dwmat(ind,:)=[];   % filter documents with length<3
doclen(ind,:)=[];
[D,W] = size(dwmat);

wdmat = dwmat'; 
doclen2 = 1./doclen./(doclen-1);   % the l(l-1) of each document
doclen3 = doclen2./(doclen-2);   % the l(l-1)(l-2) of each doc
doclen = sparse(1:D,1:D,1./doclen); % diag mat of 1/l of each doc
% diag mat of 1/l(l-1) of each doc
M1 = full(sum(doclen*dwmat, 1)')/D;  
x = sqrt(ALPHA0) * U0'*M1;
DD = D0 - x*x'; 
%U0^T*M2*U0=D0-ALPHA0*(U0'*M1)*(M1'*U0)

[U1,D1]=eig(DD); %U0^T*M2*U0=U1*D1*U1^T,M2=(U0*U1)*D1*(U0*U1)^T
D1 = diag(D1);
if any(D1 <= 0) == true
   isnegeig = true;
   return;
end
if size(D1,1) ~= k
   error(['D1 size not equal to ' num2str(k)]);
end

U1=U0*U1;

%the largest k orthonormal eigen vecs and eigen values    
W1 = diag(1./sqrt(D1))*U1'; 
% W1*M2*W1'=I,     %dim: k*W
% disp 'finish computing eigs'

punorm3 = sparse(1:D,1:D,doclen3)*dwmat;
% p(w|d)/(l-1)/(l-2)
x31 = full(sum(punorm3, 1))';
% sum(p(w|d)/(l-1)/(l-2), over d)
%         x32 = wdmat*punorm3;
% sum(p(w1,w2|d)/(l-2), over d), symetric
%         eps = mean(doclen3)
%         [i,j,s]=find(x32>eps);
%         x32 = sparse(i,j,s/sum(s),W,W); 
% x32 is normalized
%         fprintf('nnz of x32 %d\n',nnz(x32));
%         save(matfile,'-append','x31','x32');
%save(matfile,'-append','x31');
tdmat=W1*wdmat; % k by D
M3 = ktensor3(doclen3,tdmat)+ktensor3(x31*2,W1);
    % first two parts in the 3rd-order moment in Eq. (8), contribution to
    % M3(tilt)
    % note: ktensor3(v,M)_{i,j,k} = v'*M(i,:).*M(j,:).*M(k,:)
    % trick: ktensor3(v,M)(A,A,A) = ktensor3(v,A'M)
clear x31;
k2 = k*k;
W2 = zeros(k2,W); % ith column of W2 = W1_i X W1_i
for j=1:W
   W2(:,j) = kron(W1(:,j),W1(:,j));
end
cij=tdmat*(punorm3*W2'); % k by k^2
    % trick: 
    % if t32(i,i,j)=x32(i,j)
    % t32(A,A,A)_{i,j,k}=sum_{i1,i2}A(i1,i)A(i1,j)x32(i1,i2)*A(i2,k)
    % = cij(k,i,j) = permute(cij,[2 3 1])(i,j,k)
    % if t32(i,j,i)=x32(i,j)
    % t32(A,A,A)_{i,j,k}=cij(j,i,k)=permute(cij,[2,1,3])(i,j,k)
    % if t32(i,j,j)=x32(i,j)
    % t32(A,A,A)_{i,j,k}=cij(i,j,k)
cij = reshape(cij,k,k,k);
M3 = M3 - cij-permute(cij,[2 1 3])-permute(cij,[2 3 1]);
% the 3rd-order moment E(w1xw2xw3)

WM1 = W1*M1; 
% trick: v x v x v(A,A,A) = (A'v)x(A'v)x(A'v)
% trick: since W'*M2*W=I, M2=x2*(ALPHA0+1)-spM1*spM1'*ALPHA0, 
%       W'*x2*W = (I+W'*spM1*spM1'*W*ALPHA0)/(ALPHA0+1)
WM12 = WM1*WM1';
IplusWM12 = eye(k)+WM12*ALPHA0;
x2M1 = reshape(IplusWM12(:)*WM1',k,k,k);
    % trick: x2 x M1(W,W,W) = x2(W,W)x W'M1
x11 = x2M1 + permute(x2M1,[1 3 2]) + permute(x2M1, [3 1 2]);
% the 2nd part times Alpha0+1 in M3, contribution to T(tilde)
M3 = M3/D*(ALPHA0+2)*(ALPHA0+1)/2 - ALPHA0/2*x11...
        + ALPHA0*ALPHA0*reshape(WM12(:)*WM1',k,k,k);
%disp 'finish computing T(tilde)';

% the following is standard orthogonal tensor decomposition
ALPHA = zeros(1,k);
thetas = zeros(k,k);
n = options.n;
lr = options.lr;

M3 = reshape(M3,k2,k);
for kk=1:k
    disp(['computing topic ' num2str(kk)]);
    maxlambda=-1000;
    thetastar=0;
    for tao = 1:options.N
        theta = W1 * wdmat(:,ceil(rand()*D));%rand(T,1)-0.5;
        theta = theta/norm(theta);
        tmplambda = -inf;
        for t=1:n            
            update = reshape(M3*theta,k,k)*theta;
            if any(~isreal(update)) == 1
              disp 'update is imaginary'
            end
            tmplambda = dot(theta,update);
            theta=theta+lr*update;
            theta = theta./norm(theta);
        end
        if (tmplambda > maxlambda)
            maxlambda = tmplambda;
            thetastar = theta;
        end
    end
    theta = thetastar;
    lambda = maxlambda;
    ALPHA(kk)=(1/lambda)^2;
    thetas(:,kk) = theta;
    M3 = M3 - lambda*reshape(kron(kron(theta,theta),theta),k2,k);
end
wtmat = normalize_cols(U1*(diag(sqrt(D1))*thetas));
% wtmat = bsxfun(@rdivide,wtmat,sum(wtmat, 1));
%disp 'finish tensor decomposition';

function M = normalize_cols(M)  
% this function is the same as in Daniel's implementation
sp = sum(max(M,0),1);
sn = sum(max(-M,0),1);
M = max(M*diag(2*double(sp>sn)-1),0);
M = bsxfun(@rdivide, M, sum(M, 1));
