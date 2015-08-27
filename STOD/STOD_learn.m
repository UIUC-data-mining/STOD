function inferred = STOD_learn(dwmat, options)
% dwmat - document word sparse matrix
% options.ALPHA0 - summation of alpha_1, ... , alpha_T
% options.alphaiter - maximal number of shrinking alpha0 in order to find a
%                     valid alpha0
% options.learnalpha - the number of iterations to learn alpha; 0 means no
%                      learning
% options.N - number of outer iterations, for initialization times
% options.n - number of inner iterations, for fixed point iteration times
% options.lr - the learning rate, e.g., options.lr= 1
% options.epsalpha - the convergence threshold for alpha, 
%                    set to a small number,e.g., options.epsalpha = 1e-3
% options.K - the topic number if not to learn the topic number; e.g.,
%             options.K = 5
% 		    - the range of topic number if to learn the topic number,
% 			  e.g. options.K = 2:8
% options.eigen - the algorithm to compute whitening matrix, 
%                 'exact' or 'approx'

% learn k, and compute U0 and D0 by the way
tic0 = tic;
if length(options.K)>1
   % learn topic number
   disp('learning number of topics');
end
if strcmp(options.eigen,'exact')
    disp 'computing whitening matrix using exact algorithm...';
    [issmalldata,iseigsuccess,isasym,k,U0,D0]=decomp0(dwmat, ...
        options);   
else
    disp 'computing whitening matrix using approximate algorithm...';
    [issmalldata,iseigsuccess,isasym,k,U0,D0]=decomp0_lowdim(...
        dwmat,options);
end

disp(['finished computing whitening matrix in ' num2str(toc(tic0)) ' seconds']);
disp '-------------------------';
% if the data is small, or if eigen decomposition never succeeds, 
% should return directly because this dwmat can no longer be partitioned
if issmalldata == true || iseigsuccess == false || isasym == true
   disp 'data is small or cannot find a valid k or E2 asymmetric'
   return;
end

% initialize alpha as user specified, generate new alpha until this alpha 
% makes a positive D1. this process applies regardless of istrainalpha 
ALPHA0 = options.ALPHA0;
tic1 = tic;
disp 'started learning alpha...';
for iternum = 1:options.alphaiter
    [isnegeig, wtmat, alpha1] = decomp(dwmat,k,ALPHA0,U0,D0*(ALPHA0+1),...
        options);     
    if ~isnegeig
        break;
    end
    ALPHA0 = ALPHA0/2;
end

disp(['finished learning alpha in ' num2str(toc(tic1)) ' seconds']);
disp '-------------------------';
if exist('isnegeig','var') && isnegeig
    disp('cannot find a valid alpha, returning...');
    return;
end
for iternum = 1:options.learnalpha
    if abs(sum(alpha1)-ALPHA0)<options.epsalpha
        break;
    end
    lr = options.lr; % learning rate of alpha
    for inneriter = 1:options.alphaiter
        attemptALPHA0 = (ALPHA0+lr*sum(alpha1))/(1+lr);
        [isnegeig, wtmat_tmp, alpha1_tmp]=decomp(dwmat,k,attemptALPHA0, ...
            U0, D0*(attemptALPHA0+1),options); 
        if ~isnegeig
            break;
        end
        lr = lr/2;
    end
    if isnegeig
        % if any small steps cannot make a positive eigen decomposition, 
        % we cannot move forward, so take current alpha
        break;
    end
    ALPHA0 = attemptALPHA0;
    alpha1 = alpha1_tmp;
    wtmat = wtmat_tmp;
    disp(['lr ' num2str(lr) ' ALPHA0 ' num2str(ALPHA0) ' sumalpha1 ' ...
          num2str(sum(alpha1))]); 
end
disp(['tensor decomposition takes ' num2str(toc(tic0)) ' seconds']);

disp '=========================';
clear U0;
inferred.alpha0 = ALPHA0;

ALPHA0 = sum(alpha1);
[inferred.alpha,ind] = sort(alpha1/ALPHA0, 'descend');
inferred.twmat = wtmat(:,ind)';

