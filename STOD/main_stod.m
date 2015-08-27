% the main program for running STOD algorithm
% input data: folder/dataname.dp 
% function main()

path('../DataProcess/readdata/', path);

folder = '../Data/'; % the folder for input data
dataname = 'test'; % the name for the dataset

tic
disp '=========================';
disp 'loading data';
[dw, dwmat] = ReadEdge([folder dataname '.corpus']);
name = ReadName([folder dataname '.dict']);

disp(['finished loading data, loading takes ' num2str(toc) ' seconds']);
disp '=========================';

vocabulary=name{1};
[docnum, wordnum] = size(dwmat);

% parameters that are usually not tuned
options.folder = folder;    % output file location
options.dataname = dataname;    % output file prefix
options.N = 30;         % outer iteration times
options.n = 30;         % inner iteration times
options.lr = 1.0;       % learning rate


% options.K = 2:10;     % the range of topic numbers. 
options.K = 5;        % to fix it, use a single number 

% options.eigen = 'approx';   % fast algorithm to compute whitening matrix
options.eigen = 'exact';  % exact eigen decomposition

% parameters that are tuned

% parameter to tune K selection, which is the proportion of 
% sum(first K topics eigenvalues) / sum(all eigenvalues), range from 0 to 1
options.proportion = 0.95; 

options.ALPHA0 = 1.0; % initial ALPHA0
% options.ALPHA0 = 2.0235;
% options.learnalpha = 20;    % the maximal # trials for searching alpha
options.learnalpha = 0;     % fix ALPHA0 by setting it to 0
options.epsalpha = 1e-3;    % the convergence threshold for alpha

% the maximal # trials to shrink alpha by half 
% in order to have enough pos eigs
options.alphaiter = 10; 

disp 'learning STOD';
disp '-------------------------';
inferred = STOD_learn(dwmat, options);
disp 'finished constructing, saving the mat files...'

matfile = strcat(options.folder, dataname, '.stod', '.mat');
save(matfile, 'inferred', '-v7.3');

inferred.alpha = inferred.alpha * inferred.alpha0;

k = length(inferred.alpha);
twmat = inferred.twmat;
topics = cell(1, k);
for z=1:k
	phi = twmat(z,:);
    [~, ind] = sort(phi, 'descend');
    topicz = vocabulary(ind(1, 1:10));
    topics{1, z} = topicz;
end

matfile = strcat(options.folder, dataname, '.topic', '.mat');
save(matfile, 'topics', '-v7.3');