clear all;
close all;
clc;

im_name='7_9_s.bmp';

% TODO: Update library path
% Add  library paths
%basedir='UGM';
%addpath(basedir);
addpath(genpath("UGM"))

%% Set model parameters
%cluster color
K=4; % Number of color clusters (=number of states of hidden variables)

%Pair-wise parameters
smooth_term=2; % Potts Model

%% Load images
im = imread(im_name);


NumFils = size(im,1);
NumCols = size(im,2);
NumChannels = size(im,3);
NumPixels = NumFils*NumCols;

%%
%Convert to LAB colors space
% TODO: Uncomment if you want to work in the LAB space
%
im = RGB2Lab(im);


%%
%Preparing data for GMM fiting
%
% TODO: define the unary energy term: data_term
% nodePot = P( color at pixel 'x' | Cluster color 'c' )  
im = double(im);
x = reshape(im, [NumPixels, NumChannels]);
gmm_color = gmdistribution.fit(x, K);
mu_color = gmm_color.mu;
data_term=gmm_color.posterior(x);
nodePot=data_term;

%%
%Building 4-grid
%Build UGM Model for 4-connected segmentation
disp('create UGM model');

% Create UGM data
[edgePot, edgeStruct] = CreateGridUGMModel(NumFils, NumCols, K, smooth_term);

%%
if ~isempty(edgePot)

    % color clustering
    [~,c] = max(reshape(data_term,[NumFils*NumCols K]),[],2);
    im_c= reshape(mu_color(c,:),size(im));
    
    % Call different UGM inference algorithms
    display('Loopy Belief Propagation'); tic;
    % [nodeBelLBP,edgeBelLBP,logZLBP] = UGM_Infer_Chain(nodePot,edgePot,edgeStruct);toc; % métode alternatiu
    [nodeBelLBP,edgeBelLBP,logZLBP] = UGM_Infer_LBP(nodePot,edgePot,edgeStruct);toc;
    [~,im_lbp] = max(nodeBelLBP,[],2);
    im_lbp= reshape(mu_color(im_lbp,:),size(im));
    
    % Max-sum
    display('Max-sum'); tic;
    decodeLBP = UGM_Decode_LBP(nodePot,edgePot,edgeStruct);
    im_bp= reshape(mu_color(decodeLBP,:),size(im));
    toc;
    
    
    % TODO: apply other inference algorithms and compare their performance
    %
    % - Graph Cut
    % - Linear Programing Relaxation
    %
    figure
    
    subplot(2,2,1),imshow(Lab2RGB(im));xlabel('Original');
    subplot(2,2,2),imshow(Lab2RGB(im_c),[]);xlabel('Clustering without GM');
    subplot(2,2,3),imshow(Lab2RGB(im_bp),[]);xlabel('Max-Sum');
    subplot(2,2,4),imshow(Lab2RGB(im_lbp),[]);xlabel('Loopy Belief Propagation');
    saveas(figure, 'out.png')
    
else
   
    error('You have to implement the CreateGridUGMModel.m function');

end
