function fname=memmap(nam,patchf);
%% load file
%path_to_package = 'C:\Users\kraken\Documents\MATLAB\ca_source_extraction\utilities';   % path to the folder that contains the package
%addpath(genpath(path_to_package));
%close all;

magnification=1;
if ~exist([nam '_memmap.mat'],'file')
    makememmap(nam);
end
data = matfile([nam '_memmap.mat'],'Writable',true);

%% Set parameters
%tic
% try
%     sbxread(nam,1,1);
% catch
%     filenam=dir('*.sbx');
%     display(filenam(1).name);
%     sbxread(filenam(1).name(1:end-4),1,1);
% end

try
    sizY=data.sizY;
catch
    sizY = size(data,'Y');                  % size of data matrix
end
% 
% if ~exist('patchf')
%     %     mag=ceil(sqrt(sizY(3)/16000)); % if too many frames, chop the patch sizes to smaller ones
%     mag=1;
%     patchf=[ceil(sizY(1)*mag/80) ceil(sizY(2)*mag/80)];
% end
% 
% display(patchf);
% patch_size = [ceil(sizY(1)/patchf(1)*1.1),ceil(sizY(2)/patchf(2)*1.1)];                   % size of each patch along each dimension (optional, default: [32,32])
% overlap = [ceil(sizY(1)/patchf(1)*.05),ceil(sizY(2)/patchf(2)*.05)];                        % amount of overlap in each dimension (optional, default: [4,4])

patch_size = [64,64];                   % size of each patch along each dimension (optional, default: [32,32])
overlap = [8,8];                        % amount of overlap in each dimension (optional, default: [4,4])

patches = construct_patches(sizY(1:end-1),patch_size,overlap);
%K = ceil(sizY(1)*sizY(2)/10000/prod(patchf))    % number of components to be found
%K=ceil(40/prod(patchf));
K=2;
tau = 2*magnification;                                          % std of gaussian kernel (half width/height of neuron)
p = 0;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.6;                                  % merging threshold
tsub=32;
options = CNMFSetParms(...
    'd1',sizY(1),'d2',sizY(2),...
    'search_method','ellipse','dist',3,...      % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps
    'ssub',magnification,...
    'tsub',tsub,...
    'fudge_factor',0.98,...                     % bias correction for AR coefficients
    'merge_thr',merge_thr,...                    % merging threshold
    'gSig',tau...
    );
%display(sprintf('CNMFparms done %d min',round(toc/60)));

%% Run on patches
tic;

[A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(data,K,patches,tau,p,options);

display(sprintf('run patches in %d min',round(toc/60)));

%% DONnotorder , just plot
contour_threshold = 0.95;                 % amount of energy used for each component to construct contour plot
K_m = size(C,1);
%[A_or,C_or,S_or,P] = order_ROIs(A,C,S,P); % order components
try
    [sig,S_df] = extract_DF_F(data,[A_or,b],[C_or;f],K_m+1); % extract DF/F values (sig) and spiking info
catch
    [sig,S_df] = extract_DF_F(data,[A,b],[C;f],K_m+1); % extract DF/F values (sig) and spiking info
end

% try
%     V=data.V;
% catch
%     m=info.aligned.m;
%     V=zeros(sizY(1),sizY(2));
%     fac=sqrt(sizY(end));
%     for i=1:sizY(end)
%         V=V+((double(data.Y(:,:,i)-m))/fac).^2;
%     end
% end
V=reshape(P.sn,sizY(1),sizY(2));

try
    [Coor,json_file] = plot_contours(A_or,V,contour_threshold,1);
catch
    [Coor,json_file] = plot_contours(A,V,contour_threshold,1);
end
%% save signals into different files
fname=sprintf('%s_%dcell_tau%d_tsub%d',nam,K_m,tau,tsub);

sig=sig';
save([fname '.signals'],'sig','S_df');
savejson('jmesh',json_file,[fname '.jmesh']);        % optional save json file with component coordinates (requires matlab json library)
savefig(gcf,[fname '.fig']);

try
    T=0;
    for n=1:numel(data.eachsize)
        sig_chunk=sig(T+1:T+data.eachsize(1,n),:);
        save([fname '_' num2str(n) '.signals'],'sig_chunk','S_df');
        T=T+data.eachsize(1,n);
    end
catch
end