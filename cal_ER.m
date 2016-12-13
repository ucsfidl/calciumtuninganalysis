function [peak,idx,deviation]=cal_ER(sigF,sigwin)
%calculate responses of each stimulus after averaging

%sigF:seg,rep,Var,ncell
% peak 1*Var*ncell
%error:1,Var*ncell
%idx

%% filter trace : smoothing and subtract baseline
%sigF=interp1(1:seg,sigF,1:1/mag:seg);
seg=size(sigF,1);
rep=size(sigF,2);
Var=size(sigF,3);
ncell=size(sigF,4);

winSize = 5;
b = (1/winSize)*ones(1,winSize);
a = 1;
sigF = filter(b,a,sigF,[],1);
baseline=mean(sigF(sigwin(1)-2:sigwin(1)+2,:,:,:));
sigF=sigF-repmat( baseline,[seg 1 1 1 ]);

%% exclude non-significant response 

%averaging over #of repetitions
sigA= nanmean(sigF,2);  %sigA:seg,1,Var,ncell
%define response: prestim mean+2*std
threshold=mean(sigA(1:sigwin(1),:,:,:))+2*std(sigA(1:sigwin(1),:,:,:));% threshold size: 1,1,Var,ncell
threshold=max(threshold,0);
% set the  trace<threshold as NaN;
sigA(repmat(nanmean(sigA(sigwin,:,:,:))<threshold,[seg 1 1 1 ]))=0;  % seg,1,Var,ncell

%%%%OPTION A. maxium
% [peak,idx]=max(sigA(window,:,:,:),[],1);  % find peak from local maxium only 1*1*Var*ncell
% idx=idx+window(1)-1;
% peak=reshape(peak,1,Var,ncell);% peak 1*Var*ncell
% idx=reshape(idx,1,Var,ncell);  %idx 1*Var*ncell

%%%%%%%%OPTION B. mean of the response window
peak=mean(sigA(sigwin,:,:,:),1);
peak=reshape(peak,1,Var,ncell); %%just changed to match the OPTION.A dimension. 10/7
idx=[];
%figure;hold on;plot(reshape(sigA,[],Var*ncell));plot(idx(:),peak(:),'*');
%% error calculation need to be updated, because of the zeros

%error_all=std(sigF(repmat(~badmatrix,[seg 1 1 1]),2))  %error:seg,1,Var,ncell
%temp=reshape(error_all,seg,Var,ncell);%temp:seg,Var,ncell
%error=temp(:,idx);
% if rep =1 , or using mean response would cause error
try
    temp=permute(sigF,[2 1 3 4 ]);  %temp:rep,seg,Var,ncell
    pos=sub2ind([seg,Var,ncell],idx,reshape(1:prod(size(idx)),1,Var,ncell));
    deviation=nanstd(temp(:,pos)); % deviation: 1,1,Var,ncell
    deviation=reshape(deviation,1,Var,ncell);
catch
%     deviation=zeros(size(peak)); % rep=1
    deviation=std(sigA,0,1); 
    deviation=reshape(deviation,1,Var,ncell);

end
 