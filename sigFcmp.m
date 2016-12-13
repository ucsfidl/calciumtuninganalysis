function [sigR,peakR,errorR,sigS,peakS,errorS]= sigFcmp(sigF,sigwin,matrix); 
%sigF  seg,rep,Var,ncell;  
%sigR,sigS seg*Var*ncell; 
%maxtrix rep*Var
Var=size(sigF,3);
  
for k=1:Var
    tempR=sigF(:,matrix(:,k),k,:);        %tempR  seg,some rep,1,ncell
    sigR(:,k,:)=squeeze(mean(tempR,2));
%     [peakR(1,k,:),idxR(1,k,:),errorR(1,k,:)]=cal_ER(tempR,sigwin);

[peakR(1,k,:),~,errorR(1,k,:)]=cal_ER(tempR,sigwin);
    tempS=sigF(:,~matrix(:,k),k,:);
%     [peakS(1,k,:),idxS(1,k,:),errorS(1,k,:)]=cal_ER(tempS,sigwin);
    [peakS(1,k,:),~,errorS(1,k,:)]=cal_ER(tempS,sigwin);
    sigS(:,k,:)=squeeze(mean(tempS,2));
end