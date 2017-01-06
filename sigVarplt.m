function [h,Gd]=sigVarplt(sigF,window,matrix,Cor); 
%sigF seg*rep*Var*ncell
%maxtrix rep*Var
global info;
name=info.var(strcmp(info.var,'None')==0);
seg=size(sigF,1);
rep=size(sigF,2);
ncell=size(sigF,4);
Gd=[];

if numel(info.var)>1
Var(1)=info.steps(1);
Var(2)=info.steps(2);
end
if prod(Var)~=size(sigF,3)
   disp(' stimulus varation and data size donnot match!!!');
   return;
end


if ~exist('Cor')
    Cor=1:ncell;
end

if isempty(matrix)
    matrix=ones(rep,prod(Var))>0;
end

page=ceil(ncell/10);
tic
for j=1:page
    num=min(10,ncell-10*(j-1)); %plot 10 or mod(ncell,10) cells per page
    h(j)=figure('Name',['reorganized_ page#' num2str(j)], 'position',[ 200 100 1600 100*num]);
    [xpos,ypos,xwidth,yheight]=figurepara(Var(1),num,1,seg);
    for i=1:num    %%%%%cell#=i+10*(j-1)
            nth=i+10*(j-1);
            ymax=max(reshape(sigF(:,:,:,nth),[],1));
            subplot('position',[xpos(1) ypos(i) xwidth yheight/Var(2)]);
            text(0,0,['cell#' num2str(Cor(nth))],'HorizontalAlignment','right','VerticalAlignment','bottom');
            
        for k=1:Var(1)   %%%ori

            for m=1:Var(2)   %%%%%var#=m+Var(2)*(k-1)
            subplot('position',[xpos(k) ypos(i)+yheight/Var(2)*(m-1) xwidth yheight/Var(2)]);hold on;
            kth=m+(k-1)*Var(2);
            tempR=sigF(:,matrix(:,kth),kth,nth);        %tempR  seg,some rep,1,1;
            tempS=sigF(:,~matrix(:,kth),kth,nth);
            plot(tempR,'k','linewidth',.5);
            plot(tempS,'g','linewidth',.5);
            plot(1:seg,mean(tempR,2),'r','linewidth',2);
            plot(1:seg,mean(tempS,2),'b','linewidth',2);
            plot([window(1) window(1)] ,[0 ymax] ,[window(end) window(end)], [0 ymax],'--','linewidth',1);% stimulus on time
            axis([1 seg 0 ymax]);
            a=gca;
            set(a,'xticklabel',[]);axis off;
            
            end

        end
        text(seg*1.1,ymax,['Ori' num2str(k)],'VerticalAlignment','bottom');

        axis on;drawnow;
        answer = inputdlg('Enter gd cells','Pick gd cells',[1 80],{num2str(Cor(nth))});
        Gd=[Gd str2num(answer{1})];
    end

end
disp(toc)
