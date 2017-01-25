
function h=sigFplt(sigF,matrix,window,Cor)
% sigF seg*rep*Var*ncell
% maxtrix rep*Var


%global info;
%name=info.var(strcmp(info.var,'None')==0);
seg=size(sigF,1);
rep=size(sigF,2);
Var=size(sigF,3);
ncell=size(sigF,4);
%Gd=[];
if ~exist('Cor')
    Cor=1:ncell;
end

if isempty(matrix)
    matrix=logical(ones(rep,Var));
end

page=ceil(ncell/10);
tic
for j=1:page
    num=min(10,ncell-10*(j-1)); %plot 10 or mod(ncell,10) cells per page
    h(j)=figure('Name',['reorganized_ page#' num2str(j)], 'position',[ 200 100 1600 100*num]);
    [xpos,ypos,xwidth,yheight]=figurepara(Var,num,1,seg);
    for i=1:num    %%%%%cell#=i+10*(j-1)
        nth=i+10*(j-1);
        ymax=max(reshape(sigF(:,:,:,nth),[],1));  %let's just try with max for now
        subplot('position',[xpos(1) ypos(i) xwidth yheight]);
        for k=1:Var
            try
            subplot('position',[xpos(k) ypos(i) xwidth yheight]);hold on;
            tempR=sigF(:,matrix(:,k),k,nth);        %tempR  seg,some rep,1,1;
            tempS=sigF(:,~matrix(:,k),k,nth);
            plot(tempR,'k','linewidth',.5);
            plot(tempS,'g','linewidth',.5);
            plot(1:seg,mean(tempR,2),'r','linewidth',2);
            plot(1:seg,mean(tempS,2),'b','linewidth',2);
            plot([window(1) window(1)] ,[0 ymax] ,[window(end) window(end)], [0 ymax],'--','linewidth',1);% stimulus on time
            axis([1 seg 0 ymax]);
            a=gca;
            set(a,'xticklabel',[]);
            text(0,0,['Var' num2str(k)],'VerticalAlignment','top');
            end
        end
        text(seg*1.1,ymax/2,['cell#' num2str(Cor(nth))],'HorizontalAlignment','right');
        
    end
     drawnow;
%     answer = inputdlg('Enter gd cells','Pick gd cells',[1 80],{num2str(Cor(nth))});
%     Gd=[Gd str2num(answer{1})];   
end
disp(toc)
