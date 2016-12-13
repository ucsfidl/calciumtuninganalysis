function hsig=sigplt(sig,stim,Cor);

ncell=size(sig,2);
nframe=size(sig,1);
page=ceil(ncell/10);

if ~exist('stim')
    stim=[];
end
if ~exist('Cor')
    Cor=1:ncell;
end



for j=1:ceil(ncell/10)
    num=min(10,ncell-10*(j-1)); %plot 10 or mod(ncell,10) cells per page            
    hsig(j)=figure('Name',['Ca sig w stim_ page#' num2str(j)], 'position',[200 100 2400 100*num]);
    [xpos,ypos,xwidth,yheight]=figurepara(1,num);

    for k=1:num
        nth=k+10*(j-1);
        subplot('position',[xpos ypos(k) xwidth yheight]);hold on;
        plot(1:nframe,sig(:,nth),'linewidth',1);        
        plot(1:nframe,stim','linewidth',.5);
        text(nframe*1.05,0,['Cell#' num2str(Cor(nth))],'HorizontalAlignment','right');
        %xlim([1 nframe]);
        ymax=max(sig(:,nth));
        subplot('position',[xpos+xwidth*1.05 ypos(k) .1 yheight]);hold on;
        try
            histogram(sig(:,nth),[0:ymax/50:ymax],'Normalization','probability')
            
        catch
            title('noise');
        end
    end
end