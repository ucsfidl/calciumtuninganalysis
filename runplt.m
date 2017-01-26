function [matrix,h]=runplt(run); % run:rep,Var,stimON
global info;
rep=size(run,1);
Var=size(run,2);

%% define threshold for running/ still status
h(1)=figure('Name','running speed','position',[ 200 200 600 600]);
%run_diff=max(diff(run,1,3),[],3); % run_diff:rep,Var
run_max=squeeze(max(run,[],3));
subplot(3,1,1);
%plot(1:rep,run_max','*');
histogram(run_max',0:.25:10);
title('max run')
%legend('show')
xlabel('speed')
drawnow;

subplot(3,1,2)
run_avg=squeeze(mean(run,3));% run_avg:rep,Var
%plot(1:rep,run_avg','*');
histogram(run_avg,0:.25:5);
title('average speed')
xlabel('speed')
drawnow;

cutoff1=inputdlg('define states','',1,{'1'});
acc_thr=str2num(cutoff1{1});
% subplot(1,3,1);hold on;
% plot([ 1 rep],[acc_thr acc_thr],'--');
cutoff2=inputdlg('define states','',1,{'.5'});
run_thr=str2num(cutoff2{1});
% subplot(1,3,2);hold on;
% plot([ 1 rep],[run_thr run_thr],'--');

subplot(3,2,5)
ACC=run_max>acc_thr;
imagesc(ACC);

subplot(3,2,6)
RUN=run_avg>run_thr;
imagesc(RUN);
matrix=RUN&ACC;  %%%%matrix=rep*stim    %%%get two methonds, pick one to determine different states
%% to show the running/still status for each recorded segment
h(2)=figure('Name','running segment','position',[ 200 200 1200 800]);

if length(info.steps)==1
    info.steps(2)=1;
end
[xpos,ypos,xwidth,yheight]=figurepara(info.steps(1),info.steps(2));
ymax=max(run(:));
for i=1:info.steps(1)
    for j=1:info.steps(2)
        nth=j+(i-1)*info.steps(2);
        subplot('Position',[xpos(i) ypos(j) xwidth yheight]);hold on      
        plot(squeeze(run(matrix(:,nth),nth,:))','r');
        plot(squeeze(run(~matrix(:,nth),nth,:))','b');      
        ylim([0 ymax]);
        title(['ori' num2str(i) 'ctr' num2str(j)])
    end
end

if any(all(matrix)==1) || any(all(~matrix)==1)    
    choice = questdlg([' status are not covered for every type of stimulus' num2str(find(all(matrix)==1)) num2str(find(all(~matrix)==1))],'','no running','ignore','rerun','rerun' );    
    switch choice
        case 'rerun'
            [matrix,h]=runplt(run);
        case 'no running'
            matrix=[];
        case 'ignore'
    end
end


