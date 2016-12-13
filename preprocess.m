function movement=preprocess(fname)
global info;

load([fname '.align'],'m','-mat');
sbxread(fname,1,1);
if info.frame(1)==0
    info.frame=info.frame(2:end);
    info.line=info.line(2:end);
    save([fname '.mat'],'info','-append');
end

if ~isfield(info,'stimtype')   %read psychostim stim info to match TTL event recorded in
    file=dir([fname '.sbx']);
    recordedtime=file.date;
    [day,time]=strtok(recordedtime)
    [f2,p2]=uigetfile(['..\..\stim\' day],['find stim info just before' time]); %'
    %[f2,p2]=uigetfile(['\\mps-pc38\Documents\Matlab\' day],['find stim info just before' time]);%'
    load(fullfile(p2,f2), '-mat');
    sttype = intanSyncData(:,1);% the variable IntanSyncData has all of the useful information
    info.stimtype=sttype(1: find(diff(sttype),1):end);%find(diff(sttype),1) can define the length of ON and OFF, so, info.stimtype reads the first frame of each stimulus
    info.var=Var1Str{Var1Val};
    info.steps=nSteps1;
    if Var2Val>1
        info.var={info.var,Var2Str{Var2Val}};
        info.steps=[info.steps,nSteps2];
    else
        info.var={info.var,'none'};
        info.steps=[info.steps,1];
        save([fname '.mat'],'info','-append');

    end
    
end

short=abs(2*numel(info.stimtype)-numel(info.frame));
if short
    warning=sprintf(['stimtype#%d and recorded#%d do not match,check if necessary',2*numel(info.stimtype),numel(info.frame)]);
    disp(warning);
else
    save([fname '.mat'],'info','-append');disp('stim and recording matched and saved!');    
end

if exist([fname '_ball.mat'])
    movement=1;
    if ~exist([fname '.ball'])
        
        sbxballspeed([fname '_ball.mat']) ;
    end
    display('with ball info, will analyze running later');
else
    display('no ball info');
    movement=0;
end

if exist([fname '_eye.mat'])
    movement=movement+10;
    if ~exist([fname '.eye'])
        sbxeyemotion([fname '_eye.mat']) ;
    end
    display('with movement info, will analyze running later');
else
    display('no eye info');
end