function sbxballmotiondir(s)

if(nargin>=1)
    d=dir([s '_ball.mat']);
else
    d = dir('*_ball.mat');
end

for i=1:length(d)
    try
        fn = d(i).name;
        vars = whos('-file',fn);
        if (exist([fn(1:strfind(fn,'_ball')-1) '.ball'])==0 && ~ismember('ball', {vars.name}))
            display(sprintf('extract ball motion for %s',fn))
            sbxballmotion(fn);
        else
            sprintf('Ball motion %s is already extracted',fn)
        end
    catch
        display(sprintf('Could not do ballmotion for %s',fn))
    end
end
