function makememmapdir(varargin)

if(nargin>=1) % cell with filenames to be aligned
    for(i=1:length(varargin{1}))
        d(i).name = varargin{1}{i};
    end
else
    d = dir('*.sbx');
end

% Align all *.sbx files in the list

for(i=1:numel(d))
    try
        fn = strtok(d(i).name,'.')
        if(exist([fn '_memmap.mat'])==0)
            tic
            makememmap(fn);
            display(sprintf('saved memmap file in %d min',round(toc/60)));
        else
            sprintf('%s_memmap.mat is already made',fn)
        end
    catch
        sprintf('Could not make memmap file %s',fn)
    end
end