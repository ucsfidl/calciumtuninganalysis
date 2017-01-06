function sbxalign2planedir(varargin)

if (nargin>=1) % cell with filenames to be aligned
    for(i=1:length(varargin{1}))
        d(i).name = varargin{1}{i};
    end
else
    d = dir('*.sbx');
end

% Align all *.sbx files in the list

for(i=1:length(d))
    try
        fn = strtok(d(i).name,'.');
        if ~ exist([fn '.align'])
            sbxread(fn,1,1);            % read one frame to read the header of the image sequence
            global info;                % this contains the information about the structure of the image
            tic
            T=zeros([info.max_idx+1 2],'uint16');
            [m,T(1:2:info.max_idx+1,:)] = sbxalignx(fn,0:2:info.max_idx);   %
            [m(:,:,2),T(2:2:info.max_idx+1,:)] = sbxalignx(fn,1:2:info.max_idx);
            save([fn '.align'],'m','T');
            clear m T;
            display(sprintf('Done %s: Aligned %d images in %d min',fn,info.max_idx,round(toc/60)));
        else
            sprintf('File %s is already aligned',fn)
        end
    catch
        sprintf('Could not align %s',fn)
    end
end
