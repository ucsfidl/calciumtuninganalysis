function fname=memmapdir(s)
clc;close all; clearvars -except s

memmapchoice=menu('Which type?','Single file','Multiple files','Multiple planes multiple files')
%sbxballmotiondir;
%sbxeyemotiondir;


% sbxaligndir;
% global info_loaded info
% if(~isempty(info_loaded))   % try closing previous...
%     try
%         fclose(info.fid);
%     catch
%     end
% end

switch memmapchoice
    case 1  %single file
        sbxaligndir;
        makememmapdir;
        if nargin % cell with filenames to be aligned
            d = dir([s '.sbx']);
        else
            d = dir('*.sbx');
        end
        for i=1:numel(d)
            fn = strtok(d(i).name,'.');
%            sbxVariancemap(fn);
            fname{i}=memmap(fn);
            disp(sprintf('finished %s in %d seconds',fname{i},toc));
        end
        
        newname=makememmap2files;
        pos = strfind(newname,'_memmap.mat');
        fname=memmap(newname(1:pos-1));
        disp(sprintf('finished %s in %d seconds',fname,toc));
    case 2 %Multiple files
        sbxaligndir;
        newname=makememmap2files;
        pos = strfind(newname,'_memmap.mat');
        fname=memmap(newname(1:pos-1));
        disp(sprintf('finished %s in %d seconds',fname,toc));
    case 3 %Multiple planes
        sbxalign2planedir;%%%
        newnames=makememmap2plane2files;
        for i=1:2
            pos = strfind(newnames{i},'_memmap.mat');
            fname{i}=memmap(newnames{i}(1:pos-1));
            disp(sprintf('finished %s in %d seconds',fname{i},toc));
        end
end
end









