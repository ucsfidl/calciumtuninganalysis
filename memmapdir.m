function fname=memmapdir(type,s)
%type, predefine memmapchoice type: single file,multiple plane,Multiple
%files,mutiple planesmultiplefiles
%string, specify filenames to be aligned

if ~nargin   %
    memmapchoice=menu('Which type?','Single file','Multiple planes','Multiple files','Multiple planes multiple files')
else
    memmapchoice=type;
end
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
        if exist('s','var')
            d = dir([s '.sbx']);
        else
            d = dir('*.sbx');
        end
        for i=1:numel(d)
            fn = strtok(d(i).name,'.');
            %newname = makememmap(fn);
            newname = makememmap1file(fn);
            pos = strfind(newname,'_memmap.mat');
            fname=memmap(newname(1:pos-1));
            disp(sprintf('finished cell extraction for %s in %d seconds',fname,toc));
        end
    case 2 %Multiple plane
        sbxalign2planedir;
        if exist('s','var')
            d = dir([s '.sbx']);
        else
            d = dir('*.sbx');
        end
        for i=1:numel(d)
            fn = strtok(d(i).name,'.');
            newnames= makememmap2plane(fn);
%             for j=1:2
%                 pos = strfind(newnames{j},'_memmap.mat');
%                 fname= memmap(newnames{j}(1:pos-1));
%                 disp(sprintf('finished cell extraction for %s in %d seconds',fname,toc));
%             end
        end
    case 3 %Multiple files
        sbxaligndir;
        newname=makememmapnfiles;
        pos = strfind(newname,'_memmap.mat');
        fname=memmap(newname(1:pos-1));
        disp(sprintf('finished %s in %d seconds',fname,toc));
    case 4 %Multiple files with multiple planes
        sbxalign2planedir;%%%
        newnames=makememmap2plane2files;
        for i=1:2
            pos = strfind(newnames{i},'_memmap.mat');
            fname=memmap(newnames{i}(1:pos-1));
            disp(sprintf('finished %s in %d seconds',fname,toc));
        end
end

sbxballmotiondir;
sbxeyemotiondir;









