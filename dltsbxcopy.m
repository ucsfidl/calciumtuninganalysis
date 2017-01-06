function dltsbxcopy
zfs_path=pwd;
cpath=strrep(zfs_path,'\\mps-zfs\data\jsun','c:');
allfiles=ls;

for i=3:size(allfiles,1)
    if strfind(allfiles(i,:),'.sbx') | strfind(allfiles(i,:),'ball') | strfind(allfiles(i,:),'eye') 
        if exist(fullfile(cpath,allfiles(i,:)),'file')
            display(sprintf('extra copy deleted: %s',allfiles(i,:)))
            delete(allfiles(i,:));
        else
            display(sprintf('single copy moved back to C: %s',allfiles(i,:)))
            movefile(allfiles(i,:),cpath);
        end
    end                    
    if isdir(allfiles(i,:))
        cd(allfiles(i,:));dltsbxcopy; cd ..
    end
end