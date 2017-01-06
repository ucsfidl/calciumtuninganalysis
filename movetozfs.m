function movetozfs

cpath=pwd;
zfs_path=strrep(cpath,'c:','\\mps-zfs\data\jsun');
if ~exist(zfs_path)
    mkdir(zfs_path)
end
try
movefile('*tau*',zfs_path);
catch
end
try
movefile('*memmap*',zfs_path);
catch
end
allfiles=ls;
for i=3:size(allfiles,1)
    if isdir(allfiles(i,:))
        cd(allfiles(i,:));movetozfs; cd ..
    end
end
    