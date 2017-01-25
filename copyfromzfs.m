function copyfromzfs

zfs_path=pwd;
cpath=strrep(zfs_path,'\\mps-zfs\data\jsun','C:');

if ~exist(cpath)
    mkdir(cpath)
end

try
copyfile('*.signals',cpath);
catch
end

try
copyfile('*.png*',cpath);
catch
end

allfiles=ls;
for i=3:size(allfiles,1)
    if isdir(allfiles(i,:))
        cd(allfiles(i,:));copyfromzfs; cd ..
    end
end
    