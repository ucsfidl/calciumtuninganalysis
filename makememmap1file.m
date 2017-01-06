function newname=makememmap1file(fn);
%% define the name and save to MPS-ZFS directly
path=pwd;
zfs_path=strrep(path,'C:','\\mps-zfs\data\jsun');
newname=[fn '_memmap.mat'];
newname=fullfile(zfs_path,newname);

if exist(newname,'file')
    disp([newname ' memmap file already existing']);return;
elseif ~exist(zfs_path,'dir')
    mkdir(zfs_path);
end
data = matfile(newname,'Writable',true);
display(sprintf(' try saving memmapfile: %s',newname));
starttime=datetime;
%% save 1 piece
global info;            % this contains the information about the structure of the image
sbxread(fn,1,1);% read one frame to read the header of the image sequence

T0=info.max_idx+1;
m=info.aligned.m;
V=zeros(info.sz);
fac=sqrt(info.max_idx);

data.m=zeros(info.sz);
data.Y=zeros([info.sz T0],'uint16');
data.V=zeros(info.sz);

tic
for j=0:info.max_idx
    z =sbxread(fn,j,1);    z= squeeze(z);
    img = circshift(z,info.aligned.T(j+1,:)); % align the image
    data.Y(:,:,j+1) = img;
    V=V+((double(img-m))/fac).^2;
    if mod(j,1000)==0
        fprintf('File %d Frame %d/%d for %.2f seconds\n ',i,j,info.max_idx,toc);
    end
end
save([fn '.align'],'V','-append');
data.V=V;
data.m=m;
data.sizY =[info.sz T0];
data.Yr = reshape(data.Y,prod(info.sz),T0);
data.nY = min(reshape(data.Yr,prod(info.sz)*T0,1));
if mod(info.config.magnification*2,5) == 0
    data.magnification = info.config.magnification*2/5;
else
    data.magnification = info.config.magnification;
end
display(sprintf('Finished making memmap %s',newname));
display(starttime)
display(datetime)
