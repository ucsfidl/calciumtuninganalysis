function newnames=makememmap2plane(fn);
%% define the name
path=pwd;
zfs_path=strrep(path,'c:','\\mps-zfs\data\jsun');

newnames{1}=fullfile(zfs_path,[fn '_1_memmap.mat']);
newnames{2}=fullfile(zfs_path,[fn '_2_memmap.mat']);

if exist(newnames{1},'file')
    display('memmap file already exist!'); return;
elseif ~exist(zfs_path,'dir')
    mkdir(zfs_path);
end

data1 = matfile(newnames{1},'Writable',true);
data2 = matfile(newnames{2},'Writable',true);

disp(sprintf('Saving memmapfiles for 2 planes:%s,%s',newnames{1},newnames{2}));

%% save each piece
sbxread(fn,1,1);% read one frame to read the header of the image sequence
global info;

nY = Inf;
T=(info.max_idx+1)/2;
m1=info.aligned.m(:,:,1);
m2=info.aligned.m(:,:,2);
V1=zeros(info.sz);
V2=zeros(info.sz);
fac=sqrt(T-1);
data1.Y=zeros([info.sz T],'uint16');
data2.Y=zeros([info.sz T],'uint16');

for j=0:info.max_idx
    z =sbxread(fn,j,1);
    z= squeeze(z);
    if mod(j,2)==0
        img = circshift(z,info.aligned.T(j+1,:)); % align the image
        data1.Y(:,:,j/2+1) = img;
        V1=V1+((double(img-m1))/fac).^2;
    else
        img = circshift(z,info.aligned.T(j+1,:)); % align the image
        data2.Y(:,:,(j+1)/2) = img;
        V2=V2+((double(img-m2))/fac).^2;
    end
    
    if mod(j,500)==0
        fprintf('Saved Frame %d/%d for %.2f seconds\n ',j,info.max_idx,toc);
    end
end
%V=cat(3,V1,V2); save([fn '.align'],'V','-append');


data1.V=V1;
data2.V=V2;
% assert(T==(info.max_idx+1)/2,'the total image size is not matched!');
data1.sizY =[info.sz T];
data2.sizY =[info.sz T];

data1.Yr = reshape(data1.Y,prod(info.sz),T);
data2.Yr = reshape(data2.Y,prod(info.sz),T);

data1.nY = min(reshape(data1.Yr,prod(info.sz)*T,1));
data2.nY = min(reshape(data2.Yr,prod(info.sz)*T,1));

data1.m = m1;
data2.m = m2;

if mod(info.config.magnification*2,5) == 0
    magnification = info.config.magnification*2/5;
else
    magnification = info.config.magnification;
end

data1.magnification = magnification;
data2.magnification = magnification;
