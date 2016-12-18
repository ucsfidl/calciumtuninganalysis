function fullname=makememmap2files;
%% define the name
files = dir('*.sbx');
path=pwd;
if ~strfind(path,'\\MPS-ZFS')
    newpath=fullfile('\\MPS-ZFS\data\jsun',path(3:end));
else
    newpath=path;
end
newname=[files(1).name(1:end-5) 'x_memmap.mat'];
fullname=fullfile(newpath,newname);

if exist(fullname,'file')
    disp([fullname 'already existing']);
    return;
elseif ~exist(newpath,'dir')
    mkdir(newpath);
end
data = matfile(fullname,'Writable',true);
disp(sprintf(' try saving memmap file: %s',fullname));

%% save each piece
T=0;

for i=1:numel(files)
    fn=files(i).name(1:end-4);
    sbxread(fn,1,1);% read one frame to read the header of the image sequence
    global info;            % this contains the information about the structure of the image
    m(:,:,i)=info.aligned.m;
    eachsize(i)=info.max_idx;
   
    V=zeros(info.sz);
    fac=sqrt(info.max_idx);
    tic
    if i==1
        data.m=zeros(info.sz);      
        data.Y=zeros([info.sz info.max_idx],'uint16');
        data.V=zeros(info.sz);
        u=0;v=0;
    else
        [u v] = fftalign(m(:,:,i),m(:,:,1));
        figure('Name',['aligned file#' num2str(i)]);
        imshowpair(circshift(squeeze(m(:,:,i)),[u v]),squeeze(m(:,:,1)))
    end
    for j=0:info.max_idx-1
        z =sbxread(fn,j,1);
        z= squeeze(z);
        data.Y(:,:,T+j+1) = circshift(z,info.aligned.T(j+1,:)+[u v]); % align the image
        V=V+((double(data.Y(:,:,i+1)-m(:,:,i)))/fac).^2; 
        if mod(j,500)==0
            fprintf('File %d Frame %d/%d for %.2f seconds\n ',i,j,info.max_idx,toc);
        end
    end
    save([fn '.align'],'V','-append');
    ratio=T/(T+eachsize(i));
    data.V=data.V*ratio+circshift(V,[u v])*(1-ratio);
    data.m=data.m*ratio+circshift(double(m(:,:,i)),[u v])*(1-ratio)
    T=T+eachsize(i);
end
%assert(T==info.max_idx,'the total image size is not matched!')
data.eachsize=eachsize;
data.sizY =[info.sz T];
data.Yr = reshape(data.Y,prod(info.sz),T);
data.nY = min(reshape(data.Yr,prod(info.sz)*T,1));

