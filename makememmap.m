function makememmap(fn)
global info_loaded info
if(~isempty(info_loaded))   % try closing previous...
    try
        fclose(info.fid);
    catch
    end
end
sbxread(fn,1,1);        % read one frame to read the header of the image sequence
global info;            % this contains the information about the structure of the image
m=info.aligned.m;
Y=zeros([info.sz info.max_idx],'uint16');
V=zeros(info.sz);
fac=sqrt(info.max_idx);
tic
for i=0:info.max_idx-1
    z =sbxread(fn,i,1);
    z= squeeze(z);
    Y(:,:,i+1) = circshift(z,info.aligned.T(i+1,:)); % align the image
    temp=((double(Y(:,:,i+1)-m))/fac).^2;
    V=V+temp;
    if mod(i,500)==0
        fprintf('Frame %d/%d\n',i,info.max_idx);toc;
    end
end

sizY = size(Y);
Yr = reshape(Y,prod(sizY(1:end-1)),[]);
nY = min(Yr(:));

disp(' try saving memmap file');

save([fn '_memmap.mat'],'sizY','Yr','Y','nY','V','-v7.3');
save([fn '.align'],'V','-append');
toc

