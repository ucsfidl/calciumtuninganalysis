function V=sbxVariancemap2plane(fn)
n=2;

sbxread(fn,1,1);        % read one frame to read the header of the image sequence
global info;            % this contains the information about the structure of the image
V=zeros([info.sz n]);
fac=sqrt((info.max_idx+1)/n);
Y=zeros(info.sz);

tic
for i=0:info.max_idx
    z =sbxread(fn,i,1);
    z= squeeze(z);
    Y = circshift(z,info.aligned.T(i+1,:)); % align the image
    nplane=mod(i,n)+1;
    temp=((double(Y-info.aligned.m(:,:,nplane)))/fac).^2; 
    V(:,:,nplane)=V(:,:,nplane)+temp;
end
save([fn '.align'],'V','-append');