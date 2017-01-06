function DSI=calDSI(peak);

dim=ndims(peak);
ncell=size(peak,dim);
ori=size(peak,dim-1);
type=size(peak,dim-2);

if ~(dim==3)
       peak=reshape(peak,[],ori,ncell);
end
    
%peak contrast*ori*ncell
if mod(ori,2)   %if with a blank stim
    ori=ori-1;
end

[preferred,preferred_dir]=max(peak,[],2);
oppo_dir=mod(preferred_dir+ori/2-1,ori)+1;
for i=1:ncell
    A=peak(:,preferred_dir(i),i);
    B=peak(:,oppo_dir(i),i);
    DSI(i)=(A-B)./(A+B)*A;
end



