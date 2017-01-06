function OSI=calOSI(peak)
%peak 1*type(contrast/running)s*ori*ncell
dim=ndims(peak);
ncell=size(peak,dim);
ori=size(peak,dim-1);
type=size(peak,dim-2);

if ~(dim==3)
       peak=reshape(peak,[],ori,ncell);
end
if ~ mod(ori,2)==0
    ori=ori-1;
end
    
    
if  ~mod(ori,4)==0 % orientation 12 for example
    temp=permute(cat(2,peak(:,ori,:),peak),[2 3 1]);%peak_p: (0:ori)*ncell*type
    peak_p=interp1(0:ori,temp,0:1/2:ori);
    ori=ori*2;
    if ndims(peak_p)<3
        peak_p=peak_p(2:end,:);
    else
        peak_p=peak_p(2:end,:,:);
    end
else
    peak_p=permute(peak,[2 3 1]);  %peak ori*ncell*type
end

[pref,pref_ori]=max(peak_p,[],1);
orth_ori=mod(pref_ori+ori/4-1,ori)+1;
orth_ori2=mod(pref_ori-ori/4-1,ori)+1;

if ndims(peak_p)<3
    orth=peak_p(sub2ind([ori,ncell],orth_ori,1:ncell));
    orth2=peak_p(sub2ind([ori,ncell],orth_ori2,1:ncell));
else
    orth=peak_p(sub2ind([ori,ncell],orth_ori,1:ncell),:);
    orth2=peak_p(sub2ind([ori,ncell],orth_ori2,1:ncell),:);
end

orth=(orth+orth2)/2;

OSI=(pref-orth)./(pref+orth).*pref;
%OSI=(pref-orth)./(pref+orth);
OSI=permute(OSI,[3 1 2]);