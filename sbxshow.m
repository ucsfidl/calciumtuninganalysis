function sbxshow(fname,idx)
global info;
h=figure('Position', [100, 100, 2000, 2000])

for(i=1:length(idx))
    z = sbxread(fname,idx(i),1);
    z = squeeze(z);
    if ~isempty(info.aligned)
        z = circshift(z,info.aligned.T(i+1,:)); % align the image
    end
    imshow(z(:,:),'InitialMagnification',200);
    axis on;
end