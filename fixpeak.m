function peak=fixpeak(peak);
global info;

[m_ctr,m_ori]=find(isnan(peak(:,:,1)));
for k=1:numel(m_ori)
    before=m_ori(k)-1;before=mod(before-1,info.steps(1))+1;
    after=m_ori(k)+1;after=mod(after-1,info.steps(1))+1;
    peak(m_ctr(k),m_ori(k),:)= (peak(m_ctr,before,:)+peak(m_ctr,after,:))/2;
end