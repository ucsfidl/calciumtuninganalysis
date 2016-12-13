function sbxsplit(fname)
sbxread(fn,1,1);        % read one frame to read the header of the image sequence
global info;
if info.volscan=1
            info.fid = fopen([fname '.sbx']);
        d = dir([fname '.sbx']);
        info.nsamples = (info.sz(2) * info.recordsPerBuffer * 2 * info.nchan);   % bytes per record
    for     
    for nplane=1:info.otparam(3)

       
        try
            fseek(info.fid,k*info.nsamples,'bof');
            x = fread(info.fid,info.nsamples/2 * N,'uint16=>uint16');

    end
end