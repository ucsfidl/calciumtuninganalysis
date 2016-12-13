function bad=pick(f);
    %hcell=figure_enlarge(hcell,2);
    data=loadjson( [f(1:findstr(f,'cell')+3) '.jmesh']);
    circle=data.jmesh{:};
    total=numel(circle);
    
    for n=1:total
        Co(:,n)=circle{n}.bbox;
    end
    %%%%%%Define bad signals: edge and background
    ratio=[Co(2,:)-Co(1,:)]./[Co(4,:)-Co(3,:)];
    long=find(ratio>3|ratio<1/3);
    edge=([Co(3,:)+Co(3,:)])/2;
    out=find(edge<40|edge>720);
    defaultbad=cellstr(num2str([long,out]));
    
    selected=inputdlg('input bad cell#','Exclusion',[4 80],defaultbad);
    bad=unique(str2num(selected{1}));
    
%    S_df=S_df(Cor);
   
