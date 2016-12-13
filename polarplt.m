function h=polarplt(peak,Cor);  %peak 1*Var*ncell
%peak n*Var*ncell
color=[ .5 .5 0;.5 0 .5; 0 .5 .5 ;.5 0 0;0 .5 0; 0 0 .5; 0.5,0.5,0.5; 0 0 0];
dim=ndims(peak);
ncell=size(peak,dim);
Ori=size(peak,dim-1);
if ~(dim==3)
    ncolor=1;
    peak=reshape(peak,1,Ori,ncell);
    size(peak)
else
    ncolor=size(peak,1);
end


if ~exist('Cor')
    Cor=1:ncell;
end


xp=[0:2*pi/Ori:2*pi];
yp=cat(2,peak, peak(:,1,:));

%%% 10 cells per row
col=10;
[xpos,ypos,xwidth,yheight]=figurepara(col,ceil(ncell/col));
h=figure('Name','Polarmap for orientation','Position',[ 0 0 200*col 200*ceil(ncell/col)]);

for j=1:ncell
    subplot('position',[xpos(mod(j-1,col)+1) ypos(floor((j-1)/col)+1) xwidth yheight]);
    ymax=max(max(peak(:,:,j)));
    P = polar(xp, ymax * ones(1,Ori+1));
    set(P, 'Visible', 'off');
    hold on;
    for n=1:ncolor
        p2=polar(xp,yp(n,:,j),'--');
        e.Color=color(n,:);
    end
    title(['cell' num2str(Cor(j))]);
    %axis off
end

subplot('position',[xpos(mod(ncell,col)+1) ypos(floor((j-1)/col)+1) xwidth yheight]);hold on
for n=1:ncolor
    p2=plot(1,1,'--');
    e.Color=color(n,:);
end
legend('show','location','northoutside');
axis off;
