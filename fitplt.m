function h=fitplt(peak,error,Cor);
%peak n*Var*ncell
color=[ 1 0 0;   %red
        0 1 0;     %green
        0 0 1 ;     %blue
        .5 .5 .5;   %grey
        .5 0 0;
        0 .5 0;
        0 0 .5 ;
        .25 .25 0.25];

dim=ndims(peak);
ncell=size(peak,dim);
Var=size(peak,dim-1);
type=size(peak,dim-2);

if ~(dim==3)
       peak=reshape(peak,[],Var,ncell);
end


if ~exist('Cor') 
    Cor=1:ncell;
end

if isempty(error)
    error=zeros(type,Var,ncell);
end

x=1:Var;
xx=1:.1:Var;
temp=shiftdim(peak,1); %temp Var*ncell*n
temp=reshape(temp,Var,ncell,type);
temp2=interp1(x,temp,xx); %temp2 length*ncell*n
yy=shiftdim(temp2,2);  %yy n*length*ncell
yy=reshape(yy,type,[],ncell);
%10/per row
col=10;
[xpos,ypos,xwidth,yheight]=figurepara(col,ceil(ncell/col));

h=figure('Name','Linear map','Position',[ 200 100 200*col 100*ceil(ncell/col)]);

for j=1:ncell
    subplot('position',[ xpos(mod(j-1,col)+1) ypos(floor((j-1)/col)+1)   xwidth yheight]);hold on;
    temp=peak(:,:,j);
    ymax=max(temp(:));
    for n=1:type
        %e=errorbar(x,peak(n,:,j),error(n,:,j));
        e=plot(x,peak(n,:,j),'o');
       plot(xx,yy(n,:,j),'--');
        e.Color=color(n,:);
    end
    axis([ 1 Var 0 max(0.1,ymax)])
    a=gca;
    set(a,'xticklabel',[]);
    title(['cell' num2str(Cor(j))]);
    %axis off
end

%%legend
subplot('position',[xpos(mod(ncell,col)+1) ypos(floor((j-1)/col)+1) xwidth yheight]);hold on
for n=1:type
    p2=plot(1,1,'-');
    p2.Color=color(n,:);
end
legend('show');
axis off;


