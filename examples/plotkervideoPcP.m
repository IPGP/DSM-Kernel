clear;
figure(1);
clf reset;
format long; 
maxk= 1.e-25;
%wysiwyg;

HandleTitle = title('Time=');
set(HandleTitle,'Erase','xor');
   
for t=7000:7801
filename=strcat('pcp80/',sprintf('%05d',t),'.01300.trikernel')
kernel=load(filename);
coordrs=load('coordrspcp');
x=kernel(:,2);
y=kernel(:,1);
k=(kernel(:,4));

% Locations of source and receiver in a vertical view.

ndata=length(x);
nview=coordrs(1);
mx=coordrs(2);
my=coordrs(3);
srad=coordrs(4)/6371;
distan=coordrs(5);

if(nview==1)
%  y=180-y;
elseif(nview==2)
  r=y/6371;
  p=x;
  clear x;
  clear y;
  fan0=p(mx)-p(1);
  delta=(fan0/2+p(1))*2;
  p=p-delta/2;
  p=p-90;
  p=-p;
  for i=1:ndata
    x(i)=r(i)*cos(pi*p(i)/180.);
    y(i)=r(i)*sin(pi*p(i)/180.);
  end
  r670=(6371-670)/6371;
  r410=(6371-410)/6371;
  rcmb=3480/6371;
  rmoho=6340/6371;
  for i=1:mx
    xmoho(i)=rmoho*cos(pi*p(i)/180.);
    ymoho(i)=rmoho*sin(pi*p(i)/180.);
  end
  for i=1:mx
    x670(i)=r670*cos(pi*p(i)/180.);
    y670(i)=r670*sin(pi*p(i)/180.);
  end
  for i=1:mx
    x410(i)=r410*cos(pi*p(i)/180.);
    y410(i)=r410*sin(pi*p(i)/180.);
  end
  for i=1:mx
    xcmb(i)=rcmb*cos(pi*p(i)/180.);
    ycmb(i)=rcmb*sin(pi*p(i)/180.);
  end
  rtop=r(ndata);
  rbottom=r(1);
  pmax=p(ndata);
  pmin=p(1);
  for i=1:mx
    xtop(i)=rtop*cos(pi*p(i)/180.);
    ytop(i)=rtop*sin(pi*p(i)/180.);
  end
  for i=1:mx
    xbottom(i)=rbottom*cos(pi*p(i)/180.);
    ybottom(i)=rbottom*sin(pi*p(i)/180.);
  end
  xr(1)=rtop*cos(pmax*pi/180.);
  yr(1)=rtop*sin(pmax*pi/180.);
  xr(2)=rbottom*cos(pmax*pi/180.);
  yr(2)=rbottom*sin(pmax*pi/180.);
  xl(1)=rtop*cos(pmin*pi/180.);
  yl(1)=rtop*sin(pmin*pi/180.);
  xl(2)=rbottom*cos(pmin*pi/180.);
  yl(2)=rbottom*sin(pmin*pi/180.);
  xsr(1)=srad*cos((90+delta/2)*pi/180.);
  ysr(1)=srad*sin((90+delta/2)*pi/180.);
  xsr(2)=rtop*cos((90-delta/2)*pi/180.);
  ysr(2)=rtop*sin((90-delta/2)*pi/180.);
elseif(nview==3)
  r=y/6371;
  p=180-x;
  clear x;
  clear y;
  for i=1:ndata
    x(i)=r(i)*cos(pi*p(i)/180.);
    y(i)=r(i)*sin(pi*p(i)/180.);
  end 
  r670=(6371-670)/6371;
  r410=(6371-410)/6371;
  rcmb=3480/6371;
  rmoho=6340/6371;
  for i=1:mx
    xmoho(i)=rmoho*cos(pi*p(i)/180.);
    ymoho(i)=rmoho*sin(pi*p(i)/180.);
  end
  for i=1:mx
    x670(i)=r670*cos(pi*p(i)/180.);
    y670(i)=r670*sin(pi*p(i)/180.);
  end
  for i=1:mx
    x410(i)=r410*cos(pi*p(i)/180.);
    y410(i)=r410*sin(pi*p(i)/180.);
  end
  for i=1:mx
    xcmb(i)=rcmb*cos(pi*p(i)/180.);
    ycmb(i)=rcmb*sin(pi*p(i)/180.);
  end
  rtop=r(ndata);
  rbottom=r(1);
  pmax=p(ndata);
  pmin=p(1);
  for i=1:mx
    xtop(i)=rtop*cos(pi*p(i)/180.);
    ytop(i)=rtop*sin(pi*p(i)/180.);
  end
  for i=1:mx
    xbottom(i)=rbottom*cos(pi*p(i)/180.);
    ybottom(i)=rbottom*sin(pi*p(i)/180.);
  end
  xr(1)=rtop*cos(pmax*pi/180.);
  yr(1)=rtop*sin(pmax*pi/180.);
  xr(2)=rbottom*cos(pmax*pi/180.);
  yr(2)=rbottom*sin(pmax*pi/180.);
  xl(1)=rtop*cos(pmin*pi/180.);
  yl(1)=rtop*sin(pmin*pi/180.);
  xl(2)=rbottom*cos(pmin*pi/180.);
  yl(2)=rbottom*sin(pmin*pi/180.);
  xsr(1)=0;
  ysr(1)=srad;
  xsr(2)=xsr(1);
  ysr(2)=1;
end

% Radii to identify map-view plots.

clim1=[-maxk maxk];

if(nview==1)
  postk1=[0.1 0.1 0.8 0.8];
  postb1=[0.3 0.2 0.4 0.03];
elseif(nview==2)
  postk1=[0.075 0.3 0.4 0.4];
  postb=[0.325 0.3 0.4 0.02];
  postk2=[0.525 0.3 0.4 0.4];
elseif(nview==3)
  postk1=[0.075 0.3 0.4 0.4];
  postb=[0.35 0.35 0.3 0.02];
  postk2=[0.525 0.3 0.4 0.4];
end

% First plot.


%subplot(2,2,1);
    axis('off');
c=k;
for i=1:my
  m1=(i-1)*mx+1;
  m2=i*mx;
  X(i,:)=[x(m1:m2)'];
  Y(i,:)=[y(m1:m2)'];
  CP(i,:)=[c(m1:m2)'];
end
axes('position',postk1,'visible','off');
pcolor(X,Y,CP);
axis('off');
colormap(karason(128,0));
%colormap([0 0 191;0 255 255 216:191 0 0]);
%colormap(karason);
caxis(clim1);
shading flat;
axis('equal');
hold on;
if(nview==1)
  h3=colorbar('horiz');
  set(h3,'Position',postb1,'Fontsize',11);
elseif(nview==2)
   
  plot(xtop,ytop,'k-');
  plot(xbottom,ybottom,'k-'); 
  plot(xl,yl,'k-');
  plot(xr,yr,'k-');
  plot(xmoho,ymoho,'k--','linewidth',0.4);
  plot(x410,y410,'k--','linewidth',0.4);
  plot(x670,y670,'k--','linewidth',0.4);
  plot(xsr(1),ysr(1),'k*','markersize',3,'linewidth',3);
  plot(xsr(2),ysr(2),'kv','markersize',3,'linewidth',3);
  h3=colorbar('horiz');
  set(h3,'Position',postb,'Fontsize',8);
  %set(h3,'Position',
  %delete (strtmp);
  % set(figure(1),'Title', 'visible','off');
  strtmp=strcat('Time=',sprintf('%5.1f',t/10), '  sec');
  set(HandleTitle,'String',strtmp);
  legend off;
  %title(strtmp);
  %set(gca,'Title',text('String',strtmp));
  %title on;
  
else
  h3=colorbar('horiz');
  set(h3,'Position',postb,'Fontsize',11);
  plot(xtop,ytop,'k-');
  plot(xbottom,ybottom,'k-');
  plot(xl,yl,'k-');
  plot(xr,yr,'k-');
  plot(xsr(1),ysr(1),'k*','markersize',3,'linewidth',3);
  plot(xsr(2),ysr(2),'kv','markersize',3,'linewidth',3);
end
 filename=strcat('pcp80/jpg/',sprintf('%05d',t),'.00080.jpg') ;
 set(figure(1),'PaperPosition');
 saveas(figure(1),filename,'jpg');
 %filename=strcat('triplication23/eps/',sprintf('%05d',t),'.00230.eps') ;
 %saveas(figure(1),filename,'eps');
%drawnow
%pause(0.02)
%movmov(t)=getframe;

end

%
%movie2avi(movmov,'triplicationsensitivity.avi');