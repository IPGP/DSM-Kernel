filename=('kernelp');

kernel=load(filename);
coordrs=load('coordrs');
x=kernel(:,2);
y=kernel(:,1);
k=(kernel(:,5));
nview=coordrs(1);

if(nview==1)
[xi,yi]=meshgrid(-4.:.1:150.,-20.:.1:20.);
zi=griddata(x,y,k,xi,yi);
end
if(nview==2)
[xi,yi]=meshgrid(-4.:.1:150.,3551.:40.:6371.);
zi=griddata(x,y,k,xi,yi);
end

xii=xi(:);
yii=yi(:);
zii=zi(:);

fid=fopen('reshaped.kernel','w');
fprintf(fid,'%e %e %e \n', [xii';yii';zii']);
fclose(fid);
