points = 8;
sections = 40;

kamdist=kam.output_data.distr_original;

kam.plot('dist','original');
ylim([Dps(1) Dps(end)]);
caxis([1 4]);
gcc=findall(gcf,'tag','Colorbar');
delete(gcc);

[x, y] = ginput(points)
for i=1:points
A=find(kamdist(:,1) >= x(i));
x_indices(i) = A(1);
end

for i=1:points
A=find(kamdist(x_indices(i),3+sections:end)>=y(i));
y_indices(i)=A(1);
end

xs=x_indices(1):length(kamdist(:,1));
ys=interp1(x_indices,y_indices,xs);
ys=round(ys);
B=find(~isnan(ys));
ys(isnan(ys))=ys(B(end));
dist_new = kamdist;
for i=1:length(xs)
dist_new(xs(i),3:ys(i)+1)=0;
end