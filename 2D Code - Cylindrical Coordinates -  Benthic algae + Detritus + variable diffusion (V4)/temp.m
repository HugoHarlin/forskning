
colorvals = ones(p.Zn-1,p.Xn-1,1);
for j = 1:p.Xn-1
    for i =1: p.Zn-1
        colorvals(i,j,1) = D(i,j);
    end
end


 surf(p.X,p.Z,zeros(p.Xn), colorvals,'edgecolor','none');