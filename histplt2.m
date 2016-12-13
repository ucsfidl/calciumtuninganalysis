function histplt2(data1,data2);
hold on
h1 = histogram(data1(:));
h2=histogram(data2(:));
h1.FaceColor=[ 0.5 0 0];
h2.FaceColor=[0 0 0.5]; 
legend('running','still','Location','south','Orientation','horizontal')
%'Location','northoutside',