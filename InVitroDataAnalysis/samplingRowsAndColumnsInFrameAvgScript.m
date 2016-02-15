
frameAvg = frameAverage('1543',1,100,1,0);

figure;
for i=1:10
row = 10+(i-1)*50;
subplot(3,4,i); plot(frameAvg(row,:)); 
xlabel('horiz position (pixel)');
ylabel('I in frame average');
title('sampling rows in frame average');
end

figure;
for i=1:10
column = 10+(i-1)*50;
subplot(3,4,i); plot(frameAvg(:,column)); 
xlabel('vertical position (pixel)');
ylabel('I in frame average');
title('sampling columns in frame average');
end