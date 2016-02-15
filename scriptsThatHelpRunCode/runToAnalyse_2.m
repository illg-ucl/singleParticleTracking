
s498 = FindTrajects('498',7,500);
linkTrajSegments('498',7,500,s498,'cybD-mCherry-ATPase-GFp');

s500 = FindTrajects('500',5,500);
linkTrajSegments('500',5,500,s500,'cybD-mCherry-ATPase-GFp');

s509 = FindTrajects('509',22,500);
linkTrajSegments('509',22,500,s509,'cybD-mCherry-ATPase-GFp');

s513 = FindTrajects('513',11,500);
linkTrajSegments('513',11,500,s513,'cybD-mCherry-ATPase-GFp');

s515 = FindTrajects('515',10,500);
linkTrajSegments('515',10,500,s515,'cybD-mCherry-ATPase-GFp');

s518 = FindTrajects('518',7,500);
linkTrajSegments('518',7,500,s518,'cybD-mCherry-ATPase-GFp');

s522 = FindTrajects('522',2,500);
linkTrajSegments('522',2,500,s522,'cybD-mCherry-ATPase-GFp');

s524 = FindTrajects('524',20,500);
linkTrajSegments('524',20,500,s524,'cybD-mCherry-ATPase-GFp');

save 'resultStructures' 's*' % save all result structures in a .mat file.