function b = MSHwriter(ptCoordMx,faceMx,faceMxSize,ptCoordMxSize,filename )
% filename = 'ianCoWv2.msh';   
b = 1;
d = faceMxSize.in;
c = max(faceMx.in(d-20:d,5:6));
allFace = faceMxSize.in+faceMxSize.inlet+faceMxSize.sur;
fid1 = fopen(filename,'w');
fprintf(fid1, '(0 " Created by : Fluent_V6 Interface Vers. 14.0.3")\n');
% fprintf(fid1, '(0 "Structured mesh for vascular network - Model Generator Mahsa Ghaffari, 4/21/2014")\n\n');
fprintf(fid1, '(2 3)\n');
fprintf(fid1, '(0 "Node Section")\n');
fprintf(fid1, '(10 (0 1 %X 0 3))\n',ptCoordMxSize);
fprintf(fid1, '(10 (f 1 %X 1 3)\n(\n',ptCoordMxSize);
for i = 1:ptCoordMxSize   
    fprintf(fid1, '%d %d %d\n',ptCoordMx(i,1),ptCoordMx(i,2),ptCoordMx(i,3));
end
fprintf(fid1, '))\n');
fprintf(fid1, '(12 (0 1 %X 0 0))\n',c(1)); 
fprintf(fid1, '(12 (10 1 %X 1 4))\n',c(1));
fprintf(fid1, '(13 (0 1 %X 0 0))\n',(allFace)); 
fprintf(fid1, '(0 "Interior faces of zone VESSEL_INTERIOR")\n');
fprintf(fid1, '(13(11 1 %X 2 4)(\n',faceMxSize.in);
for iFace = 1:faceMxSize.in
    fprintf(fid1, '%X %X %X %X %X %X\n',faceMx.in(iFace,1), faceMx.in(iFace,2), faceMx.in(iFace,3), faceMx.in(iFace,4),faceMx.in(iFace,5), faceMx.in(iFace,6));
end
fprintf(fid1, ')\n)\n');
fprintf(fid1, '(0 "Faces of zone INLET")\n');
fprintf(fid1, '(13 (12 %X %X a 4)(\n',faceMxSize.in+1,(faceMxSize.in+faceMxSize.inlet));
for iFace = 1:length(faceMx.inlet)
    fprintf(fid1, '%X %X %X %X %X %X\n',faceMx.inlet(iFace,1), faceMx.inlet(iFace,2), faceMx.inlet(iFace,3), faceMx.inlet(iFace,4),faceMx.inlet(iFace,5), faceMx.inlet(iFace,6));
end
fprintf(fid1, ')\n)\n');
% fprintf(fid1, '(0 "Faces of zone OUTLET")\n');
% fprintf(fid1, '(13 (13 %X %X 5 4)(\n',(length(faceMx.in)+length(faceMx.inlet))+1,(length(faceMx.in)+length(faceMx.inlet))+length(faceMx.outlet));
% for iFace = 1:length(faceMx.outlet)
%     fprintf(fid1, '%X %X %X %X %X %X\n',faceMx.outlet(iFace,1), faceMx.outlet(iFace,2), faceMx.outlet(iFace,3), faceMx.outlet(iFace,4),faceMx.outlet(iFace,5), faceMx.outlet(iFace,6));
% end
% fprintf(fid1, ')\n)\n');
fprintf(fid1, '(0 "Faces of zone VESSEL_WALL")\n');
fprintf(fid1, '(13 (13 %X %X 3 4)(\n',(faceMxSize.in+faceMxSize.inlet)+1,allFace);
for iFace = 1:faceMxSize.sur
    fprintf(fid1, '%X %X %X %X %X %X\n',faceMx.sur(iFace,1), faceMx.sur(iFace,2), faceMx.sur(iFace,3), faceMx.sur(iFace,4),faceMx.sur(iFace,5), faceMx.sur(iFace,6));
end
fprintf(fid1, ')\n)\n');


fprintf(fid1, '(0 "Zone Sections")\n');
fprintf(fid1, '(39 (16 fluid VESSEL_INTERIOR)())\n');
fprintf(fid1, '(39 (17 interior int_VESSEL_INTERIOR)())\n');
fprintf(fid1, '(39 (18 velocity-inlet INLET)())\n');
% fprintf(fid1, '(39 (19 outlet-vent OUTLET)())\n');
% fprintf(fid1, '(39 (20 wall VESSEL_WALL)())\n');
fprintf(fid1, '(39 (19 wall VESSEL_WALL)())\n');
fclose(fid1);