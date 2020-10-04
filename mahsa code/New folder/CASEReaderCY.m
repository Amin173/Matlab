%% C.Y. Hsu & T.J. Marrinan BIOE310
function [diameterVector, ptCoordMx, faceMx, groupMx, pointMx,viscosity,P1,PN] = CASEReaderCY(fileName)
% tic
Header = 0;
Constants = 1;
Dia = 2;
other = 3;
fid = fopen(fileName);
state = Header;
while ~feof(fid)
    switch state
        case Header
            nwkFileName = readNWKFileHeader(fid);
            state = Constants;
        case Constants
            [sysViscosity, P1, PN] = readConstants(fid);
            state = Dia;
        case Dia                        
            [ptCoordMx, faceMx, groupMx, pointMx] = NWKReaderCY(nwkFileName);
            diameterVector = readDiameterVector(fid,groupMx);             
            state = other;
        case other
            fgetl(fid);
    end    
end
fclose(fid);
viscosity = computeViscosity(diameterVector, sysViscosity);
% toc
% memory
%%
function nwkFileName = readNWKFileHeader(fid)
foundHeader = 0;
while foundHeader == 0
   line = fgetl(fid);
   if strncmp(line,'meshfile=',9) == 1
       nwkFileName = sscanf(line,'meshfile=%s');
       foundHeader = 1;
   end
end
%%
function [sysViscosity, P1, PN] = readConstants(fid)
line = fgetl(fid);
if strncmp(line,'<const',6) == 1
    line = fgetl(fid);
    sysViscosity = 0;
    P1 = 0;
    PN = 0;
    while strncmp(line,'</const',7) ~= 1
        if strncmp(line,'MY=',3) == 1
            sysViscosity = sscanf(line,'MY=%s');
            sysViscosity = str2double(sysViscosity);
        elseif strncmp(line,'P1=',3) == 1
            P1 = sscanf(line,'P1=%d');
        elseif strncmp(line,'PN=',3) == 1
            PN = sscanf(line,'PN=%d');
        end
        line = fgetl(fid);
    end
end

%%
function diameterVector = readDiameterVector(fid,groupMx)
for iGroup = 1:size(groupMx,1)
    [fromIdx toIdx foundHeader value] = readDiameterHeader(fid);
    if foundHeader == 1
        for i = fromIdx:toIdx
            line = fgetl(fid);
            diameterVector(i,1) = str2num(sscanf(line,'%s'));
        end
    elseif foundHeader == 2
        diameterVector(fromIdx:toIdx,1) = value;    
    end
end
%%
function [fromIdx toIdx foundHeader value] = readDiameterHeader(fid)
foundHeader = 0;
while foundHeader == 0
    line = fgetl(fid);
    if strncmp(line, '(vector=Dia',11) == 1
        groupID = sscanf(line,'(vector=Dia group=%x');
        line = fgetl(fid);
        if strfind(line,'type=variable') ~= 0
            NFaces = sscanf(line,'(N=%d type=variable)');
            line = fgetl(fid);
            groupHeader = sscanf(line,'(13(%x %x %x %*x %*x)(');
            fromIdx = groupHeader(2);
            toIdx = groupHeader(3);
            if groupHeader(1) ~= groupID 
                error('GroupID Mismatch');
            end
            if toIdx-fromIdx+1 ~= NFaces
                error('NumberOfFaces Mismatch');                
            end
            foundHeader = 1;
            value = 0;
        elseif strfind(line, 'type=const') ~= 0
            NFaces = sscanf(line,'(N=%d type=const)');
            value = sscanf(line,'(N=%*d type=const value=%f)');
            line = fgetl(fid);
            groupHeader = sscanf(line,'(13(%x %x %x %*x %*x)(');
            fromIdx = groupHeader(2);
            toIdx = groupHeader(3);         
            foundHeader = 2;
        end
    end
end
%%
function viscosity = computeViscosity(dd,sysViscosity)
aD = dd*1e6; 
beta = 1./((1+aD.^12*10^-11));
C = (0.8+exp(-0.075*aD)).*(-1+beta)+beta;
MYpt45 = 3.2+220*exp(-1.3*aD)-2.44*exp(-0.06*aD.^0.645);
num = (1-sysViscosity).^C-1; 
den = (1-0.45).^C-1;
viscosity = (1+(MYpt45-1).*num./den)*0.001;
