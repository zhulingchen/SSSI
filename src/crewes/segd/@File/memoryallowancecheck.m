function mryflag=memoryallowancecheck(szneeded,type,allowance)
% mryflag=memoryallowancecheck(szneeded,type,allowance)
%
% memoryallowancecheck returns a flag stating whenether or not there is
% enough memory availble to create an array.
%
% Inputs:
%    szneeded is a vector containing the [length width] of the size of the
%      array desired
%    type is a string indicating the type the array will be.  acceptable
%      values are: 'uint8','int8','uint16','int16','uint32','int32',
%                'char','uchar','double','single'
%    allowance is the fraction of the maximum memory the array is being 
%      compared to. allowance=1 wouldbe 100% of the maximum memory
%      available.  Please note that this will cause the operating system to
%      become sluggish.
%      **************************Default 0.25*****************************
%
% Outputs:
%    mryflag is a logical flag indicating if sufficient memory is
%    available.  1 indicates that there is enough memory and 0 indicates
%    that there is not


if nargin<3
    allowance=0.25;
end

[~, sysV]=memory;
mryflag=true;
% get number of bytes each number in array requires
num=1;
numty=typecast(num,type);
sznum=size(numty);sznum=8/sznum(2);

%calculate size of array required
arraysz=szneeded(1)*szneeded(2)*sznum;

%compare total memory to arraysz

if (allowance*sysV.PhysicalMemory.Total)<arraysz
    warndlg('Array requires more memory then is available',...
        'Insufficent Memory');
    mryflag=false;
    return
end

%compare avaliable memory to arraysz
if (allowance*sysV.PhysicalMemory.Available)<arraysz
    response=questdlg(char('Creating Array requires saving data to disk.',...
        'This could make the reading process slower','Would You Like To Continue'),...
        'Approaching Maximum Memory','Yes','No','No');
    if strcmp(response,'Yes')
        mryflag=True;
    else
        mryflag=false;
    end
end


end