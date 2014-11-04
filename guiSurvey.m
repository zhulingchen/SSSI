function varargout = guiSurvey(varargin)
% GUISURVEY MATLAB code for guiSurvey.fig
%      GUISURVEY, by itself, creates a new GUISURVEY or raises the existing
%      singleton*.
%
%      H = GUISURVEY returns the handle to a new GUISURVEY or the handle to
%      the existing singleton*.
%
%      GUISURVEY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUISURVEY.M with the given input arguments.
%
%      GUISURVEY('Property','Value',...) creates a new GUISURVEY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before guiSurvey_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guiSurvey_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guiSurvey

% Last Modified by GUIDE v2.5 04-Nov-2014 12:32:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @guiSurvey_OpeningFcn, ...
    'gui_OutputFcn',  @guiSurvey_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before guiSurvey is made visible.
function guiSurvey_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guiSurvey (see VARARGIN)

% Choose default command line output for guiSurvey
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes guiSurvey wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = guiSurvey_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btn_shot.
function btn_shot_Callback(hObject, eventdata, handles)
% hObject    handle to btn_shot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btn_loadPWaveVelocityModel.
function btn_loadPWaveVelocityModel_Callback(hObject, eventdata, handles)
% hObject    handle to btn_loadPWaveVelocityModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file, path] = uigetfile('*.mat', 'Select a velocity model');
load(fullfile(path, file));

dim = length(size(velocityModel));
vmin = min(velocityModel(:));
vmax = max(velocityModel(:));

% enable edit texts
set(handles.edit_dx, 'Enable', 'on');
set(handles.edit_dz, 'Enable', 'on');
set(handles.edit_dt, 'Enable', 'on');
set(handles.edit_boundary, 'Enable', 'on');
set(handles.ppm_approxOrder, 'Enable', 'on');
set(handles.edit_centerFreq, 'Enable', 'on');
set(handles.edit_sx, 'Enable', 'on');
set(handles.edit_sz, 'Enable', 'on');
set(handles.ppm_sweepAll, 'Enable', 'on');
set(handles.ppm_receiveAll, 'Enable', 'on');
% set finite difference setting values
dx = 10;
dz = 10;
dt = 0.5*(min([dx, dz])/vmax/sqrt(2));
[nz, nx, ~] = size(velocityModel);
nt = round((sqrt((dx*nx)^2 + (dz*nz)^2)*2/vmin/dt + 1));
x = (1:nx) * dx;
z = (1:nz) * dz;

nBoundary = 20;
nDiffOrder = 2;
f = 20;
% set values in edit texts
set(handles.edit_dx, 'String', num2str(dx));
set(handles.edit_dz, 'String', num2str(dz));
set(handles.edit_dt, 'String', num2str(dt));
set(handles.edit_nx, 'String', num2str(nx));
set(handles.edit_nz, 'String', num2str(nz));
set(handles.edit_nt, 'String', num2str(nt));
set(handles.edit_boundary, 'String', num2str(nBoundary));
set(handles.ppm_approxOrder, 'Value', nDiffOrder);
set(handles.edit_centerFreq, 'String', num2str(f));

% 3D case
if (dim > 2)
    % enable edit texts for y-axis
    set(handles.edit_dy, 'Enable', 'on');
    set(handles.edit_sy, 'Enable', 'on');
    % set finite difference setting values
    dy = 10;
    dt = 0.5*(min([dx, dy, dz])/vmax/sqrt(3));
    [~, ~, ny] = size(velocityModel);
    nt = round((sqrt((dx*nx)^2 + (dy*ny)^2 + (dz*nz)^2)*2/vmin/dt + 1));
    % set values in edit texts for y-axis
    set(handles.edit_dy, 'String', num2str(dy));
    set(handles.edit_dt, 'String', num2str(dt));
    set(handles.edit_ny, 'String', num2str(ny));
    set(handles.edit_nt, 'String', num2str(nt));
end

axes(handles.axes_velocityModel);
imagesc(x, z, velocityModel);
xlabel('Distance (m)'); ylabel('Depth (m)');
title('Velocity Model');

% share variables among callback functions
handles.velocityModel = velocityModel;
guidata(hObject, handles);



function edit_dx_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dx as text
%        str2double(get(hObject,'String')) returns contents of edit_dx as a double

velocityModel = handles.velocityModel;
dim = length(size(velocityModel));
vmin = min(velocityModel(:));
vmax = max(velocityModel(:));

dx = str2double(get(handles.edit_dx, 'String'));
dz = str2double(get(handles.edit_dz, 'String'));
dt = 0.5*(min([dx, dz])/vmax/sqrt(2));
[nz, nx, ~] = size(velocityModel);
nt = round((sqrt((dx*nx)^2 + (dz*nz)^2)*2/vmin/dt + 1));
set(handles.edit_dt, 'String', num2str(dt));
set(handles.edit_nt, 'String', num2str(nt));

% 3D case
if (dim > 2)
    dt = 0.5*(min([dx, dy, dz])/vmax/sqrt(3));
    [~, ~, ny] = size(velocityModel);
    nt = round((sqrt((dx*nx)^2 + (dy*ny)^2 + (dz*nz)^2)*2/vmin/dt + 1));
    set(handles.edit_dt, 'String', num2str(dt));
    set(handles.edit_nt, 'String', num2str(nt));
end


% --- Executes during object creation, after setting all properties.
function edit_dx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dy_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dy as text
%        str2double(get(hObject,'String')) returns contents of edit_dy as a double


% --- Executes during object creation, after setting all properties.
function edit_dy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dz_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dz as text
%        str2double(get(hObject,'String')) returns contents of edit_dz as a double

velocityModel = handles.velocityModel;
dim = length(size(velocityModel));
vmin = min(velocityModel(:));
vmax = max(velocityModel(:));

dx = str2double(get(handles.edit_dx, 'String'));
dz = str2double(get(handles.edit_dz, 'String'));
dt = 0.5*(min([dx, dz])/vmax/sqrt(2));
[nz, nx, ~] = size(velocityModel);
nt = round((sqrt((dx*nx)^2 + (dz*nz)^2)*2/vmin/dt + 1));
set(handles.edit_dt, 'String', num2str(dt));
set(handles.edit_nt, 'String', num2str(nt));

% 3D case
if (dim > 2)
    dt = 0.5*(min([dx, dy, dz])/vmax/sqrt(3));
    [~, ~, ny] = size(velocityModel);
    nt = round((sqrt((dx*nx)^2 + (dy*ny)^2 + (dz*nz)^2)*2/vmin/dt + 1));
    set(handles.edit_dt, 'String', num2str(dt));
    set(handles.edit_nt, 'String', num2str(nt));
end


% --- Executes during object creation, after setting all properties.
function edit_dz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dt_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dt as text
%        str2double(get(hObject,'String')) returns contents of edit_dt as a double


% --- Executes during object creation, after setting all properties.
function edit_dt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_nx_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nx as text
%        str2double(get(hObject,'String')) returns contents of edit_nx as a double


% --- Executes during object creation, after setting all properties.
function edit_nx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ny_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ny as text
%        str2double(get(hObject,'String')) returns contents of edit_ny as a double


% --- Executes during object creation, after setting all properties.
function edit_ny_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_nz_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nz as text
%        str2double(get(hObject,'String')) returns contents of edit_nz as a double


% --- Executes during object creation, after setting all properties.
function edit_nz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_nt_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nt as text
%        str2double(get(hObject,'String')) returns contents of edit_nt as a double


% --- Executes during object creation, after setting all properties.
function edit_nt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function menuHelp_Callback(hObject, eventdata, handles)
% hObject    handle to menuHelp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuHelpAbout_Callback(hObject, eventdata, handles)
% hObject    handle to menuHelpAbout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guiAbout;



function edit_boundary_Callback(hObject, eventdata, handles)
% hObject    handle to edit_boundary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_boundary as text
%        str2double(get(hObject,'String')) returns contents of edit_boundary as a double


% --- Executes during object creation, after setting all properties.
function edit_boundary_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_boundary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ppm_approxOrder.
function ppm_approxOrder_Callback(hObject, eventdata, handles)
% hObject    handle to ppm_approxOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ppm_approxOrder contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ppm_approxOrder


% --- Executes during object creation, after setting all properties.
function ppm_approxOrder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ppm_approxOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_centerFreq_Callback(hObject, eventdata, handles)
% hObject    handle to edit_centerFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_centerFreq as text
%        str2double(get(hObject,'String')) returns contents of edit_centerFreq as a double


% --- Executes during object creation, after setting all properties.
function edit_centerFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_centerFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_sx_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sx as text
%        str2double(get(hObject,'String')) returns contents of edit_sx as a double


% --- Executes during object creation, after setting all properties.
function edit_sx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_sy_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sy as text
%        str2double(get(hObject,'String')) returns contents of edit_sy as a double


% --- Executes during object creation, after setting all properties.
function edit_sy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_sz_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sz as text
%        str2double(get(hObject,'String')) returns contents of edit_sz as a double


% --- Executes during object creation, after setting all properties.
function edit_sz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ppm_sweepAll.
function ppm_sweepAll_Callback(hObject, eventdata, handles)
% hObject    handle to ppm_sweepAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ppm_sweepAll contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ppm_sweepAll

velocityModel = handles.velocityModel;
dim = length(size(velocityModel));

isSweepAll = get(hObject, 'Value');
if (isSweepAll == 1) % Yes
    set(handles.edit_sx, 'Enable', 'off');
    set(handles.edit_sz, 'Enable', 'off');
end
if (isSweepAll == 2) % No
    set(handles.edit_sx, 'Enable', 'on');
    set(handles.edit_sz, 'Enable', 'on');
end

% 3D case
if (dim > 2)
    if (isSweepAll == 1) % Yes
        set(handles.edit_sy, 'Enable', 'off');
    end
    if (isSweepAll == 2) % No
        set(handles.edit_sy, 'Enable', 'on');
    end
end


% --- Executes during object creation, after setting all properties.
function ppm_sweepAll_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ppm_sweepAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ppm_receiveAll.
function ppm_receiveAll_Callback(hObject, eventdata, handles)
% hObject    handle to ppm_receiveAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ppm_receiveAll contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ppm_receiveAll

velocityModel = handles.velocityModel;
[nz, nx, ~] = size(velocityModel);
dim = length(size(velocityModel));

isReceiveAll = get(hObject, 'Value');
if (isReceiveAll == 1) % Yes
    set(handles.edit_rx, 'Enable', 'off');
    set(handles.edit_rz, 'Enable', 'off');
end
if (isReceiveAll == 2) % No
    set(handles.edit_rx, 'Enable', 'on');
    set(handles.edit_rz, 'Enable', 'on');
    str_rx = sprintf('1:%d', nx);
    str_rz = sprintf('%d', 1);
    set(handles.edit_rx, 'String', str_rx);
    set(handles.edit_rz, 'String', str_rz);
end

% 3D case
if (dim > 2)
    [~, ~, ny] = size(velocityModel);
    if (isReceiveAll == 1) % Yes
        set(handles.edit_ry, 'Enable', 'off');
    end
    if (isReceiveAll == 2) % No
        set(handles.edit_ry, 'Enable', 'on');
        str_ry = sprintf('1:%d', ny);
        set(handles.edit_ry, 'String', str_ry);
    end
end


% --- Executes during object creation, after setting all properties.
function ppm_receiveAll_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ppm_receiveAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_rx_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rx as text
%        str2double(get(hObject,'String')) returns contents of edit_rx as a double


% --- Executes during object creation, after setting all properties.
function edit_rx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ry_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ry as text
%        str2double(get(hObject,'String')) returns contents of edit_ry as a double


% --- Executes during object creation, after setting all properties.
function edit_ry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_rz_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rz as text
%        str2double(get(hObject,'String')) returns contents of edit_rz as a double


% --- Executes during object creation, after setting all properties.
function edit_rz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
