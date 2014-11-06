function varargout = guiAbout(varargin)
% GUIABOUT MATLAB code for guiAbout.fig
%      GUIABOUT, by itself, creates a new GUIABOUT or raises the existing
%      singleton*.
%
%      H = GUIABOUT returns the handle to a new GUIABOUT or the handle to
%      the existing singleton*.
%
%      GUIABOUT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIABOUT.M with the given input arguments.
%
%      GUIABOUT('Property','Value',...) creates a new GUIABOUT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before guiAbout_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guiAbout_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guiAbout

% Last Modified by GUIDE v2.5 04-Nov-2014 00:07:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guiAbout_OpeningFcn, ...
                   'gui_OutputFcn',  @guiAbout_OutputFcn, ...
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


% --- Executes just before guiAbout is made visible.
function guiAbout_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guiAbout (see VARARGIN)

% Choose default command line output for guiAbout
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% display logos
imshow('logo_csip.jpg', 'Parent', handles.axes_logo1);

imshow('logo_cegp.jpg', 'Parent', handles.axes_logo2);

% UIWAIT makes guiAbout wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = guiAbout_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btnOK.
function btnOK_Callback(hObject, eventdata, handles)
% hObject    handle to btnOK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.figure1);
