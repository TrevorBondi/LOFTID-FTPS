function varargout = LOFTID_FTPS_GUI(varargin)
% LOFTID_FTPS_GUI MATLAB code for LOFTID_FTPS_GUI.fig
%      LOFTID_FTPS_GUI, by itself, creates a new LOFTID_FTPS_GUI or raises the existing
%      singleton*.
%
%      H = LOFTID_FTPS_GUI returns the handle to a new LOFTID_FTPS_GUI or the handle to
%      the existing singleton*.
%
%      LOFTID_FTPS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOFTID_FTPS_GUI.M with the given input arguments.
%
%      LOFTID_FTPS_GUI('Property','Value',...) creates a new LOFTID_FTPS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LOFTID_FTPS_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LOFTID_FTPS_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LOFTID_FTPS_GUI

% Last Modified by GUIDE v2.5 01-Aug-2019 13:07:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LOFTID_FTPS_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @LOFTID_FTPS_GUI_OutputFcn, ...
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


% --- Executes just before LOFTID_FTPS_GUI is made visible.
function LOFTID_FTPS_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LOFTID_FTPS_GUI (see VARARGIN)

% Choose default command line output for LOFTID_FTPS_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes LOFTID_FTPS_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LOFTID_FTPS_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function torus_n_text_Callback(hObject, eventdata, handles)
% hObject    handle to torus_n_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of torus_n_text as text
%        str2double(get(hObject,'String')) returns contents of torus_n_text as a double
Run_button_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function torus_n_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to torus_n_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function d1_text_Callback(hObject, eventdata, handles)
% hObject    handle to d1_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d1_text as text
%        str2double(get(hObject,'String')) returns contents of d1_text as a double
Run_button_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function d1_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to d1_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function d2_text_Callback(hObject, eventdata, handles)
% hObject    handle to d2_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d2_text as text
%        str2double(get(hObject,'String')) returns contents of d2_text as a double
Run_button_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function d2_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to d2_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function dp_text_Callback(hObject, eventdata, handles)
% hObject    handle to dp_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dp_text as text
%        str2double(get(hObject,'String')) returns contents of dp_text as a double
Run_button_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function dp_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dp_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function angle_c_deg_text_Callback(hObject, eventdata, handles)
% hObject    handle to angle_c_deg_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of angle_c_deg_text as text
%        str2double(get(hObject,'String')) returns contents of angle_c_deg_text as a double
Run_button_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function angle_c_deg_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angle_c_deg_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function angle_s_deg_vertical_text_Callback(hObject, eventdata, handles)
% hObject    handle to angle_s_deg_vertical_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of angle_s_deg_vertical_text as text
%        str2double(get(hObject,'String')) returns contents of angle_s_deg_vertical_text as a double
Run_button_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function angle_s_deg_vertical_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angle_s_deg_vertical_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ply_text_Callback(hObject, eventdata, handles)
% hObject    handle to ply_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ply_text as text
%        str2double(get(hObject,'String')) returns contents of ply_text as a double
Run_button_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function ply_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ply_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function c_factor_text_Callback(hObject, eventdata, handles)
% hObject    handle to c_factor_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c_factor_text as text
%        str2double(get(hObject,'String')) returns contents of c_factor_text as a double
Run_button_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function c_factor_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c_factor_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SiC_text_Callback(hObject, eventdata, handles)
% hObject    handle to SiC_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SiC_text as text
%        str2double(get(hObject,'String')) returns contents of SiC_text as a double
Run_button_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function SiC_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SiC_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in geometry_popup.
function geometry_popup_Callback(hObject, eventdata, handles)
% hObject    handle to geometry_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns geometry_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from geometry_popup
Run_button_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function geometry_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to geometry_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Run_button.
function Run_button_Callback(hObject, eventdata, handles)
% hObject    handle to Run_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Poll GUI for inputs and change values if outside of working bounds
clc;
cla(handles.cross_section_plot,'reset');

n = abs(floor(str2double(get(handles.torus_n_text,'string'))));
if n < 2
    n = 2;
end
set(handles.torus_n_text,'string', n)

d1 = abs(str2double(get(handles.d1_text,'string')));
set(handles.d1_text,'string', d1)

d2 = abs(str2double(get(handles.d2_text,'string')));
if d2 >= d1
    d2 = d1/2;
end
set(handles.d2_text,'string', d2)

dp = abs(str2double(get(handles.dp_text,'string')));
set(handles.dp_text,'string', dp)

ply = abs(floor(str2double(get(handles.ply_text,'string'))));
set(handles.ply_text,'string', ply)

angle_c_deg = str2double(get(handles.angle_c_deg_text,'string'));
if angle_c_deg >= 90
    angle_c_deg = 89.9999;
    set(handles.angle_c_deg_text,'string',90)
end

angle_s_deg_vertical = str2double(get(handles.angle_s_deg_vertical_text,'string'));
if angle_s_deg_vertical > 90 
    angle_s_deg_vertical = 90;
elseif angle_s_deg_vertical < 0
    angle_s_deg_vertical = 0;
end
set(handles.angle_s_deg_vertical_text,'string',angle_s_deg_vertical)

c_factor = abs(str2double(get(handles.c_factor_text,'string')));
set(handles.c_factor_text,'string', c_factor)

SiC_density = abs(str2double(get(handles.SiC_text,'string')));
set(handles.SiC_text,'string',SiC_density);

popup_array = get(handles.geometry_popup,'string');
index = get(handles.geometry_popup,'Value');
if strcmp(popup_array{index},'Spherical') == true
    nose_geometry = 's';
elseif strcmp(popup_array{index},'Conical') == true
    nose_geometry = 'c';
else
    nose_geometry = 'X';
end


%% Run function call
[L, A, r_in, r_out, d_frac, angle_c_rad, angle_s_rad] = FTPS_Mass(n, d1, d2, dp, ply, angle_c_deg, angle_s_deg_vertical, c_factor, nose_geometry, SiC_density);

%% Post processing
A_outer = sum(A(1:5))/36^2;
A_inner = (A(1)+A(2)+A(6))/36^2;
A_total = (A_outer+(ply-1)*A_inner);
A_corrected = c_factor*A_total;
SiC_mass = A_corrected*SiC_density;

%% Display outputs in GUI
set(handles.A_outer_output, 'string',A_outer)
set(handles.A_inner_output, 'string',A_inner)
set(handles.A_total_output, 'string',A_total)
set(handles.A_corrected_output, 'string',A_corrected)
set(handles.SiC_output, 'string',SiC_mass)

%% Plot cross section in GUI
hold on;
axis equal;
grid on;
set(gca,'DefaultLineLineWidth',ply,'color',1/255*[50, 53, 56],'GridAlpha',0.2,'GridColor',1/255*[242, 242, 242]);
xlim([-r_out(6)-10, r_out(6)+10]);

% Plots front conical cross section
P1 = [dp/2,0];
P2 = [dp/2+(L(1)+L(5))*sin(angle_c_rad),(L(1)+L(5))*cos(angle_c_rad)];
plot([P1(1), P2(1)], [P1(2), P2(2)],'color',1/255*[222, 157, 27])
plot(-[P1(1), P2(1)], [P1(2), P2(2)],'color',1/255*[222, 157, 27])

% Plots outer shoulder
y_range = [P2(2):0.1:(L(1)+L(5))*cos(angle_c_rad)+d2/2*sin(angle_c_rad)+d2/2*cos(2*asin(d_frac)-pi/2+angle_c_rad)];
x_offset = dp/2+(L(1)+L(5))*sin(angle_c_rad)-d2/2*cos(angle_c_rad);
y_offset = (L(1)+L(5))*cos(angle_c_rad)+d2/2*sin(angle_c_rad);
x_function = sqrt((d2/2)^2-(y_range-y_offset).^2)+x_offset;
x_reflect = -x_function;
plot(x_function, y_range,'color',1/255*[222, 157, 27])
plot(x_reflect, y_range,'color',1/255*[222, 157, 27])

% Plots inner shoulder
x_offset = dp/2+(L(5)+(n-1.5)*d1)*sin(angle_c_rad)-d1/2*cos(angle_c_rad);
y_offset =(L(5)+(n-1.5)*d1)*cos(angle_c_rad)+d1/2*sin(angle_c_rad);
x_range = [x_offset-d1/2*sin(pi/2-angle_c_rad+angle_s_rad):0.1:r_in(4)];
x_reflect = -x_range;
y_function = sqrt((d1/2)^2-(x_range-x_offset).^2)+y_offset;
plot(x_range, y_function,'color',1/255*[222, 157, 27])
plot(x_reflect, y_function,'color',1/255*[222, 157, 27])

% Plots flat between shoulders
P3 = [r_out(4), max(y_range)];
P4 = [r_in(4), y_function(end)];
plot([P3(1), P4(1)], [P3(2), P4(2)],'color',1/255*[222, 157, 27])
plot(-[P3(1), P4(1)], [P3(2), P4(2)],'color',1/255*[222, 157, 27])

% Plot nose geometry
if nose_geometry == 's'
    x_range = [0:0.1:dp/2];
    x_reflect = -x_range;
    y_offset = dp/2*tan(angle_c_rad);
    r = dp/(2*cos(angle_c_rad));
    y_function = -sqrt(r^2-x_range.^2)+y_offset;
    plot(x_range, y_function,'color',1/255*[255, 255, 0])
    plot(x_reflect, y_function,'color',1/255*[255, 255, 0])
    % Plot packing enclosure
    x = [dp/2, dp/2, -dp/2, -dp/2];
    y = [0, dp, dp, 0];
    plot(x,y,'--w','linewidth',1)
    hold off;
elseif nose_geometry == 'c'
    P5 = [0,-dp/2*tan(pi/2-angle_c_rad)];
    plot([P5(1),P1(1)],[P5(2),P1(2)],'color',1/255*[255, 255, 0])
    plot(-[P5(1),P1(1)],[P5(2),P1(2)],'color',1/255*[255, 255, 0])
    % Plot packing enclosure
    x = [dp/2, dp/2, -dp/2, -dp/2];
    y = [0, dp, dp, 0];
    plot(x,y,'--w','linewidth',1)
    hold off;
end

% Plot tori
for i = 1:(n-1)
    hold on;
    set(gca, 'linewidth', 1)
    x_cent = dp/2+(L(5)+(i-0.5)*d1)*sin(angle_c_rad)-d1/2*cos(angle_c_rad);
    y_cent = (L(5)+(i-0.5)*d1)*cos(angle_c_rad)+d1/2*sin(angle_c_rad);
    r = d1/2;
    circle(x_cent, y_cent, r);
    circle(-x_cent, y_cent, r);
end

% Plot shoulder torus
x_cent = dp/2+(L(1)+L(5))*sin(angle_c_rad)-d2/2*cos(angle_c_rad);
y_cent = (L(1)+L(5))*cos(angle_c_rad)+d2/2*sin(angle_c_rad);
r = d2/2;
circle(x_cent, y_cent, r);
circle(-x_cent, y_cent, r);




function circle(x,y,r)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
ang=0:0.1:2*pi;
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp,'--w','LineWidth',1);
