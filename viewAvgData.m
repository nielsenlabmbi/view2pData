function varargout = viewAvgData(varargin)
%%%%% this works with the realtime files (i.e. cells selected during
%%%%% acquisition)


% VIEWAVGDATA MATLAB code for viewAvgData.fig
%      VIEWAVGDATA, by itself, creates a new VIEWAVGDATA or raises the existing
%      singleton*.
%
%      H = VIEWAVGDATA returns the handle to a new VIEWAVGDATA or the handle to
%      the existing singleton*.
%
%      VIEWAVGDATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEWAVGDATA.M with the given input arguments.
%
%      VIEWAVGDATA('Property','Value',...) creates a new VIEWAVGDATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before viewAvgData_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to viewAvgData_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help viewAvgData

% Last Modified by GUIDE v2.5 11-May-2015 18:28:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @viewAvgData_OpeningFcn, ...
                   'gui_OutputFcn',  @viewAvgData_OutputFcn, ...
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


% --- Executes just before viewAvgData is made visible.
function viewAvgData_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to viewAvgData (see VARARGIN)

%read file open icon
foIcon=imread('fileopen.jpg');
foIcon=imresize(foIcon,[20 20]);
foIcon=double(foIcon)./255;

set(handles.b2Pfilep,'CData',foIcon);
set(handles.bAnafilep,'CData',foIcon);


% Choose default command line output for viewAvgData
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes viewAvgData wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = viewAvgData_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%% FILE DETAILS PANEL

% --- Executes during object creation, after setting all properties.
function eAnimalID_CreateFcn(hObject, eventdata, handles) %#ok<DEFNU,*INUSD>
% hObject    handle to eAnimalID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function eUnitID_CreateFcn(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to eUnitID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function eExpID_CreateFcn(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to eExpID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function e2Pfilep_CreateFcn(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to e2Pfilep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function eAnafilep_CreateFcn(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to eAnafilep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in b2Pfilep.
function b2Pfilep_Callback(hObject, eventdata, handles)  %#ok<DEFNU>
% hObject    handle to b2Pfilep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fn=uigetdir;
if ~isequal(fn,0)
    set(handles.e2Pfilep,'String',fn);
end


% --- Executes on button press in bAnafilep.
function bAnafilep_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to bAnafilep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fn=uigetdir;
if ~isequal(fn,0)
    set(handles.eAnafilep,'String',fn);
end

% --- Executes on button press in bLoad.
function bLoad_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to bLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



%get basic components for file names
animalID=get(handles.eAnimalID,'String');
unitID=get(handles.eUnitID,'String');
expID=get(handles.eExpID,'String');
base2P=get(handles.e2Pfilep,'String');
baseAna=get(handles.eAnafilep,'String');

%make file names
handles.filename2P=fullfile(base2P,animalID,[animalID '_' unitID '_' expID]);
filenameAna=fullfile(baseAna,animalID,[animalID '_u' unitID '_' expID '.analyzer']);



%disp(filename2P)

%check whether both files exist and throw error message if not
exist2P=exist([handles.filename2P '.mat'],'file');
existAna=exist(filenameAna,'file');

if exist2P~=2
    errordlg('2P file not found!');
end
if existAna~=2
    errordlg('Analyzer file not found!');
end

%proceed if both files exist
if exist2P==2 && existAna==2
    %load 2P file - this generates roipix, ttl_log, rtdata
    load([handles.filename2P '_realtime.mat']);
    
    %also load info file for timing info
    load(handles.filename2P); 
    
    %get times per frame
    handles.tPerFrame=info.config.lines/info.resfreq;
    
    %get lines per frame
    handles.nrLines=info.config.lines;
    handles.nrPixels=796;
    
    %load analyzer file - this generates Analyzer
    load(filenameAna,'-mat');
    handles.Analyzer=Analyzer;
    
    %from Analyzer file, extract list of trials and conditions, as well as
    %nr of trials
    handles.triallist=getcondtrial(Analyzer);
    [handles.domains,handles.domval,handles.blankid]=getdomainvalue(Analyzer);
    handles.Ntrial=getnotrials(Analyzer);
    handles.Nrepeat=getnorepeats(1,Analyzer); %1 is always one of the stimulus cond, not the blank
    if ~isempty(handles.blankid)
        handles.NrepeatBlank=getnorepeats(handles.blankid,Analyzer);
    else
        handles.NrepeatBlank = NaN;
    end

    %extract stimulus on events
    handles.stimOn=find(diff(double(ttl_log))>0)+1;
    
    %determine nr of cells and get cell ROIs
    handles.Ncell=size(rtdata,2);
    handles.roipix=roipix;
    
    %transform rtdata
    handles.rtdata=double(intmax('uint16'))-rtdata;
    %handles.rtdata=rtdata;
    
    %set up cell selection part
    handles.inclCell=ones(handles.Ncell,1);
    
    %data loaded, update time windows in GUI and enable next button
    pv=getparam('predelay',Analyzer);
    set(handles.tBase,'String',['max: ' num2str(pv*1000) ' ms']);
    pv=getparam('stim_time',Analyzer);
    set(handles.tStim,'String',['max: ' num2str(pv*1000) ' ms']);
    
    set(handles.bTimeWindow,'Enable','on');
    
    %also reset all other selections in case the GUI was used before
    set(handles.sCond1,'String',' ','Value',1);
    set(handles.sCond2,'String',' ','Value',1);
    set(handles.sCirc,'Value',0);
    set(handles.sMod,'Value',0);
    set(handles.sAvg,'Value',1);
    set(handles.sAll,'Value',0);
    set(handles.sVar,'Value',0);
    set(handles.sValVar2,'String',' ','Value',1);
    set(handles.bPlotFrames,'Enable','off');
    set(handles.bPlotAvg,'Enable','off');
    set(handles.bRet,'Enable','off');
    cla(handles.axes1);
    cla(handles.axes2);
    cla(handles.axes3);
    
    
    %notify that done with import
    helpdlg('Finished loading data!')
end %if
guidata(hObject, handles);


%%%%%GET TIME COURSES PART OF GUI

function eBaseLow_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to eBaseLow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eBaseLow as text
%        str2double(get(hObject,'String')) returns contents of eBaseLow as a double


% --- Executes during object creation, after setting all properties.
function eBaseLow_CreateFcn(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to eBaseLow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eStimLow_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to eStimLow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eStimLow as text
%        str2double(get(hObject,'String')) returns contents of eStimLow as a double


% --- Executes during object creation, after setting all properties.
function eStimLow_CreateFcn(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to eStimLow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function eRange1_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to eRange1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eRange1 as text
%        str2double(get(hObject,'String')) returns contents of eRange1 as a double


% --- Executes during object creation, after setting all properties.
function eRange1_CreateFcn(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to eRange1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function eRange2_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to eRange2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eRange2 as text
%        str2double(get(hObject,'String')) returns contents of eRange2 as a double


% --- Executes during object creation, after setting all properties.
function eRange2_CreateFcn(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to eRange2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in bTimeWindow.
function bTimeWindow_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to bTimeWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get duration of intervals in frames
baseSec=str2num(get(handles.eBaseLow,'String')); %#ok<*ST2NM>
baseFrameLow=round((baseSec/1000)/handles.tPerFrame);
    
baseSec=str2num(get(handles.eBaseHigh,'String'));
baseFrameHigh=round((baseSec/1000)/handles.tPerFrame);

stimSec=str2num(get(handles.eStimLow,'String'));
stimFrameLow=round((stimSec/1000)/handles.tPerFrame);
  
stimSec=str2num(get(handles.eStimHigh,'String'));
stimFrameHigh=round((stimSec/1000)/handles.tPerFrame);

tcSec=str2num(get(handles.eRange1,'String'));
tcFrame1=round((tcSec/1000)/handles.tPerFrame);

tcSec=str2num(get(handles.eRange2,'String'));
tcFrame2=round((tcSec/1000)/handles.tPerFrame);

%get response data for every trial - this computes dF/F for every trial
handles.trialResp=[];
handles.trialAvg=[];
for i=1:handles.Ntrial
    %get F0 first
    F0=mean(handles.rtdata(handles.stimOn(i)-abs(baseFrameLow):handles.stimOn(i)-abs(baseFrameHigh),:),1);
    
    %now compute dF/F for every cell
    for c=1:handles.Ncell
        %get time course - this is determined by the selected range
        tmpCell=(handles.rtdata(handles.stimOn(i)+tcFrame1:handles.stimOn(i)+tcFrame2,c)-F0(c))/F0(c);
        handles.trialResp{i}(:,c)=tmpCell;
        
        %also compute average response in specified response period
        handles.trialAvg(i,c)=mean(tmpCell(stimFrameLow:stimFrameHigh));
    end
end    
handles.NTimepoints=abs(tcFrame1)+tcFrame2+1;
handles.timevec=([1:handles.NTimepoints]+tcFrame1-1)*handles.tPerFrame;

%now that data has been loaded, set the plot part of the GUI up
%correctly
    
%populate selection for conditions
set(handles.sCond1,'String',handles.domains);
if length(handles.domains)==1
    set(handles.sCond2,'Enable','off');
    set(handles.sAvg,'Enable','off');
    set(handles.sAll,'Enable','off');
    set(handles.sVar,'Enable','off');
    set(handles.sValVar2,'Enable','off');
else
    %enable part just to make sure that switching between files doesn't
    %generate problems
    set(handles.sCond2,'Enable','on');
    set(handles.sAvg,'Enable','on');
    set(handles.sAll,'Enable','on');
    set(handles.sCond2,'String',handles.domains,'Value',2);
    set(handles.sVar,'Enable','on');
    set(handles.sValVar2,'Enable','on');
    set(handles.sValVar2,'String',num2str(unique(handles.domval(:,2))));
end
    
%enable plotting
set(handles.bPlotAvg,'Enable','on');
set(handles.bPlotFrames,'Enable','on');

if all(strcmp(handles.domains,{'a' 'b'}))
    set(handles.bRet,'Enable','on');
end

%notify that done with import
helpdlg('Finished loading data!')

guidata(hObject, handles);



%%%%% PLOT PART OF THE GUI

% --- Executes during object creation, after setting all properties.
function sCond1_CreateFcn(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to sCond1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in sCond1.
function sCond1_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to sCond1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sCond1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sCond1

if length(handles.domains)==2 %coupling two parameters together like doesn't make sense if there are 3 or more
    varid=get(hObject,'Value'); %is either 1 or 2
    set(handles.sCond2,'Value',3-varid);
    set(handles.sValVar2,'String',num2str(unique(handles.domval(:,3-varid))));
end

% --- Executes during object creation, after setting all properties.
function sCond2_CreateFcn(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to sCond2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in sCond2.
function sCond2_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to sCond2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sCond2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sCond2

varid=get(hObject,'Value');

if length(handles.domains)==2
    set(handles.sCond1,'Value',3-varid);
end

set(handles.sValVar2,'String',num2str(unique(handles.domval(:,varid))));

% --- Executes during object creation, after setting all properties.
function eMod_CreateFcn(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to eMod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in sAvg.
function sAvg_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to sAvg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sAvg
if get(hObject,'Value')==1
    set(handles.sAll,'Value',0);
    set(handles.sVar,'Value',0);
else
    set(handles.sAll,'Value',1);
    set(handles.sVar,'Value',0);
end



% --- Executes on button press in sAll.
function sAll_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to sAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sAll
if get(hObject,'Value')==1
    set(handles.sAvg,'Value',0);
    set(handles.sVar,'Value',0);
else
    set(handles.sAvg,'Value',0);
    set(handles.sVar,'Value',1);
end


% --- Executes on button press in sVar.
function sVar_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to sVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sVar
if get(hObject,'Value')==1
    set(handles.sAvg,'Value',0);
    set(handles.sAll,'Value',0);
else
    set(handles.sAvg,'Value',1);
    set(handles.sAll,'Value',0);
end



% --- Executes during object creation, after setting all properties.
function sValVar2_CreateFcn(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to sValVar2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%%% PLOT INDIVIDUAL FRAMES

% --- Executes on button press in bPlotFrames.
function bPlotFrames_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to bPlotFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%show individual frames
global info

%first open figure
if ~isfield(handles,'plotFHandle') || ~ishandle(handles.plotFHandle)  %make sure figure isn't already open
    handles.plotFHandle=figure;
end

handles.plotFrame=0;

tmp=sbxread(handles.filename2P,handles.plotFrame,1);

if info.nchan==1
    imshow(squeeze(tmp(1,:,:)),[]);
    title(['Frame: ' num2str(handles.plotFrame)]);
else
    for i=1:2
        subplot(2,1,i)
        imshow(squeeze(tmp(i,:,:)),[]);
        if i==1
            title(['Frame: ' num2str(handles.plotFrame)]);
        end
    end
end


set(handles.slFrames,'Enable','on');
set(handles.slFrames,'Min',0);
set(handles.slFrames,'Max',info.max_idx);
set(handles.slFrames,'SliderStep',[1/info.max_idx 1])
set(handles.slFrames,'Value',1);

guidata(hObject, handles);



% --- Executes on slider movement.
function slFrames_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to slFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global info;

handles.plotFrame=round(get(hObject,'Value'));

figure(handles.plotFHandle);

tmp=sbxread(handles.filename2P,handles.plotFrame,1);

if info.nchan==1
    imshow(squeeze(tmp(1,:,:)),[]);
    title(['Frame: ' num2str(handles.plotFrame)]);
else
    for i=1:2
        subplot(2,1,i)
        imshow(squeeze(tmp(i,:,:)),[]);
        if i==1
            title(['Frame: ' num2str(handles.plotFrame)]);
        end
    end
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slFrames_CreateFcn(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to slFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




%%% PLOT AVERAGE

% --- Executes on button press in bPlotAvg.
function bPlotAvg_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to bPlotAvg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get primary condition and its values
primCondID=get(handles.sCond1,'Value'); %this tells you which domain (i.e. which looper parameter) was selected
%get the values for that parameter, either directly, or modulo whatever
%value was set
if get(handles.sMod,'Value')==0
    handles.primCondVal=unique(handles.domval(:,primCondID)); 
else
    handles.modVal=str2num(get(handles.eMod,'String'));
    handles.primCondVal=unique(mod(handles.domval(:,primCondID),handles.modVal));
end


%get secondary condition and its values
secCondID=get(handles.sCond2,'Value');
handles.secCondVal=unique(handles.domval(:,secCondID));

%disp(secCondID)

%now compute averages across trials
handles.avgTimeCond=[];
handles.sDevTimeCond=[];
handles.avgCond=[];
handles.sDevCond=[];
    
handles.nCond1=length(handles.primCondVal);
handles.nCond2=length(handles.secCondVal);
    
if get(handles.sAvg,'Value')==1 %collapse across condition 2
    for c=1:handles.nCond1
                
        %get list of conditions that have the correct parameters
        if get(handles.sMod,'Value')==0 
            condID=find(handles.domval(:,primCondID)==handles.primCondVal(c)); 
        else
            condID=find(mod(handles.domval(:,primCondID),handles.modVal)==handles.primCondVal(c));
        end
            
        trialID=find(ismember(handles.triallist,condID)); %find the matching trials
        %necessary because if there are multiple parameters, the same
        %parameter value shows up multiple times in the condition list
        
        %compute average timecourse
        tmpAvg=zeros(handles.NTimepoints,handles.Ncell);
        tmpSDev=zeros(handles.NTimepoints,handles.Ncell);
        for i=1:length(trialID)
            tmpAvg=tmpAvg+handles.trialResp{trialID(i)}/length(trialID);
            stemp=tmpSDev+handles.trialResp{trialID(i)}/length(trialID).^2/length(trialID); %for standard deviation
        end
        handles.avgTimeCond(:,:,c)=tmpAvg;
        handles.sDevTimeCond(:,:,c)=sqrt(stemp-tmpAvg.^2);
        
        handles.avgCond(c,:)=mean(handles.trialAvg(trialID,:),1);
        handles.sDevCond(c,:)=std(handles.trialAvg(trialID,:),0,1);
    end %for cond
    
    %get max and min across conditions
    [maxResp,handles.maxCond]=max(handles.avgCond,[],1);
    [minResp,handles.minCond]=min(handles.avgCond,[],1);
    
    %get individual timecourses for the best  condition
    for c=1:handles.Ncell
       bc=handles.maxCond(c); %cond id for best condition along dim 1
       
       if get(handles.sMod,'Value')==0
           condID=find(handles.domval(:,primCondID)==handles.primCondVal(bc));
       else
           condID=find(mod(handles.domval(:,primCondID),handles.modVal)==handles.primCondVal(bc));
       end
            
       trialID=find(ismember(handles.triallist,condID)); %find the matching trials
       handles.bestResp{c}=[];
       for i=1:length(trialID)
           handles.bestResp{c}(:,i)= handles.trialResp{trialID(i)}(:,c);
       end   
    end
    
    
elseif get(handles.sAll,'Value')==1 %keep all conditions for 2nd parameter
    
    for c1=1:handles.nCond1
        for c2=1:handles.nCond2  
        
            %get list of conditions that have the correct parameters
            if get(handles.sMod,'Value')==0
                condID=find(handles.domval(:,primCondID)==handles.primCondVal(c1) & ...
                    handles.domval(:,secCondID)==handles.secCondVal(c2));
            else
                condID=find(mod(handles.domval(:,primCondID),handles.modVal)==handles.primCondVal(c1) & ...
                     handles.domval(:,secCondID)==handles.secCondVal(c2));
            end
            
            trialID=find(ismember(handles.triallist,condID)); %find the matching trials
    
            
            tmpAvg=zeros(handles.NTimepoints,handles.Ncell);
            tmpSDev=zeros(handles.NTimepoints,handles.Ncell);
            for i=1:length(trialID)
                tmpAvg=tmpAvg+handles.trialResp{trialID(i)}/length(trialID);
                stemp=tmpSDev+handles.trialResp{trialID(i)}/length(trialID).^2/length(trialID); %for standard deviation
            end
            handles.avgTimeCond(:,:,c1,c2)=tmpAvg;
            handles.sDevTimeCond(:,:,c1,c2)=sqrt(stemp-tmpAvg.^2);
        
            handles.avgCond(c1,c2,:)=mean(handles.trialAvg(trialID,:),1);
            handles.sDevCond(c1,c2,:)=std(handles.trialAvg(trialID,:),0,1);
            
        end
    end
    
    
    %get max and min across conditions 
    [mTemp,maxCond2]=max(handles.avgCond,[],2);
    mTemp=squeeze(mTemp);
    handles.maxCond2=squeeze(maxCond2);
    [maxResp,handles.maxCond]=max(mTemp,[],1);
        
    [mTemp,minCond2]=min(handles.avgCond,[],2);
    mTemp=squeeze(mTemp);
    handles.minCond2=squeeze(minCond2);
    [minResp,handles.minCond]=min(mTemp,[],1);
           
    %get individual timecourses for the best  condition
    for c=1:handles.Ncell
        bc1=handles.maxCond(c); %cond id for best condition along dim 1
        bc2=handles.maxCond2(bc1,c); %cond id for best condition along dim 2
        
        if get(handles.sMod,'Value')==0
            condID=find(handles.domval(:,primCondID)==handles.primCondVal(bc1) & ...
                handles.domval(:,secCondID)==handles.secCondVal(bc2));
        else
            condID=find(mod(handles.domval(:,primCondID),handles.modVal)==handles.primCondVal(bc1) & ...
                handles.domval(:,secCondID)==handles.secCondVal(bc2));
        end
        
        trialID=find(ismember(handles.triallist,condID)); %find the matching trials
        handles.bestResp{c}=[];
        for i=1:length(trialID)
            handles.bestResp{c}(:,i)= handles.trialResp{trialID(i)}(:,c);
        end
    end
    
elseif get(handles.sVar,'Value')==1 %select one particular condition for the 2nd parameter
    
    c2=get(handles.sValVar2,'Value');
   
    for c=1:handles.nCond1
                
        %get list of conditions that have the correct parameters
        if get(handles.sMod,'Value')==0
            condID=find(handles.domval(:,primCondID)==handles.primCondVal(c) & ...
                    handles.domval(:,secCondID)==handles.secCondVal(c2)); 
        else
            condID=find(mod(handles.domval(:,primCondID),handles.modVal)==handles.primCondVal(c) & ...
                    handles.domval(:,secCondID)==handles.secCondVal(c2));
        end
            
        trialID=find(ismember(handles.triallist,condID)); %find the matching trials
        %necessary because if there are multiple parameters, the same
        %parameter value shows up multiple times in the condition list
        
        %compute average timecourse
        tmpAvg=zeros(handles.NTimepoints,handles.Ncell);
        tmpSDev=zeros(handles.NTimepoints,handles.Ncell);
        for i=1:length(trialID)
            tmpAvg=tmpAvg+handles.trialResp{trialID(i)}/length(trialID);
            stemp=tmpSDev+handles.trialResp{trialID(i)}/length(trialID).^2/length(trialID); %for standard deviation
        end
        handles.avgTimeCond(:,:,c)=tmpAvg;
        handles.sDevTimeCond(:,:,c)=sqrt(stemp-tmpAvg.^2);
        
        handles.avgCond(c,:)=mean(handles.trialAvg(trialID,:),1);
        handles.sDevCond(c,:)=std(handles.trialAvg(trialID,:),0,1);
    end %for cond
    
    %get max and min across conditions
    [maxResp,handles.maxCond]=max(handles.avgCond,[],1);
    [minResp,handles.minCond]=min(handles.avgCond,[],1);
    
    %get individual timecourses for the best  condition
    for c=1:handles.Ncell
       bc=handles.maxCond(c); %cond id for best condition along dim 1
       
       if get(handles.sMod,'Value')==0
           condID=find(handles.domval(:,primCondID)==handles.primCondVal(bc) & ...
                    handles.domval(:,secCondID)==handles.secCondVal(c2));
       else
           condID=find(mod(handles.domval(:,primCondID),handles.modVal)==handles.primCondVal(bc) & ...
                    handles.domval(:,secCondID)==handles.secCondVal(c2));
       end
            
       trialID=find(ismember(handles.triallist,condID)); %find the matching trials
       handles.bestResp{c}=[];
       for i=1:length(trialID)
           handles.bestResp{c}(:,i)= handles.trialResp{trialID(i)}(:,c);
       end   
    end
    
end % if average


%compute response strength - difference between best and worst
%condition
mag=maxResp-minResp;

 
%build image
maxCondImg=zeros(handles.nrLines,handles.nrPixels);
cellMask=zeros(handles.nrLines,handles.nrPixels);
cellMag=zeros(handles.nrLines,handles.nrPixels);
for i=1:handles.Ncell
    maxCondImg(handles.roipix{i})=handles.maxCond(i);
    cellMask(handles.roipix{i})=1;
    cellMag(handles.roipix{i})=mag(i);
end

cellMag = cellMag-min(cellMag(:));
cellMag = cellMag/max(cellMag(:));


%convert into color
if get(handles.sCirc,'Value') || get(handles.sMod,'Value')
    handles.nColors=handles.nCond1+1;
else
    handles.nColors=handles.nCond1;
end
maxCondImg=(maxCondImg-1)/(handles.nColors-1); %scale between 0 and 1, making sure that the first condition repeats
maxCondImg=round(maxCondImg*63+1);


%determine which plot options were selected
scaleMag=get(handles.sScale,'Value');
anatFlag=get(handles.sAnatomy,'Value');

if anatFlag==1 %compute anatomy picture 
    imAnat=sbxread(handles.filename2P,0,10);    
    imAnat=squeeze(mean(imAnat(1,:,:,:),4));
    
    imAnat=imAnat-min(imAnat(:));
    imAnat=imAnat/max(imAnat(:));
    
    imAnat(:,:,2) = imAnat;
    imAnat(:,:,3) = imAnat(:,:,1);
end
%disp(size(imAnat))

% main image with cells
if get(handles.sCirc,'Value') && get(handles.sMod,'Value')==0
    ytlabel=[handles.primCondVal;handles.primCondVal(1)];
elseif get(handles.sMod,'Value')
    ytlabel=[handles.primCondVal;handles.modVal];
else
    ytlabel=handles.primCondVal;
end

if anatFlag==0 
    axes(handles.axes1)
    if scaleMag==0
        image(1:length(maxCondImg(1,:)),1:length(maxCondImg(:,1)),maxCondImg,...
            'CDataMapping','direct','AlphaData',cellMask,'AlphaDataMapping','none');
    else
        image(1:length(maxCondImg(1,:)),1:length(maxCondImg(:,1)),maxCondImg,...
        'CDataMapping','direct','AlphaData',cellMag,'AlphaDataMapping','none');
    end
    axis image;
    set(gca,'XTick',[],'YTick',[])
    if ~get(handles.sCirc,'Value')
        colormap('jet');
    else
        colormap('hsv');
    end
    cbar=colorbar;
    yt=linspace(1,64,handles.nColors);
    set(cbar,'YTick',yt,'YTicklabel',ytlabel)
else 
    if ~get(handles.sCirc,'Value')
        cid = jet;
    else
        cid = hsv;
    end
    
    imout = zeros([size(maxCondImg) 3]);
    for i = 1:size(maxCondImg,1)
        for j = 1:size(maxCondImg,2)
            if cellMask(i,j)==1
                if scaleMag==1
                    imout(i,j,:) = cellMag(i,j)*cid(maxCondImg(i,j),:);
                else
                    imout(i,j,:) = cid(maxCondImg(i,j),:);
                end
            else
                imout(i,j,:)=[0 0 0];
            end
        end
    end
    
    imout = imout + imAnat;
    imout = imout/max(imout(:));
    axes(handles.axes1)
    image(imout,'CDataMapping','direct','AlphaDataMapping','none');
    set(gca,'XTick',[],'YTick',[])
    cbar=colorbar;
    yt=linspace(1,64,handles.nColors);
    set(cbar,'YTick',yt,'YTicklabel',ytlabel)
end


%add marker for first cell
locxdum = ceil(handles.roipix{1}/handles.nrLines);
locydum = handles.roipix{1}-handles.nrLines*(locxdum-1);
cellCenterX = round(mean(locxdum));
cellCenterY = round(mean(locydum));
handles.textMarker=text(cellCenterX,cellCenterY,'+',...
    'Color',[1 0.78 0.80],'HorizontalAlignment','center','FontWeight','bold');


%initialize selected cell to 1 and plot cell data
handles.cellSelected = 1;
set(handles.tCell,'String',num2str(handles.cellSelected));
set(handles.rInclude,'Value',handles.inclCell(handles.cellSelected));
plotData(handles)

guidata(hObject, handles);






%%%%% PLOT 2D RETINOTOPY
% --- Executes on button press in bRet.
function bRet_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to bRet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%this assumes that there are two conditions, named 'a' and 'b'
%a hast the steps, b distinguishes between 0 and 1
%conditions: b=0: vertical stripes -> measures along x
% b=1: horizontal stripes -> measures along y

%error checking:
if all(strcmp(handles.domains,{'a' 'b'}))==0
    errordlg('Incorrect analyzer structure for retinotopy');
else
    %determine numbers of bins
    primCondVal=unique(handles.domval(:,1));
    secCondVal=unique(handles.domval(:,2));
    nCond1=length(primCondVal);
    nCond2=length(secCondVal);
    
    if nCond2~=2
        errordlg('Incorrect analyzer structure for retinotopy');
    end
    
    
    %we need to get the exact settings for the stimuli, which are in the
    %formula in the analyzer file
    f=handles.Analyzer.loops.formula;
    f2=strsplit(f,';');
    for i=1:length(f2)
        tmp=strsplit(f2{i},'=');
        part1{i}=tmp{1}(~isspace(tmp{1}));
        part2{i}=tmp{2}(~isspace(tmp{2}));
    end
    
    %find certain parameters in the formula
    x0=0;
    y0=0;
    for i=1:length(part1)
        if strfind(part1{i},'N')
            N=str2num(part2{i});
        end
        if strfind(part1{i},'xW')
            xW=str2num(part2{i});
        end
        if strfind(part1{i},'yW')
            yW=str2num(part2{i});
        end
        if strfind(part1{i},'x0')
            x0=str2num(part2{i});
        end
        if strfind(part1{i},'y0')
            y0=str2num(part2{i});
        end
    end
    
    %also get screen size from parameters
    pcmX=handles.Analyzer.M.xpixels/handles.Analyzer.M.screenXcm;
    pcmY=handles.Analyzer.M.ypixels/handles.Analyzer.M.screenYcm;
    x_size=360*(xW/N)/(2*pi*handles.Analyzer.M.screenDist*pcmX);
    y_size=360*(yW/N)/(2*pi*handles.Analyzer.M.screenDist*pcmY);
    
    %get average response per condition
    avgCond=[];
    for c1=1:nCond1 %this is a
        for c2=1:2 % this is b
        
            %get list of conditions that have the correct parameters
            condID=find(handles.domval(:,1)==primCondVal(c1) & ...
                    handles.domval(:,2)==secCondVal(c2));
         
            trialID=find(ismember(handles.triallist,condID)); %find the matching trials
            avgCond(c1,c2,:)=mean(handles.trialAvg(trialID,:),1);        
        end
    end
    
    %get maximum bins for horizontal (b=1) and vertical (b=0) bars
    [~,maxB0]=max(squeeze(avgCond(:,1,:)),[],1);
    [~,maxB1]=max(squeeze(avgCond(:,2,:)),[],1);
    
    %also compute CoM for each tuning function
    %for this part, set all negative values to 0
    recCond=avgCond;
    recCond(recCond<0)=0;
    
    stepvec=linspace(0,180,nCond1);
    for i=1:handles.Ncell
        for c2=1:2
            %set sum of tuning function for each neuron to 1
            sumR=sum(squeeze(recCond(:,c2,i)));
            normCond=squeeze(recCond(:,c2,i))./sumR;
            
            %compute angle of vector sum
            CoM(i,c2)=angle(sum(normCond'.*exp(1i*stepvec*pi/180)));
        end
    end
    CoM=CoM/pi*(nCond1-1); %this is back in bin space 
    CoMPixX=CoM(:,1)*xW/N+xW/(2*N)+x0;
    CoMPixY=CoM(:,2)*yW/N+yW/(2*N)+y0;
    
    
    %plot 1: surface plot showing nr cells that have maximum response for a particular bin
    retbins=zeros(nCond1);
    for i=1:handles.Ncell
        retbins(maxB1(i),maxB0(i))=retbins(maxB1(i),maxB0(i))+1;
    end
    
    %scale for plotting
    maxCells=max(retbins(:));
    retbins=retbins/maxCells;
    retbins=round(retbins*64);
    
    figure
    
    %coordinates for plotting
    x=[0:4]*xW/N+xW/(2*N)+x0;
    y=[0:4]*yW/N+yW/(2*N)+y0;    
    xedge=[1:4]*xW/N+x0;
    yedge=[1:4]*yW/N+y0;
    xtick=[x0 xedge x0+xW];
    ytick=[y0 yedge y0+yW];
    
    %plot distribution of max
    subplot(1,2,1)
    image(x,y,retbins,'CDataMapping','direct');
    hold on
    for i=1:4
        plot([xedge(i) xedge(i)],[y0 y0+yW],'w')
        plot([x0 x0+xW],[yedge(i) yedge(i)],'w')
    end
    hold off
    colormap('jet')
    set(gca,'XTick',xtick,'YTick',ytick);
    
    h=colorbar;
    ytlabel=linspace(1,maxCells,5);
    yt=linspace(1,64,5);
    set(h,'YTick',yt,'YTickLabel',ytlabel);
    h.Label.String='N Neurons';
    title('Maximum responses')
    
    %plot distribution of CoM
    subplot(1,2,2)
    plot(CoMPixX,CoMPixY,'k+','LineWidth',2);
    hold on
    for i=1:4
        plot([xedge(i) xedge(i)],[y0 y0+yW],'k--')
        plot([x0 x0+xW],[yedge(i) yedge(i)],'k--')
    end
    hold off
    set(gca,'XLim',[x0 x0+xW],'YLim',[y0 y0+yW],'XTick',xtick,'YTick',ytick)
    xlabel(['x, 1 bin = ' sprintf('%0.1f',x_size) 'deg'])
    ylabel(['y, 1 bin = ' sprintf('%0.1f',y_size) 'deg'])
    axis ij
    title('Center of mass')
    
    
end %error checking


%%%% CELL SELECTION
% --- Executes on button press in bDecCell.
function bDecCell_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to bDecCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cellSelected=max(1,handles.cellSelected-1);
set(handles.tCell,'String',num2str(handles.cellSelected));
set(handles.rInclude,'Value',handles.inclCell(handles.cellSelected));
plotData(handles)

guidata(hObject, handles);

% --- Executes on button press in bIncCell.
function bIncCell_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to bIncCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.cellSelected=min(handles.Ncell,handles.cellSelected+1);
set(handles.tCell,'String',num2str(handles.cellSelected));
set(handles.rInclude,'Value',handles.inclCell(handles.cellSelected));
plotData(handles)

guidata(hObject, handles);

% --- Executes on button press in rInclude.
function rInclude_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to rInclude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rInclude


% --- Executes on button press in bUpdate.
function bUpdate_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to bUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
v=get(handles.rInclude,'Value');
set(handles.rInclude,'Value',1-v);
handles.inclCell(handles.cellSelected)=1-v;
guidata(hObject, handles);




%%%%%   SAVE DATA
% --- Executes on button press in bSave.
function bSave_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to bSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ncond = size(handles.domval,1);

resp = zeros(ncond,handles.Ncell,handles.Nrepeat,handles.NTimepoints);
blankResp = zeros(1,handles.Ncell,handles.NrepeatBlank,handles.NTimepoints);
times = zeros(ncond,handles.Ncell,handles.Nrepeat,handles.NTimepoints);
blankTimes = zeros(1,handles.Ncell,handles.NrepeatBlank,handles.NTimepoints);

for i=1:ncond
    idx = find(handles.triallist==i);
    for j=1:length(idx)
        for c=1:handles.Ncell
            resp(i,c,j,:) = handles.trialResp{idx(j)}(:,c);
            times(i,c,j,:) = handles.timevec;
        end
    end
end

idx = find(handles.triallist == handles.blankid);
for j=1:length(idx)
    for c=1:handles.Ncell
        blankResp(1,c,j,:) = handles.trialResp{idx(j)}(:,c);
        blankTimes(i,c,j,:) = handles.timevec;
    end
end

unitStat = ones(1,handles.Ncell); 
respRange = [str2num(get(handles.eStimLow,'String')) str2num(get(handles.eStimHigh,'String'))];

fn = uigetdir;

if ~isequal(fn,0)
    save(fullfile(fn, 'resp.mat'),'resp','blankResp','unitStat','times','blankTimes','respRange');
    msgbox('Responses saved.','Success','help');
end


%%%%% plotting function
function plotData(handles)

%get center of mass for the selected cell
locxdum = ceil(handles.roipix{handles.cellSelected}/handles.nrLines);
locydum = handles.roipix{handles.cellSelected}-handles.nrLines*(locxdum-1);
cellCenterX = round(mean(locxdum));
cellCenterY = round(mean(locydum));
set(handles.textMarker,'Position',[cellCenterX cellCenterY 0]);

%left subplot: tuning curve
xval=handles.primCondVal;
if get(handles.sMod,'Value')
    xval=[handles.primCondVal;handles.modVal];
end        
lval=handles.secCondVal;
primCondID=get(handles.sCond1,'Value');

if get(handles.sAll,'Value')==1 %collapse across condition 2  
    xmat=repmat(xval,1,handles.nCond2);    
    ymat=squeeze(handles.avgCond(:,:,handles.cellSelected));
    ermat=squeeze(handles.sDevCond(:,:,handles.cellSelected));
    if get(handles.sMod,'Value')
        ymat=[ymat;ymat(1,:)];
        ermat=[ermat;ermat(1,:)];
    end
    axes(handles.axes2)
    errorbar(xmat,ymat,ermat,'o-')
    set(gca,'XLim',[xval(1) xval(end)])
    ylabel('DF/F')
    xlabel(handles.domains{primCondID})
    legend(num2str(lval),'location','best')
   
else 
    yval=squeeze(handles.avgCond(:,handles.cellSelected));
    erval=squeeze(handles.sDevCond(:,handles.cellSelected));
    if get(handles.sMod,'Value')
        yval=[yval;yval(1)];
        erval=[erval;erval(1)];
    end
       
    axes(handles.axes2)
    errorbar(xval,yval,erval,'o-r')
    set(gca,'XLim',[xval(1) xval(end)])
    ylabel('DF/F')
    xlabel(handles.domains{primCondID})
end
    
%right subplot: timecourse for best and worst condition
bc1=handles.maxCond(handles.cellSelected); %cond id for best condition along dim 1
mc1=handles.minCond(handles.cellSelected);
txtStr=['C1: ' num2str(xval(bc1))];
if get(handles.sAll,'Value')==1
    bc2=handles.maxCond2(bc1,handles.cellSelected); %cond id for best condition along dim 2
    txtStr=[txtStr ', C2: ' num2str(lval(bc2))];
    mc2=handles.minCond2(mc1,handles.cellSelected); %cond id for best condition along dim 2
    tcbest=squeeze(handles.avgTimeCond(:,handles.cellSelected,bc1,bc2));
    tcworst=squeeze(handles.avgTimeCond(:,handles.cellSelected,mc1,mc2));
else
    tcbest=squeeze(handles.avgTimeCond(:,handles.cellSelected,bc1));
    tcworst=squeeze(handles.avgTimeCond(:,handles.cellSelected,mc1));
end
tctrial=squeeze(handles.bestResp{handles.cellSelected});

set(handles.eBestCond,'String',txtStr);
axes(handles.axes3)
plot(handles.timevec*1000,tcbest,'-r','LineWidth',2)
hold on
plot(handles.timevec*1000,tcworst,'-b','LineWidth',2)
plot(handles.timevec*1000,tctrial,'-r')
hold off
set(gca,'XLim',[handles.timevec(1)*1000 handles.timevec(end)*1000])
ylabel('DF/F')
xlabel('time from stim on (ms)')
legend({'best','worst'},'location','best')



function eBaseHigh_Callback(hObject, eventdata, handles)
% hObject    handle to eBaseHigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eBaseHigh as text
%        str2double(get(hObject,'String')) returns contents of eBaseHigh as a double


% --- Executes during object creation, after setting all properties.
function eBaseHigh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eBaseHigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eStimHigh_Callback(hObject, eventdata, handles)
% hObject    handle to eStimHigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eStimHigh as text
%        str2double(get(hObject,'String')) returns contents of eStimHigh as a double


% --- Executes during object creation, after setting all properties.
function eStimHigh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eStimHigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
