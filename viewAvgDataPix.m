function varargout = viewAvgDataPix(varargin)
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

% Last Modified by GUIDE v2.5 27-Jan-2015 15:53:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @viewAvgDataPix_OpeningFcn, ...
                   'gui_OutputFcn',  @viewAvgDataPix_OutputFcn, ...
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
function viewAvgDataPix_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
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
function varargout = viewAvgDataPix_OutputFcn(hObject, eventdata, handles) 
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
function b2Pfilep_Callback(hObject, eventdata, handles)
% hObject    handle to b2Pfilep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fn=uigetdir;
if ~isequal(fn,0)
    set(handles.e2Pfilep,'String',fn);
end


% --- Executes on button press in bAnafilep.
function bAnafilep_Callback(hObject, eventdata, handles)
% hObject    handle to bAnafilep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fn=uigetdir;
if ~isequal(fn,0)
    set(handles.eAnafilep,'String',fn);
end

% --- Executes on button press in bLoad.
function bLoad_Callback(hObject, eventdata, handles)
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
filename2P=fullfile(base2P,animalID,[animalID '_' unitID '_' expID '_trials.mat']);
filenameAna=fullfile(baseAna,animalID,[animalID '_u' unitID '_' expID '.analyzer']);

%disp(filename2P)

%check whether both files exist and throw error message if not
exist2P=exist(filename2P,'file');
existAna=exist(filenameAna,'file');

if exist2P~=2
    errordlg('2P file not found!');
end
if existAna~=2
    errordlg('Analyzer file not found!');
end

%proceed if both files exist
if exist2P==2 && existAna==2
    %load 2P file - this generates trial_acc and trial_n
    load(filename2P);
    
    %load analyzer file - this generates Analyzer
    load(filenameAna,'-mat');
    handles.Analyzer=Analyzer;
    
    %from Analyzer file, extract list of trials and conditions
    handles.triallist=getcondtrial(Analyzer);
    [handles.domains,handles.domval,handles.blankid]=getdomainvalue(Analyzer);
    
    %get nr of trials, and figure out whether baseline period was saved for
    %every trial
    Nt=getnotrials(Analyzer);
    handles.periodPerTrial=2;
    if length(trial_acc)~=2*Nt
        errordlg('Missing baseline period!')
    end
    %disp(handles.periodPerTrial)   
    
    %2P data needs to be inverted, and divided by number of frames
    for i=1:length(trial_acc)/handles.periodPerTrial
        td1=-double(trial_acc{(i-1)*2+1}/trial_n((i-1)*2+1));
        td2=-double(trial_acc{i*2}/trial_n(i*2));
        
        df=(td2-td1)./abs(td1);
        %additional hack (bug)
        df=fliplr(imrotate(df,-90));
        
        %compute dF/F for every trial
        handles.trialdF{i}=df;
    end
    
    %save one picture for anatomy (just average over a few images)
    handles.imAnat=zeros(size(handles.trialdF{1}));
    for i=1:6
        %hack
        tmp=-double(trial_acc{i}/trial_n(i));
        tmp=fliplr(imrotate(tmp,-90));
        handles.imAnat=handles.imAnat+tmp/6;
    end
    
    %populate selection for conditions
    set(handles.sCond1,'String',handles.domains);
    
    %enable plotting
    set(handles.bPlotAvg,'Enable','on');
    
    %notify that done with import
    helpdlg('Finished loading data!')
    
end %if
guidata(hObject, handles);


    


%%%%% PLOT PART OF THE GUI 

% --- Executes during object creation, after setting all properties.
function sCond1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sCond1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function eMod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eMod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





%%% PLOT AVERAGE


% --- Executes during object creation, after setting all properties.
function eLowpass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eLowpass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%plot average

% --- Executes on button press in bPlotAvg.
function bPlotAvg_Callback(hObject, eventdata, handles)
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
    modVal=str2num(get(handles.eMod,'String'));
    handles.primCondVal=unique(mod(handles.domval(:,primCondID),modVal));
end
handles.nCond=length(handles.primCondVal);


%filter data if selected
handles.FtrialdF=handles.trialdF; %need to do this so that filtering doesn't mess with the original data
if get(handles.sLowpass,'Value')==1
    Lwidth=str2double(get(handles.eLowpass,'String'));
    L = fspecial('gaussian',size(handles.FtrialdF{1}),Lwidth);
    hh = ifft2(abs(fft2(L)));   %Eliminate phase information

    for i=1:length(handles.FtrialdF)
        valmap=handles.trialdF{i};
        id = find(isnan(valmap));
        valmap(id) = 0;
        handles.FtrialdF{i} = ifft2(abs(fft2(hh)).*fft2(valmap));
    end
end

%now compute averages across conditions

handles.dFData=zeros([size(handles.FtrialdF{1}) handles.nCond]);
 

for c=1:handles.nCond
    %get condition indices for each condition
    if get(handles.sMod,'Value')==0
        condID=find(handles.domval(:,primCondID)==handles.primCondVal(c));
    else
        condID=find(mod(handles.domval(:,primCondID),modVal)==handles.primCondVal(c));
    end
    
    trialID=find(ismember(handles.triallist,condID)); %find the matching trials
    %necessary because if there are multiple parameters, the same
    %parameter value shows up multiple times in the condition list
        
    %compute dF/F average
    stemp=zeros(size(handles.FtrialdF{1}));
    for i=1:length(trialID)
        handles.dFData(:,:,c)=handles.dFData(:,:,c)+handles.FtrialdF{trialID(i)}/length(trialID);
        stemp=stemp+handles.FtrialdF{trialID(i)}.^2/length(trialID); %for standard deviation
    end
    handles.sdevdFData(:,:,c)=sqrt(stemp-squeeze(handles.dFData(:,:,c)).^2);
        
end %for cond
    
    
%determine value in blank condition
handles.dFBlank=zeros(size(handles.FtrialdF{1}));
trialID=find(ismember(handles.triallist,handles.blankid)); %find blank trials
for i=1:length(trialID)
    handles.dFBlank=handles.dFBlank+handles.FtrialdF{trialID(i)}/length(trialID);
end


%now find the maximum condition at every pixel 
[maxResp,handles.maxCond]=max(handles.dFData,[],3);
[minResp,minCond]=min(handles.dFData,[],3);


%also compute response strength if - difference between best and worst
%condition
mag=maxResp-minResp;
mag = mag-min(mag(:));
mag = mag/max(mag(:));

%convert into color
if get(handles.sCirc,'Value') | get(handles.sMod,'Value')
    nColors=handles.nCond+1;
else
    nColors=handles.nCond;
end
maxCondImg=(handles.maxCond-1)/(nColors-1); %scale between 0 and 1, making sure that the first condition repeats
maxCondImg=round(maxCondImg*63+1);

%plot
if ~isfield(handles,'plotAHandle') || ~ishandle(handles.plotAHandle)  %make sure figure isn't already open
    handles.plotAHandle=figure;
end

%determine which plot options were selected
scaleMag=get(handles.sScale,'Value');
anatFlag=get(handles.sAnatomy,'Value');



if anatFlag==1 %compute anatomy picture 
    imAnat=handles.imAnat;
    imAnat=imAnat-min(imAnat(:));
    imAnat=imAnat/max(imAnat(:));
    
    imAnat(:,:,2) = imAnat;
    imAnat(:,:,3) = imAnat(:,:,1);
end


figure(handles.plotAHandle);

if get(handles.sCirc,'Value') & get(handles.sMod,'Value')==0
    ytlabel=[handles.primCondVal;handles.primCondVal(1)];
elseif get(handles.sMod,'Value')
    ytlabel=[handles.primCondVal;modVal];
else
    ytlabel=handles.primCondVal;
end

if anatFlag==0 & scaleMag==0
    image(1:length(maxCondImg(1,:)),1:length(maxCondImg(:,1)),maxCondImg,'CDataMapping','direct');
    if ~get(handles.sCirc,'Value')
        colormap('jet');
    else
        colormap('hsv');
    end
    cbar=colorbar;
    yt=linspace(1,64,nColors);
    set(cbar,'YTick',yt,'YTicklabel',ytlabel)
elseif anatFlag==0 & scaleMag==1
    image(1:length(maxCondImg(1,:)),1:length(maxCondImg(:,1)),maxCondImg,'CDataMapping','direct','AlphaData',mag,'AlphaDataMapping','none');
    axis image;
    if ~get(handles.sCirc,'Value')
        colormap('jet');
    else
        colormap('hsv');
    end
    cbar=colorbar;
    yt=linspace(1,64,nColors);
    set(cbar,'YTick',yt,'YTicklabel',ytlabel)
elseif anatFlag==1 
    if ~get(handles.sCirc,'Value')
        cid = jet;
    else
        cid = hsv;
    end
    
    imout = zeros([size(maxCondImg) 3]);
    for i = 1:size(maxCondImg,1)
        for j = 1:size(maxCondImg,2)   
            if scaleMag==1
                imout(i,j,:) = mag(i,j)*cid(maxCondImg(i,j),:);
            else
                imout(i,j,:) = cid(maxCondImg(i,j),:);
            end
        end
    end
    
    imout = imout + imAnat;
    imout = imout/max(imout(:));
    image(imout,'CDataMapping','direct','AlphaDataMapping','none');
    cbar=colorbar;
    yt=linspace(1,64,nColors);
    set(cbar,'YTick',yt,'YTicklabel',ytlabel)
end
    
%add user data to the figure so that we can retrieve it using the data
%cursor
set(handles.plotAHandle,'UserData',handles);

% enable data cursor
datacursormode on
dcm_obj=datacursormode(handles.plotAHandle);
set(dcm_obj,'DisplayStyle','window','SnapToDataVertex','on','UpdateFcn',@myupdatefcn);

guidata(hObject, handles);



%%%%% callback function for mouse clicks on average data
function txt = myupdatefcn(~,event_obj)
%get pointer position
pos = round(get(event_obj,'Position'));


%get data
h=get(gcf,'UserData');


%generate figure
fh=figure(99);



%plot responses for the conditions
if get(h.sMod,'Value')
    modVal=str2num(get(h.eMod,'String'));
    xval=[h.primCondVal;modVal];
    yval=squeeze(h.dFData(pos(2),pos(1),:));
    yval=[yval;yval(1)];
    erval=squeeze(h.sdevdFData(pos(2),pos(1),:));
    erval=[erval;erval(1)];
else
    xval=h.primCondVal;
    yval=squeeze(h.dFData(pos(2),pos(1),:));
    erval=squeeze(h.sdevdFData(pos(2),pos(1),:));
end

primCondID=get(h.sCond1,'Value');


errorbar(xval,yval,erval,'o-r')
hold on
plot([xval(1) xval(end)],ones(1,2).*squeeze(h.dFBlank(pos(2),pos(1))),'k--')
hold off
legend({'stim';'blank'});
title(['X: ' num2str(pos(1)) ', Y: ' num2str(pos(2))])
set(gca,'XLim',[xval(1) xval(end)])
ylabel('DF/F')
xlabel(h.domains{primCondID})


%update text window
txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))],...
       ['Pref. Cond: ',num2str(mc)]};
