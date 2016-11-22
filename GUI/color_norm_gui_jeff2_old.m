function varargout = color_norm_gui_jeff2(varargin)
% COLOR_NORM_GUI MATLAB code for color_norm_gui.fig
%      COLOR_NORM_GUI, by itself, creates a new COLOR_NORM_GUI or raises the existing
%      singleton*.
%
%      H = COLOR_NORM_GUI returns the handle to a new COLOR_NORM_GUI or the handle to
%      the existing singleton*.
%
%      COLOR_NORM_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COLOR_NORM_GUI.M with the given input arguments.
%
%      COLOR_NORM_GUI('Property','Value',...) creates a new COLOR_NORM_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before color_norm_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to color_norm_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help color_norm_gui

% Last Modified by GUIDE v2.5 22-Jul-2016 09:53:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @color_norm_gui_jeff2_OpeningFcn, ...
                   'gui_OutputFcn',  @color_norm_gui_jeff2_OutputFcn, ...
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


% --- Executes just before color_norm_gui is made visible.
function color_norm_gui_jeff2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to color_norm_gui (see VARARGIN)

% Choose default command line output for color_norm_gui
handles.output = hObject;
handles.count = 1;
norm_pairs = readtable('GUI color norm pair.txt', 'Delimiter', ',', 'ReadVariableNames', true);
rand_index = randperm(height(norm_pairs));
handles.norm_pairs = norm_pairs(rand_index, :);
handles.count_colors = {'g','r','y','m','c','b','w'};
handles.color_index = repmat([1:7],1,5);

im2 = imread(char(norm_pairs{handles.count, 'Target'})); 
im3 = imread(char(norm_pairs{handles.count, 'm1'})); 
im4 = imread(char(norm_pairs{handles.count, 'm2'}));
axes(handles.axes3);imshow(im3);
axes(handles.axes4);imshow(im4);
handles.axes14.Color = handles.count_colors{handles.color_index(handles.count)};

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes color_norm_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = color_norm_gui_jeff2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if strcmp(eventdata.Key,'numpad1')
    handles.axes5.Color=[0.69,0.69,0.973];
    handles.axes6.Color=[0.941,0.941,0.941];
    handles.axes15.Color=[0.941,0.941,0.941];
    handles.user_select = 1;
    handles.user_choices{handles.count} = handles.user_select;
    handles.count = handles.count + 1;
    handles.axes14.Color = handles.count_colors{handles.color_index(handles.count)};
    set(handles.countlabel, 'String', strcat('#',num2str(handles.count),' out of 300'));
    if handles.count > height(handles.norm_pairs)
        % stop the GUI
        handles.norm_pairs.Results = handles.user_choices';
        writetable(handles.norm_pairs,'user_results_jeff.txt');
        [msg_text, table_rank] = rank_methods();
        avg_inconst = prcm_incons()
        h = msgbox(sprintf('You have finished the test with inconsistency score of %.2f\n %s',...
            avg_inconst,msg_text),'Done');        
        save('table_rank.mat', 'table_rank');
        uiwait(h);
        close all;
        return;
    end
    %im1 = imread(char(handles.norm_pairs{handles.count, 'Source'}));
    im2 = imread(char(handles.norm_pairs{handles.count, 'Target'}));
    im3 = imread(char(handles.norm_pairs{handles.count, 'm1'}));
    im4 = imread(char(handles.norm_pairs{handles.count, 'm2'}));
    %axes(handles.axes1);imshow(im1);
    %axes(handles.axes2);imshow(im2);
    axes(handles.axes3);imshow(im3);
    axes(handles.axes4);imshow(im4);
elseif strcmp(eventdata.Key,'numpad2')
    handles.axes6.Color=[0.69,0.69,0.973];
    handles.axes5.Color=[0.941,0.941,0.941];
    handles.axes15.Color=[0.941,0.941,0.941];
    handles.user_select = 0;
    handles.user_choices{handles.count} = handles.user_select;
    handles.count = handles.count + 1;
    handles.axes14.Color = handles.count_colors{handles.color_index(handles.count)};
    set(handles.countlabel, 'String', strcat('#',num2str(handles.count),' out of 300'));
    if handles.count > height(handles.norm_pairs)
        % stop the GUI
        handles.norm_pairs.Results = handles.user_choices';
        writetable(handles.norm_pairs,'user_results_jeff.txt');
        [msg_text, table_rank] = rank_methods();
        avg_inconst = prcm_incons()
        h = msgbox(sprintf('You have finished the test with inconsistency score of %.2f\n %s',...
            avg_inconst,msg_text),'Done');        
        save('table_rank.mat', 'table_rank');
        uiwait(h);
        close all;
        return;
    end
    %im1 = imread(char(handles.norm_pairs{handles.count, 'Source'}));
    im2 = imread(char(handles.norm_pairs{handles.count, 'Target'}));
    im3 = imread(char(handles.norm_pairs{handles.count, 'm1'}));
    im4 = imread(char(handles.norm_pairs{handles.count, 'm2'}));
    %axes(handles.axes1);imshow(im1);
    %axes(handles.axes2);imshow(im2);
    axes(handles.axes3);imshow(im3);
    axes(handles.axes4);imshow(im4);
elseif strcmp(eventdata.Key,'numpad3')
    handles.axes15.Color=[0.69,0.69,0.973];
    handles.axes6.Color=[0.941,0.941,0.941];
    handles.axes5.Color=[0.941,0.941,0.941];
    handles.user_select = -1;
    handles.user_choices{handles.count} = handles.user_select;
    handles.count = handles.count + 1;
    handles.axes14.Color = handles.count_colors{handles.color_index(handles.count)};
    set(handles.countlabel, 'String', strcat('#',num2str(handles.count),' out of 300'));
    if handles.count > height(handles.norm_pairs)
        % stop the GUI
        handles.norm_pairs.Results = handles.user_choices';
        writetable(handles.norm_pairs,'user_results_jeff.txt');
        [msg_text, table_rank] = rank_methods();
        avg_inconst = prcm_incons()
        h = msgbox(sprintf('You have finished the test with inconsistency score of %.2f\n %s',...
            avg_inconst,msg_text),'Done');
        save('table_rank.mat', 'table_rank');
        uiwait(h);
        close all;
        return;
    end
    %im1 = imread(char(handles.norm_pairs{handles.count, 'Source'}));
    im2 = imread(char(handles.norm_pairs{handles.count, 'Target'}));
    im3 = imread(char(handles.norm_pairs{handles.count, 'm1'}));
    im4 = imread(char(handles.norm_pairs{handles.count, 'm2'}));
    %axes(handles.axes1);imshow(im1);
    %axes(handles.axes2);imshow(im2);
    axes(handles.axes3);imshow(im3);
    axes(handles.axes4);imshow(im4);
elseif strcmp(eventdata.Key, 'backspace')
    if handles.count>1
        handles.count = handles.count - 1;
        %im1 = imread(char(handles.norm_pairs{handles.count, 'Source'}));
        %im2 = imread(char(handles.norm_pairs{handles.count, 'Target'}));
        im3 = imread(char(handles.norm_pairs{handles.count, 'm1'}));
        im4 = imread(char(handles.norm_pairs{handles.count, 'm2'}));
        %axes(handles.axes1);imshow(im1);
        %axes(handles.axes2);imshow(im2);
        axes(handles.axes3);imshow(im3);
        axes(handles.axes4);imshow(im4);
    end

else
    return;
end
guidata(hObject, handles);


function [msg_txt, table_rank] = rank_methods()

norm_pairs = readtable('user_results_jeff.txt'); 
% user_choices = randi(3,size(norm_pairs,1),1) - 2; 
% norm_pairs.user_choices = user_choices; 

count_win = zeros(6,1);

for mm = 1:6
    count_win(mm) = sum(norm_pairs.m1_num == mm & norm_pairs.Results == 1) + ...
        sum(norm_pairs.m2_num == mm & norm_pairs.Results == 0) + ...
        0.5*sum((norm_pairs.m1_num == mm | norm_pairs.m2_num == mm) ...
        & norm_pairs.Results == 0);
end

avg_count_win = count_win'/20; 
method_names = {'Luong', 'Macenko', 'Reinhard','Khan', 'Vahadane', 'Vahadane_Fast'};
rank_method = [(1:6); avg_count_win];
[~, indx] = sort(rank_method(2,:), 'descend');
table_rank = rank_method(:,indx);
% table_rank = table(avg_count_win(1), avg_count_win(2),avg_count_win(3), avg_count_win(4), avg_count_win(5),avg_count_win(1,6)); 
% table_rank.Properties.VariableNames = {method_names{rank_method(1,1)}, ...
%     method_names{rank_method(1,2)}, method_names{rank_method(1,3)},...
%     method_names{rank_method(1,4)}, method_names{rank_method(1,5)},...
%     method_names{rank_method(1,6)}};
%message_box = strcat('Here is the ranking:\n', table_rank();
msg_txt= sprintf('Here is the ranking\n Method %s with score %.2f\n Method %s with score %.2f\n Method %s with score %.2f\n Method %s with score %.2f\n Method %s with score %.2f\n Method %s with score %.2f\n', ...
    method_names{table_rank(1,1)},table_rank(2,1),...
    method_names{table_rank(1,2)},table_rank(2,2),...
    method_names{table_rank(1,3)},table_rank(2,3),...
    method_names{table_rank(1,4)},table_rank(2,4),...
    method_names{table_rank(1,5)},table_rank(2,5),...
    method_names{table_rank(1,6)},table_rank(2,6));


function avg_inconst = prcm_incons()
T = readtable('user_results_jeff.txt'); 
% test the logic
count_inconst = 0;
for tt = 1:2
    for ss = 1:10
       indx_pair = (T.source_num == ss) & (T.target_num == tt);
       pairwise_comp_mat = sparse([T.m1_num(indx_pair) T.m2_num(indx_pair)],...
          [T.m2_num(indx_pair) T.m1_num(indx_pair)], [T.Results(indx_pair) T.Results(indx_pair)]);
       pairwise_comp_mat = full(pairwise_comp_mat);
       scores = zeros(6,1);
       for mm = 1:6
           %scores(mm) = sum(pairwise_comp_mat(:,mm) == 1) + 0.5*sum(pairwise_comp_mat(:,mm) == 0);
           scores(mm) = sum(T.m1_num(indx_pair) == mm & T.Results(indx_pair) == 1) + ...
               sum(T.m2_num(indx_pair) == 1 & T.Results(indx_pair) == -1) + ...
               0.5* sum(T.m1_num(indx_pair) == mm & T.Results(indx_pair) == 0); 
       end
       [sorted_scores, sorted_indx] = sort(scores,'descend');
       for m1 = 1:5
           for m2 = (m1+1):6
               m1_num = sorted_indx(m1);
               m2_num = sorted_indx(m2);
               if m1_num > max(size(pairwise_comp_mat)) || m2_num > max(size(pairwise_comp_mat))
                  disp('something wrong'); 
                  %uiwait;
                  break;
               end
               comp_result = pairwise_comp_mat(m1_num,m2_num);
               if (comp_result == -1) || (comp_result == 0  && ...
                       sorted_scores(m1_num) ~= sorted_scores(m2_num))
                  %fprintf('method %d score is %.2f, method %d score is %.2f\n', ...
                  %    m1_num, sorted_scores(m1_num),m2_num,  sorted_scores(m2_num));
                  count_inconst = count_inconst + 1; 
               end
           end
       end
    end 
end
avg_inconst = count_inconst/300;