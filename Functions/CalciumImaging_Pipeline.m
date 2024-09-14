classdef CalciumImaging_Pipeline < handle
    %% CalciumImaging_Pipeline -
    %
    % TODO LIST
    %    - Perform some safety checks for potential frames dropped just in case
    %    (we are using the timing infos rather than the pure frames number
    %    and period to deduce the stamps, but... just in case)
    %    - Create a static function to rederive the whole processing tree
    %    onto the F: drive from the archived raw data + final metadata file
    %
    %
    %     Copyright (C) 2024 Jérémy Signoret-Genest, DefenseCircuitsLab
    %     Original version: -/-/2020
    %     First release:    07/01/2021
    %     Current version:  21/02/2024
    %
    %     Changelog:
    %       - 22/03/2020:
    %           . Added options for archiving
    %       - 21/03/2021:
    %           . Added a GUI element to validate one extraction as the one
    %           to be used; in pratice it exports the calcium data to a
    %           _CalciumData.mat file in the session folder, adds set a
    %           tag to true in the corresponding CNMFE file, and delete the
    %           other "processing branches"
    %           . Added a function to be called externally, to check all
    %           the sessions present down a selected folder: if the data
    %           was exported/validated, everything but the motion corrected
    %           movie will be archived
    %       - 10/03/2021:
    %           . added a GUI element to change the ring_radius value that
    %           is used for background estimation (normally derived from
    %           gSiz, but in some circumstances, it works better with
    %           different values: now easier to try out)
    %       - 21/01/2021:
    %           . completed missing operations in callbacks for the
    %           post-processing (traces were not updated properly)
    %           . added a "back-up" for the initial CNMFE output that is
    %           preserved despite post-processing and can be restored
    %           (added the UI elements to restore it)
    %       - 25/10/2023:
    %           . added support for open miniscope data
    %
    %     This program is free software: you can redistribute it and/or modify
    %     it under the terms of the GNU General Public License as published by
    %     the Free Software Foundation, either version 3 of the License, or
    %     (at your option) any later version.
    %
    %     This program is distributed in the hope that it will be useful,
    %     but WITHOUT ANY WARRANTY; without even the implied warranty of
    %     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %     GNU General Public License for more details.
    %
    %     You should have received a copy of the GNU General Public License
    %     along with this program.  If not, see <https://www.gnu.org/licenses/>.

    properties(SetAccess = private, GetAccess = public, Hidden = true) % We might want to access the object(s) manually... But only if we know what we're doing!

    end


    properties(SetAccess = private, GetAccess = public, Hidden = false)
        DefaultParameters
        Handles
        Parameters
        PlayRate = 1;
        StartPath
    end

    properties(SetAccess = private, GetAccess = private, Hidden = true)
        CurrentPlayers
        CurrentPlayingStatus
        CurrentTab
        CurrentTime = 0;
        Dragging = false;
        JustLoaded = false;
        LastOperation      % To hold the last operation done in post-processing step (to allow for cancellation)
        MergeIndex         % Matrix to hold the indexes of merge neurons, to delete them from the actual data at the very end (until then they are just disabled and the merged neurons are added at the end)
        Playing = false;
        PlayingSingle = false;
        UIResponsiveness
    end

    methods
        % Constructor
        function obj = CalciumImaging_Pipeline
            %CalciumImaging_Pipeline - Construct an instance of this class
            addpath(genpath('F:\MATLAB\Common\CalciumImaging'))
            addpath(genpath('C:\Program Files\Inscopix'))
            %% Parameters
            obj.DefaultParameters.PreProcessing = struct(...
                'Basename',                 [],...
                'Cropping',                 [],...      % Rectangle coordinates used to crop during preprocessing (handled by API)
                'd1',                       [],...      % d1 for easy access
                'd2',                       [],...      % d2 for easy access
                'Date',                     [],...      % Date of the recording
                'Duration',                 [],...      % Duration (deduced from frames number and period)
                'Exposure',                 [],...      % Exposure time
                'Fix_Pixels',               true,...    % To fix defective pixels (careful, median filtering, so not several times)
                'Focus',                    [],...
                'ForkNumber',               0,...       % In case of several preprocessings
                'FramesNum',                [],...      % Number of frames for easy access
                'Gain',                     [],...      % Gain
                'MainFolder',               [],...
                'Miniscope',                [],...
                'Period',                   [],...      % Period for acquisition
                'Power',                    [],...      % Led power (mW)
                'PreProcessingFolder',      [],...      % Folder containg the preprocessed files (.isxd and .h5)
                'PreProcessingFile_isxd',   [],...      % .isxd output file from the inscopix API (inscopix recordings)
                'PreProcessingFile',        [],...      % .h5 file, the .isxd won't be used afterwards
                'RawFile',                  [],...      % Complete filepath to the raw file
                'Spatial_Downsampling',     2,...
                'Temporal_Downsampling',    1,...
                'TimeStamps',               []);        % Exact timestamps for each frame, from the API

            obj.DefaultParameters.PreFiltering = struct(...
                'Cropping',                 [],...      % Coordinates that can be used to derive the cropped mask
                'd1',                       [],...      % d1 for easy access
                'd2',                       [],...      % d2 for easy access
                'Duration',                 [],...      % Duration (deduced from frames number and period)
                'ForkNumber',               0,...       % In case of several prefilterings
                'FramesNum',                [],...      % Number of frames for easy access
                'Kernel',                   [],...      % Filter kernel derived from Sig1/2 and WindowSize
                'PreFilteringFolder',       [],...      % Folder containg the prefiltered file
                'PreFilteringFile',         [],...
                'SamplingRate',             [],...
                'Sig1',                     0.5,...     % Sigma for the first gaussian*
                'Sig2',                     2,...       % Sigma for the second gaussian*
                'Threshold',                [0],...      % Anything below Threshold low is set to NaN, and above to threshold high
                'TimeStamps',               [],...      % Exact timestamps for each frame, from the API
                'WindowSize',               30);        % Window size for the filter kernel*
            % * depends of course on the subsampling, and different FOV/populations can
            % have very different neurons' sizes

            obj.DefaultParameters.MotionCorrection = struct(...
                'bin_width',                250,...         % width of each bin
                'ForkNumber',               0,...           % In case of several motion corrections
                'grid_size',                [192,192,1],... % size of non-overlapping regions
                'init_batch',               250,...         % length of initial batch
                'iter',                     1,...           % number of data passes
                'max_dev',                  [32,32,1],...   % maximum deviation of patch shift from rigid shift
                'max_shift',                [80,80,1],...   % maximum rigid shift in each direction
                'MotionCorrectionFolder',   [],...          % Folder containg the motion corrected (raw) file
                'MotionCorrectionFile',     [],...          % .h5 file
                'overlap_pre',              [96,96,1],...   % size of overlapping region
                'overlap_post',             [96,96,1]);     % size of overlapping region after upsampling


            obj.DefaultParameters.CNMFE = struct(...
                'maxRAM',                   64,...              % maximum RAM to use during CNMFE
                'max_tau',                  10,...              % maximum decay time (seconds)
                'detrend_nk',               1,...               % detrending the slow fluctuation. usually 1 is fine (no detrending) / Integers <= TotalDuration/30
                'merge_thr_spatial',        [0.8, 0.1, -inf],...% merge components with highly correlated spatial shapes (corr=0.8) and small temporal correlations (corr=0.1)
                'merge_thr',                0.4,...             % thresholds for merging neurons; [spatial overlap ratio, temporal correlation of calcium traces, spike correlation]
                'CorrSubsampling',          1,...               % Subsampling to process the correlation/PNR pictures (every nth frame is taken to speed up)
                'min_corr',                 0.9,...             % minimum local correlation for a seeding pixel
                'min_pnr',                  14,...              % minimum peak-to-noise ratio for a seeding pixel
                'Crop',                     [],...              % rectangle coordinates to crop the movie before CNMFE (to keep only what's relevant and get rid of motion correction border artifacts)
                'Cell',                     [],...              % coordinates of the ROI used to get the average cell size - more for convenience and replotting than real analytical use
                'MP_Sig',                   [],...              % sigma value used to filter the data before processing the max projection - again, more for convenience and replotting
                'MP_Sub',                   [],...              % factor for frame skipping that was used to process maximum projection - again for convenience and reloading
                'DataFile',                 [],...              % file to use CNMFE on; basically if we crop, it's the cropped motion corrected movie, otherwise, just the motion corrected movie
                'CNMFEFolder',              [],...              % Folder containing CNMFE results; everything is generated by CNMFE itself
                'Constrained',              true,...            % If set to false, no negativity constraint (deconvolution is disabled)
                'Processed',                false);             % Flag to know when CNMFE was already run
            % Most variables can be set to the CNMFE parameters structure directly, but keeping those for which we want to change the default values here might be easier than going to CNMFE

            %% Main GUI
            set(0,'units','pixels')
            Scrsz = get(0,'ScreenSize');
            obj.Handles.MainFigure = figure('Position',[1 40 Scrsz(3) Scrsz(4)-70],...
                'MenuBar','none','ToolBar','none','Color',[0.25,0.25,0.25],'Renderer','painters');
            h = zoom(obj.Handles.MainFigure);
            h.ActionPostCallback = {@(f,src)obj.ZoomPostCallback(f,src)};

            obj.Handles.TreePanel.Axes = axes('Position',[0.01 0.025 0.145 0.95],...
                'Color',[0.45,0.45,0.45],'Box','on','BoxStyle','full',...
                'LineWidth',3);
            obj.Handles.TreePanel.Axes.XColor = [0.25 0.25 0.25];
            obj.Handles.TreePanel.Axes.YColor = [0.25 0.25 0.25];
            obj.Handles.TreePanel.Axes.XTick = [];
            obj.Handles.TreePanel.Axes.YTick = [];

            obj.Handles.MainDisplay = axes('Position',[0.17 0.025 0.82 0.918],...
                'Color',[0.2,0.2,0.2],'Box','on','BoxStyle','full',...
                'LineWidth',3);
            obj.Handles.MainDisplay.XColor = 'k';
            obj.Handles.MainDisplay.YColor = 'k';
            obj.Handles.MainDisplay.XTick = '';
            obj.Handles.MainDisplay.YTick = '';

            % Tree panel
            obj.Handles.TreePanel.LoadSession_Button = uicontrol('Style','pushbutton',...
                'Units','Normalized','Position',[0.035 0.035 0.095 0.04],...
                'String','Load session',...
                'FontSize',14,'FontName','Arial','FontWeight','b',...
                'BackgroundColor',[0.15 0.15 0.15],'ForegroundColor',[0.6 0.6 0.6],...
                'HorizontalAlignment','center',...
                'Callback',{@(~,~)obj.LoadSession},...
                'Enable','on');
            obj.Handles.TreePanel.SessionName_Legend = uicontrol('Style','text',...
                'Units','Normalized','Position',[0.015 0.925 0.135 0.04],...
                'String','',...
                'FontSize',13,'FontName','Arial','FontWeight','b',...
                'BackgroundColor',[0.45,0.45,0.45],'ForegroundColor','k',...
                'HorizontalAlignment','center');


            obj.Handles.Tabs.PreProcessing = uicontrol('Style','text','Units','normalized',...
                'Position',[0.17 0.94 0.12 0.035],'BackgroundColor',[0.25,0.25,0.25],...
                'String','PreProcessing','FontSize',16,'FontWeight','b',...
                'UserData','1','Enable','Inactive',...
                'ButtonDownFcn',@(src,evt)obj.TabCallback(src,evt));
            obj.Handles.Tabs.PreFiltering = uicontrol('Style','text','Units','normalized',...
                'Position',[0.17+0.12 0.94 0.12 0.035],'BackgroundColor',[0.25,0.25,0.25],...
                'String','PreFiltering','FontSize',16,'FontWeight','b',...
                'UserData','0','Enable','Inactive',...
                'ButtonDownFcn',@(src,evt)obj.TabCallback(src,evt));
            obj.Handles.Tabs.MotionCorrection = uicontrol('Style','text','Units','normalized',...
                'Position',[0.17+0.12*2 0.94 0.12 0.035],'BackgroundColor',[0.25,0.25,0.25],...
                'String','MotionCorrection','FontSize',16,'FontWeight','b',...
                'UserData','0','Enable','Inactive',...
                'ButtonDownFcn',@(src,evt)obj.TabCallback(src,evt));
            obj.Handles.Tabs.CNMFE = uicontrol('Style','text','Units','normalized',...
                'Position',[0.17+0.12*3 0.94 0.12 0.035],'BackgroundColor',[0.25,0.25,0.25],...
                'String','CNMF-E','FontSize',16,'FontWeight','b',...
                'UserData','0','Enable','Inactive',...
                'ButtonDownFcn',@(src,evt)obj.TabCallback(src,evt));
            obj.Handles.Tabs.PostProcessing = uicontrol('Style','text','Units','normalized',...
                'Position',[0.17+0.12*4 0.94 0.12 0.035],'BackgroundColor',[0.25,0.25,0.25],...
                'String','PostProcessing','FontSize',16,'FontWeight','b',...
                'UserData','0','Enable','Inactive',...
                'ButtonDownFcn',@(src,evt)obj.TabCallback(src,evt));
            obj.Handles.Tabs.Visualisation = uicontrol('Style','text','Units','normalized',...
                'Position',[0.872 0.94 0.12 0.035],'BackgroundColor',[0.25,0.25,0.25],...
                'String','Visualisation','FontSize',16,'FontWeight','b',...
                'UserData','0','Enable','Inactive',...
                'ButtonDownFcn',@(src,evt)obj.TabCallback(src,evt));
            obj.TabAppearence;

            % Load default preprocessing parameters - will be overwritten,
            % just for the GUI initialization
            obj.Parameters.PreProcessing = obj.DefaultParameters.PreProcessing;

            % Set default layout (Preprocessing)
            obj.PreProcessing_Layout;
        end

    end

    methods(Hidden)
        % Starting path
        function SetStartPath(obj)



        end

        % Disable all uicontrols
        function DisableAll(obj)
            % Enable
            AllEnable = findall(obj.Handles.MainFigure,'-property', 'Enable');
            while strcmpi(class(AllEnable),'matlab.ui.control.WebComponent')
                % There is a ROI menu last - that's a problem
                AllEnable = AllEnable(1:end-1);
            end
            obj.UIResponsiveness.Enable = {AllEnable.Enable};
            set(AllEnable,'Enable','inactive')

            % Taking care of the other callbacks
            AllEnable = findall(obj.Handles.MainFigure.Children(:),'-property', 'ButtonDownFcn');
            obj.UIResponsiveness.CB = cell(numel(AllEnable),1);
            for C = 1 : numel(AllEnable),
                obj.UIResponsiveness.CB{C} = AllEnable(C).ButtonDownFcn;
                set(AllEnable(C),'ButtonDownFcn',[])
            end
        end

        % Reenable all uicontrols (as they were before)
        function EnableAll(obj)
            % Enable
            AllEnable = findall(obj.Handles.MainFigure,'-property', 'Enable');
            arrayfun(@(x) set(AllEnable(x),'Enable',obj.UIResponsiveness.Enable{x}),1:numel(AllEnable))

            % Taking care of the other callbacks
            AllEnable = findall(obj.Handles.MainFigure.Children(:),'-property', 'ButtonDownFcn');
            for C = 1 : numel(AllEnable),
                set(AllEnable(C),'ButtonDownFcn',obj.UIResponsiveness.CB{C})
            end
        end

        % Session/Step choice (file choice)
        function LoadSession(obj)
            % Load file
            if isempty(obj.StartPath),
                Experimenters = DataBase.Lists.GetList('Experimenters');
                CurrentExperimenter = Experimenters(strcmpi(Experimenters(:,2),getenv('username')),1);
                [File, Path] = uigetfile({'*_StepMeta.mat','Meta file';
                    },'Please choose a file...',['T:\' CurrentExperimenter{1} '\Data\CalciumImaging\']);
            else
                [File, Path] = uigetfile({'*_StepMeta.mat','Meta file';
                    },'Please choose a file...','Please choose a file...',[obj.StartPath]);
            end
            if isempty(File) | File == 0,
                return
            end

            % Clean general properties (also stop movies)
            obj.PropertiesReset;

            % Read the metafile and set as parameters
            obj.Parameters = load(fullfile(Path,File));

            % Check whether files/folders are still in the same drive
            Fields = fieldnames(obj.Parameters);
            for Fi = 1  :numel(Fields)
                FieldsFi = fieldnames(obj.Parameters.(Fields{Fi}));
                FilesListFi = find(contains(FieldsFi,'File'));
                FoldersListFi = find(contains(FieldsFi,'Folder'));
                
                for FlFi = 1 : numel(FilesListFi)
                    FileFlFi = obj.Parameters.(Fields{Fi}).(FieldsFi{FilesListFi(FlFi)});
                    if exist(FileFlFi,'file')~=2
                        Drives = {'F','G','T','Z'};
                        for Dd = 1 : numel(Drives)
                            FileFlFi(1) = Drives{Dd};
                            if exist(FileFlFi,'file')==2
                                obj.Parameters.(Fields{Fi}).(FieldsFi{FilesListFi(FlFi)})(1) = Drives{Dd};
                                break
                            end
                        end
                    end
                end
                for FlFi = 1 : numel(FoldersListFi)
                    FolderFlFi = obj.Parameters.(Fields{Fi}).(FieldsFi{FoldersListFi(FlFi)});
                    if exist(FolderFlFi,'dir')~=7
                        Drives = {'F','G','T','Z'};
                        for Dd = 1 : numel(Drives)
                            FolderFlFi(1) = Drives{Dd};
                            if exist(FolderFlFi,'dir')==7
                                obj.Parameters.(Fields{Fi}).(FieldsFi{FoldersListFi(FlFi)})(1) = Drives{Dd};
                                break
                            end
                        end
                    end
                end
            end


            % Refresh tree

            % Reset the temporary plots for CNMFE
            if isfield(obj.Handles,'Values')
                obj.Handles = rmfield(obj.Handles,'Values');
            end

            % Simply reload the first tab: clicking different tabs will run
            % the normal checks
            obj.CurrentTab = 'ForceReload';
            obj.PreProcessing_Layout;
            obj.JustLoaded = true;

            % Update the tab display
            Tabs = fieldnames(obj.Handles.Tabs);
            for F = 1 : numel(Tabs),
                obj.Handles.Tabs.(Tabs{F}).UserData = '0';
            end
            obj.Handles.Tabs.PreProcessing.UserData = '1';
            obj.CurrentTab = 'PreProcessing';
            obj.TabAppearence;

            % Update the tree display
            obj.UpdateTree;
        end

        % Update tree
        function UpdateTree(obj)
            obj.Handles.TreePanel.SessionName_Legend.String = obj.Parameters.PreProcessing.Basename;
        end

        % Clean players (close readers properly if needed)
        function CleanMovies(obj)
            % As of now, everything can be deleted the same way, but this
            % might change in the future so let's keep the option from now
            for M = 1 : numel(obj.CurrentPlayers)
                switch class(obj.CurrentPlayers(M).Readers)
                    case 'isx.Movie'
                        % Apparently no special methods in the class
                        delete(obj.CurrentPlayers(M).Readers);
                    case 'VideoReader'
                        delete(obj.CurrentPlayers(M).Readers);
                    case 'h5Reader'
                        % Very simple class to read h5 files
                        % (Mainly to have the possibility to treat it the
                        % same way as avi or isxd)
                        delete(obj.CurrentPlayers(M).Readers);
                end
            end
        end

        % Load movie(s) to player(s)
        function LoadMovies(obj)
            if ~isempty(obj.CurrentPlayers)
                obj.CurrentTime = 0;
                MaxTimes = zeros(numel(obj.CurrentPlayers),1);
                for M = 1 : numel(obj.CurrentPlayers)
                    % Prepare reader and store in structure
                    switch obj.CurrentPlayers(M).Data
                        case 'Raw'
                            if contains(obj.Parameters.PreProcessing.RawFile,'.isxd')
                                obj.CurrentPlayers(M).Readers = isx.Movie.read(obj.Parameters.PreProcessing.RawFile);
                            elseif contains(obj.Parameters.PreProcessing.RawFile,'.h5')
                                obj.CurrentPlayers(M).Readers = h5Reader(obj.Parameters.PreProcessing.RawFile);
                            end

                            % Plot first frame
                            delete(obj.CurrentPlayers(M).Axes.Children(:))
                            if contains(obj.Parameters.PreProcessing.RawFile,'.isxd')
                                obj.CurrentPlayers(M).Plot = imagesc(obj.CurrentPlayers(M).Readers.get_frame_data(0),'Parent',obj.CurrentPlayers(M).Axes);
                            elseif contains(obj.Parameters.PreProcessing.RawFile,'.h5')
                                obj.CurrentPlayers(M).Plot = imagesc(obj.CurrentPlayers(M).Readers.GetFrame(1,[obj.Parameters.PreProcessing.d1 obj.Parameters.PreProcessing.d2]),'Parent',obj.CurrentPlayers(M).Axes);
                            end
                            MaxTimes = obj.Parameters.PreProcessing.Duration;


                            % Adjust limits to preserve ratio
                            if obj.Parameters.PreProcessing.d2/obj.Parameters.PreProcessing.d1 >= obj.CurrentPlayers(M).AbsolutePosition(3)/obj.CurrentPlayers(M).AbsolutePosition(4)
                                % Y needs to be adjusted
                                YDelta = obj.Parameters.PreProcessing.d2*obj.CurrentPlayers(M).AbsolutePosition(4)/obj.CurrentPlayers(M).AbsolutePosition(3) - obj.Parameters.PreProcessing.d1;
                                obj.CurrentPlayers(M).Axes.YLim = [1-0.5*YDelta obj.Parameters.PreProcessing.d1+0.5*YDelta];
                                obj.CurrentPlayers(M).Axes.XLim = [1 obj.Parameters.PreProcessing.d2];
                                obj.CurrentPlayers(M).Axes.YTick = obj.CurrentPlayers(M).Axes.YTick(obj.CurrentPlayers(M).Axes.YTick>=0 & obj.CurrentPlayers(M).Axes.YTick<=obj.Parameters.PreProcessing.d1);
                            else
                                % X needs to be adjusted
                                XDelta = obj.Parameters.PreProcessing.d1*obj.CurrentPlayers(M).AbsolutePosition(3)/obj.CurrentPlayers(M).AbsolutePosition(4) - obj.Parameters.PreProcessing.d2;
                                obj.CurrentPlayers(M).Axes.XLim = [1-0.5*XDelta obj.Parameters.PreProcessing.d2+0.5*XDelta];
                                obj.CurrentPlayers(M).Axes.YLim = [1 obj.Parameters.PreProcessing.d1];
                                obj.CurrentPlayers(M).Axes.XTick = obj.CurrentPlayers(M).Axes.XTick(obj.CurrentPlayers(M).Axes.XTick>=0 & obj.CurrentPlayers(M).Axes.XTick<=obj.Parameters.PreProcessing.d2);
                            end

                            %
                            %                             obj.CurrentPlayers(M).Axes.XLim = [1 obj.Parameters.PreProcessing.d2];
                            %                             obj.CurrentPlayers(M).Axes.YLim = [1 obj.Parameters.PreProcessing.d1];
                            obj.Handles.UIElements.MovieDisplay(1).TimeLine.XData = [0 0];
                            [~, RawFilename, ExtRaw] = fileparts(obj.Parameters.PreProcessing.RawFile);
                            obj.Handles.UIElements.BasenameLegend.String = [RawFilename,ExtRaw];
                        case 'PreProcessed'
                            if exist(obj.Parameters.PreProcessing.PreProcessingFile,'file')~=2
                                obj.Parameters.PreProcessing.PreProcessingFile = CalciumImaging_Pipeline.ToArchivePath(obj.Parameters.PreProcessing.PreProcessingFile_isxd);
                            end
                            %                                 obj.CurrentPlayers(M).Readers = isx.Movie.read(obj.Parameters.PreProcessing.PreProcessingFile_isxd);
                            obj.CurrentPlayers(M).Readers = h5Reader(obj.Parameters.PreProcessing.PreProcessingFile);
                            MaxTimes = obj.Parameters.PreFiltering.Duration;
                            % Plot first frame
                            delete(obj.CurrentPlayers(M).Axes.Children(:))
                            if obj.CurrentPlayers(M).Filter
                                if strcmpi(obj.CurrentPlayers(M).Kernel,'PreFiltering')
                                    obj.CurrentPlayers(M).Plot = imagesc(imfilter(obj.CurrentPlayers(M).Readers.GetFrame(1,[obj.Parameters.PreFiltering.d1 obj.Parameters.PreFiltering.d2]),obj.Parameters.PreFiltering.Kernel,'symmetric'),'Parent',obj.CurrentPlayers(M).Axes);
                                else
                                    obj.CurrentPlayers(M).Plot = imagesc(imfilter(obj.CurrentPlayers(M).Readers.GetFrame(1,[obj.Parameters.PreFiltering.d1 obj.Parameters.PreFiltering.d2]),obj.Parameters.Visualization.Kernel,'symmetric'),'Parent',obj.CurrentPlayers(M).Axes);
                                end
                            else
                                obj.CurrentPlayers(M).Plot = imagesc(obj.CurrentPlayers(M).Readers.GetFrame(1,[obj.Parameters.PreFiltering.d1 obj.Parameters.PreFiltering.d2]),'Parent',obj.CurrentPlayers(M).Axes);
                            end


                            % Adjust limits to preserve ratio
                            if obj.Parameters.PreFiltering.d2/obj.Parameters.PreFiltering.d1 >= obj.CurrentPlayers(M).AbsolutePosition(3)/obj.CurrentPlayers(M).AbsolutePosition(4),
                                % Y needs to be adjusted
                                YDelta = obj.Parameters.PreFiltering.d2*obj.CurrentPlayers(M).AbsolutePosition(4)/obj.CurrentPlayers(M).AbsolutePosition(3) - obj.Parameters.PreFiltering.d1;
                                obj.CurrentPlayers(M).Axes.YLim = [1-0.5*YDelta obj.Parameters.PreFiltering.d1+0.5*YDelta];
                                obj.CurrentPlayers(M).Axes.XLim = [1 obj.Parameters.PreFiltering.d2];
                                obj.CurrentPlayers(M).Axes.YTick = obj.CurrentPlayers(M).Axes.YTick(obj.CurrentPlayers(M).Axes.YTick>=0 & obj.CurrentPlayers(M).Axes.YTick<=obj.Parameters.PreFiltering.d1);
                            else
                                % X needs to be adjusted
                                XDelta = obj.Parameters.PreFiltering.d1*obj.CurrentPlayers(M).AbsolutePosition(3)/obj.CurrentPlayers(M).AbsolutePosition(4) - obj.Parameters.PreFiltering.d2;
                                obj.CurrentPlayers(M).Axes.XLim = [1-0.5*XDelta obj.Parameters.PreFiltering.d2+0.5*XDelta];
                                obj.CurrentPlayers(M).Axes.YLim = [1 obj.Parameters.PreFiltering.d1];
                                obj.CurrentPlayers(M).Axes.XTick = obj.CurrentPlayers(M).Axes.XTick(obj.CurrentPlayers(M).Axes.XTick>=0 & obj.CurrentPlayers(M).Axes.XTick<=obj.Parameters.PreFiltering.d2);
                            end

                            obj.Handles.UIElements.MovieDisplay(1).TimeLine.XData = [0 0];
                            [~, RawFilename, ExtRaw] = fileparts(obj.Parameters.PreProcessing.PreProcessingFile_isxd);
                            obj.Handles.UIElements.BasenameLegend.String = [RawFilename,ExtRaw];
                        case 'MotionCorrected'
                            %                             if exist(obj.Parameters.PreProcessing.PreProcessingFile_isxd,'file')~=2
                            %                                 obj.Parameters.MotionCorrection.MotionCorrectionFile = CalciumImaging_Pipeline.ToArchivePath(obj.Parameters.MotionCorrection.MotionCorrectionFile);
                            %                             end

                            obj.CurrentPlayers(M).Readers = h5Reader(obj.Parameters.MotionCorrection.MotionCorrectionFile);
                            MaxTimes = obj.Parameters.PreFiltering.Duration;
                            % Plot first frame
                            delete(obj.CurrentPlayers(M).Axes.Children(:))

                            if obj.CurrentPlayers(M).Filter
                                if strcmpi(obj.CurrentPlayers(M).Kernel,'PreFiltering')
                                    obj.CurrentPlayers(M).Plot = imagesc(imfilter(obj.CurrentPlayers(M).Readers.GetFrame(1,[obj.Parameters.PreFiltering.d1 obj.Parameters.PreFiltering.d2]),obj.Parameters.PreFiltering.Kernel,'symmetric'),'Parent',obj.CurrentPlayers(M).Axes);
                                elseif strcmpi(obj.CurrentPlayers(M).Kernel,'Visualization')
                                    obj.CurrentPlayers(M).Plot = imagesc(imfilter(obj.CurrentPlayers(M).Readers.GetFrame(1,[obj.Parameters.PreFiltering.d1 obj.Parameters.PreFiltering.d2]),obj.Parameters.Visualization.Kernel,'symmetric'),'Parent',obj.CurrentPlayers(M).Axes);
                                elseif strcmpi(obj.CurrentPlayers(M).Kernel,'PostProcessing')
                                    obj.CurrentPlayers(M).Plot = imagesc(imfilter(obj.CurrentPlayers(M).Readers.GetFrame(1,[obj.Parameters.PreFiltering.d1 obj.Parameters.PreFiltering.d2]),obj.Handles.Values.PostProcessing.Filtering.Kernel,'symmetric'),'Parent',obj.CurrentPlayers(M).Axes);
                                    if strcmpi(obj.Handles.UIElements.MovieDisplay(1).Plot.CLimMode,'auto')
                                        obj.Handles.UIElements.ThresholdLowPostProc_Edit.String = num2str(obj.CurrentPlayers(M).Axes.CLim(1));
                                        obj.Handles.UIElements.ThresholdHighPostProc_Edit.String = num2str(obj.CurrentPlayers(M).Axes.CLim(2));
                                        obj.Handles.Values.CaMovieFilteredCLim = obj.CurrentPlayers(M).Axes.CLim;
                                    end
                                end
                            else
                                obj.CurrentPlayers(M).Plot = imagesc(obj.CurrentPlayers(M).Readers.GetFrame(1, [obj.Parameters.PreFiltering.d1 obj.Parameters.PreFiltering.d2]),'Parent',obj.CurrentPlayers(M).Axes);
                                if strcmpi(obj.CurrentPlayers(M).Kernel,'PostProcessing')
                                    if strcmpi(obj.Handles.UIElements.MovieDisplay(1).Plot.CLimMode,'auto')
                                        obj.Handles.UIElements.ThresholdLowPostProc_Edit.String = num2str(obj.CurrentPlayers(M).Axes.CLim(1));
                                        obj.Handles.UIElements.ThresholdHighPostProc_Edit.String = num2str(obj.CurrentPlayers(M).Axes.CLim(2));
                                        obj.Handles.Values.CaMovieFilteredCLim = obj.CurrentPlayers(M).Axes.CLim;
                                    end
                                end
                            end


                            % Adjust limits to preserve ratio
                            if obj.Parameters.PreFiltering.d2/obj.Parameters.PreFiltering.d1 >= obj.CurrentPlayers(M).AbsolutePosition(3)/obj.CurrentPlayers(M).AbsolutePosition(4),
                                % Y needs to be adjusted
                                YDelta = obj.Parameters.PreFiltering.d2*obj.CurrentPlayers(M).AbsolutePosition(4)/obj.CurrentPlayers(M).AbsolutePosition(3) - obj.Parameters.PreFiltering.d1;
                                obj.CurrentPlayers(M).Axes.YLim = [1-0.5*YDelta obj.Parameters.PreFiltering.d1+0.5*YDelta];
                                obj.CurrentPlayers(M).Axes.XLim = [1 obj.Parameters.PreFiltering.d2];
                                obj.CurrentPlayers(M).Axes.YTick = obj.CurrentPlayers(M).Axes.YTick(obj.CurrentPlayers(M).Axes.YTick>=0 & obj.CurrentPlayers(M).Axes.YTick<=obj.Parameters.PreFiltering.d1);
                            else
                                % X needs to be adjusted
                                XDelta = obj.Parameters.PreFiltering.d1*obj.CurrentPlayers(M).AbsolutePosition(3)/obj.CurrentPlayers(M).AbsolutePosition(4) - obj.Parameters.PreFiltering.d2;
                                obj.CurrentPlayers(M).Axes.XLim = [1-0.5*XDelta obj.Parameters.PreFiltering.d2+0.5*XDelta];
                                obj.CurrentPlayers(M).Axes.YLim = [1 obj.Parameters.PreFiltering.d1];
                                obj.CurrentPlayers(M).Axes.XTick = obj.CurrentPlayers(M).Axes.XTick(obj.CurrentPlayers(M).Axes.XTick>=0 & obj.CurrentPlayers(M).Axes.XTick<=obj.Parameters.PreFiltering.d2);
                            end

                            drawnow
                            try
                            obj.Handles.UIElements.MovieDisplay(1).TimeLine.XData = [0 0];
                            end

                        otherwise
                            % Directly the file (e.g. behaviour)
                            % Should be .avi (previous UI), so we won't
                            % check here
                            obj.CurrentPlayers(M).Readers = VideoReader(obj.CurrentPlayers(M).Data);
                    end
                end
                drawnow
                for MN = 1 : numel(obj.Handles.UIElements.MovieDisplay)
                    obj.Handles.UIElements.MovieDisplay(MN).SourceXLim = obj.Handles.UIElements.MovieDisplay(MN).Plot.XLim;
                    obj.Handles.UIElements.MovieDisplay(MN).SourceYLim = obj.Handles.UIElements.MovieDisplay(MN).Plot.YLim;
                    if isfield(obj.Handles.UIElements.MovieDisplay(MN),'Box')
                        obj.Handles.UIElements.MovieDisplay(MN).Box.XLim = obj.Handles.UIElements.MovieDisplay(MN).Plot.XLim;
                        obj.Handles.UIElements.MovieDisplay(MN).Box.YLim = obj.Handles.UIElements.MovieDisplay(MN).Plot.YLim;
                    end
                end

                obj.Handles.UIElements.MovieDisplay(1).SliderTickAxes.XLim = [0 max(MaxTimes)];
                obj.Handles.UIElements.MovieDisplay(1).TimeLine.XData = [0 max(MaxTimes)];
                obj.Handles.UIElements.MovieDisplay(1).Slider.XData = obj.CurrentTime * [1 1];
                obj.Handles.UIElements.CurrentTimeEdit.String = num2str(1/100 * round(obj.CurrentTime*100),'%.2f');

            end
        end

        % Play movies
        function Play_CB(obj)
            if ~isempty(obj.CurrentPlayers)
                if obj.Playing,
                    obj.Playing = false;
                    obj.Handles.UIElements.PlayButton.String = 'Play';
                    drawnow
                else
                    obj.Playing = true;
                    obj.Handles.UIElements.PlayButton.String = 'Pause';
                    drawnow;
                    obj.PlayMovies;
                end
            end
        end

        function Slower_CB(obj)
            if (obj.PlayRate / 2) >= 0.125,
                obj.PlayRate = obj.PlayRate / 2;
                obj.Handles.UIElements.PlayRateLegend.String = ['x' num2str(1/1000 * round(obj.PlayRate*1000))];
            end
        end

        function Faster_CB(obj)
            if (obj.PlayRate / 2) <= 128,
                obj.PlayRate = obj.PlayRate * 2;
                obj.Handles.UIElements.PlayRateLegend.String  = ['x' num2str(1/100 * round(obj.PlayRate*100))];
            end
        end

        function PlayMovies(obj)
            if ~isempty(obj.CurrentPlayers)
                TicPlayer = tic;
                Ended = false(numel(obj.CurrentPlayers),1);
                while (obj.Playing || obj.PlayingSingle) && ~all(Ended)
                    TocPlayer = toc(TicPlayer);
                    TicPlayer = tic;
                    if ~obj.PlayingSingle
                        obj.CurrentTime = obj.CurrentTime + TocPlayer * obj.PlayRate;
                    end
                    for M = 1 : numel(obj.CurrentPlayers)
                        % Prepare reader and store in structure
                        switch obj.CurrentPlayers(M).Data
                            case 'Raw'
                                % Find index
                                IndxFrame = find(obj.Parameters.PreProcessing.TimeStamps >= obj.CurrentTime,1,'first');
                                if ~isempty(IndxFrame)
                                    if contains(obj.Parameters.PreProcessing.RawFile,'.isxd')
                                        obj.CurrentPlayers(M).Plot.CData = obj.CurrentPlayers(M).Readers.get_frame_data(IndxFrame-1);
                                    elseif contains(obj.Parameters.PreProcessing.RawFile,'.h5')
                                        obj.CurrentPlayers(M).Plot.CData = obj.CurrentPlayers(M).Readers.GetFrame(IndxFrame,[obj.Parameters.PreProcessing.d1 obj.Parameters.PreProcessing.d2]);
                                    end
                                else
                                    Ended(M) = true;
                                end
                            case 'PreProcessed'
                                % Find index
                                IndxFrame = find(obj.Parameters.PreFiltering.TimeStamps >= obj.CurrentTime,1,'first');
                                if ~isempty(IndxFrame)
                                    if obj.CurrentPlayers(M).Filter
                                        if strcmpi(obj.CurrentPlayers(M).Kernel,'PreFiltering'),
                                            obj.CurrentPlayers(M).Plot.CData = imfilter(obj.CurrentPlayers(M).Readers.GetFrame(IndxFrame,[obj.Parameters.PreFiltering.d1 obj.Parameters.PreFiltering.d2]),obj.Parameters.PreFiltering.Kernel,'symmetric');
                                        elseif strcmpi(obj.CurrentPlayers(M).Kernel,'Visualization')
                                            obj.CurrentPlayers(M).Plot.CData = imfilter(obj.CurrentPlayers(M).Readers.GetFrame(IndxFrame,[obj.Parameters.PreFiltering.d1 obj.Parameters.PreFiltering.d2]),obj.Parameters.Visualization.Kernel,'symmetric');
                                        elseif strcmpi(obj.CurrentPlayers(M).Kernel,'PostProcessing')
                                            obj.CurrentPlayers(M).Plot.CData = imfilter(obj.CurrentPlayers(M).Readers.GetFrame(IndxFrame,[obj.Parameters.PreFiltering.d1 obj.Parameters.PreFiltering.d2]),obj.Handles.Values.PostProcessing.Filtering.Kernel,'symmetric');
                                        end
                                    else
                                        obj.CurrentPlayers(M).Plot.CData = obj.CurrentPlayers(M).Readers.GetFrame(IndxFrame,[obj.Parameters.PreFiltering.d1 obj.Parameters.PreFiltering.d2]);
                                    end
                                else
                                    Ended(M) = true;
                                end
                            case 'MotionCorrected'
                                % Find index
                                IndxFrame = find(obj.Parameters.PreFiltering.TimeStamps >= obj.CurrentTime,1,'first');
                                if ~isempty(IndxFrame)
                                    if obj.CurrentPlayers(M).Filter
                                        if strcmpi(obj.CurrentPlayers(M).Kernel,'PreFiltering')
                                            obj.CurrentPlayers(M).Plot.CData = imfilter(obj.CurrentPlayers(M).Readers.GetFrame(IndxFrame,[obj.Parameters.PreFiltering.d1 obj.Parameters.PreFiltering.d2]),obj.Parameters.PreFiltering.Kernel,'symmetric');
                                        elseif strcmpi(obj.CurrentPlayers(M).Kernel,'Visualization')
                                            obj.CurrentPlayers(M).Plot.CData = imfilter(obj.CurrentPlayers(M).Readers.GetFrame(IndxFrame,[obj.Parameters.PreFiltering.d1 obj.Parameters.PreFiltering.d2]),obj.Parameters.Visualization.Kernel,'symmetric');
                                        elseif strcmpi(obj.CurrentPlayers(M).Kernel,'PostProcessing')
                                            obj.CurrentPlayers(M).Plot.CData = imfilter(obj.CurrentPlayers(M).Readers.GetFrame(IndxFrame,[obj.Parameters.PreFiltering.d1 obj.Parameters.PreFiltering.d2]),obj.Handles.Values.PostProcessing.Filtering.Kernel,'symmetric');
                                            if strcmpi(obj.Handles.UIElements.MovieDisplay(1).Plot.CLimMode,'auto')
                                                obj.Handles.UIElements.ThresholdLowPostProc_Edit.String = num2str(obj.CurrentPlayers(M).Axes.CLim(1));
                                                obj.Handles.UIElements.ThresholdHighPostProc_Edit.String = num2str(obj.CurrentPlayers(M).Axes.CLim(2));
                                                obj.Handles.Values.CaMovieFilteredCLim = obj.CurrentPlayers(M).Axes.CLim;
                                            end
                                        end
                                    else
                                        obj.CurrentPlayers(M).Plot.CData = obj.CurrentPlayers(M).Readers.GetFrame(IndxFrame, [obj.Parameters.PreFiltering.d1 obj.Parameters.PreFiltering.d2]);
                                        if strcmpi(obj.CurrentPlayers(M).Kernel,'PostProcessing')
                                            if strcmpi(obj.Handles.UIElements.MovieDisplay(1).Plot.CLimMode,'auto')
                                                obj.Handles.UIElements.ThresholdLowPostProc_Edit.String = num2str(obj.CurrentPlayers(M).Axes.CLim(1));
                                                obj.Handles.UIElements.ThresholdHighPostProc_Edit.String = num2str(obj.CurrentPlayers(M).Axes.CLim(2));
                                                obj.Handles.Values.CaMovieCLim = obj.CurrentPlayers(M).Axes.CLim;
                                            end
                                        end
                                    end
                                else
                                    Ended(M) = true;
                                end
                            otherwise
                                % Directly the file (e.g. behaviour)
                                % Should be .avi (previous UI), so we won't
                                % check here
                                obj.CurrentPlayers(M).Readers = VideoReader(obj.CurrentPlayers(M).Data);
                        end
                    end
                    obj.Handles.UIElements.MovieDisplay(1).Slider.XData = obj.CurrentTime * [1 1];
                    if isfield(obj.Handles.UIElements,'TimeLine')
                        obj.Handles.UIElements.TimeLine.XData = obj.CurrentTime * [1 1];
                    end
                    obj.Handles.UIElements.CurrentTimeEdit.String = num2str(1/100 * round(obj.CurrentTime*100),'%.2f');
                    drawnow
                    if obj.PlayingSingle
                        obj.PlayingSingle = false;
                    end
                end
                if all(Ended)
                    obj.Playing = false;
                    obj.Handles.UIElements.PlayButton.String = 'Play';
                    obj.CurrentTime = 0;
                    obj.Handles.UIElements.MovieDisplay(1).Slider.XData = obj.CurrentTime * [1 1];
                    obj.Handles.UIElements.CurrentTimeEdit.String = num2str(1/100 * round(obj.CurrentTime*100),'%.2f');
                end
            end
        end

        function CurrentTime_CB(obj)
            InputTime = str2double(obj.Handles.UIElements.CurrentTimeEdit.String);
            if ~isempty(InputTime) & InputTime <= obj.Handles.UIElements.MovieDisplay(1).SliderTickAxes.XLim(2),
                obj.CurrentTime = InputTime;
                obj.PlayingSingle = true;
                obj.PlayMovies;
            end
        end

        %% Processing functions

        % Preprocessing (Inscopix API)
        function PreProcess(obj)
            % Stop movies
            obj.Playing = false;
            pause(0.02)

            % Check whether we already have a preprocessing session with
            % the same parameters
            MatchingFolders = dir([obj.Parameters.PreProcessing.MainFolder 'PP-*']);


            % Special case: already preprocessed outside MATLAB
            % (e.g. for multiplane imaging)
            if contains(obj.Parameters.PreProcessing.Basename,'-PP')
                TrueBasename = strsplit(obj.Parameters.PreProcessing.Basename,'-PP');
                obj.Parameters.PreProcessing.Basename = TrueBasename{1};
                % Already preprocessed: we jump directly to the end of
                % preprocessing
                obj.Parameters.PreProcessing.ForkNumber = obj.GetFork('PreProcessing')+1;
                % Set the preprocessing files/folders names
                obj.Parameters.PreProcessing.PreProcessingFolder = ...
                    [obj.Parameters.PreProcessing.MainFolder 'PP-' num2str(obj.Parameters.PreProcessing.ForkNumber) filesep];
                mkdir(obj.Parameters.PreProcessing.PreProcessingFolder)
                obj.Parameters.PreProcessing.PreProcessingFile_isxd = ...
                    [obj.Parameters.PreProcessing.PreProcessingFolder obj.Parameters.PreProcessing.Basename '_PP-' num2str(obj.Parameters.PreProcessing.ForkNumber) '.isxd'];

                % Disable all interactions
                obj.DisableAll;
                obj.Parameters.PreProcessing.PreProcessingFile_isxd = obj.Parameters.PreProcessing.RawFile;
            else

                if ~isempty(MatchingFolders),
                    MatchingNoCropping = false;
                    MatchingAndCropping = false;
                    MatchingFolders = MatchingFolders([MatchingFolders.isdir]);
                    for M = 1 : numel(MatchingFolders)
                        MetaFiles = dir([MatchingFolders.folder '\*_StepMeta.mat']);
                        if ~isempty(MetaFiles),
                            MetaTemp = load(fullfile(MetaFiles.folder,MetaFiles.name));
                            if ...
                                    MetaTemp.PreProcessing.Spatial_Downsampling == obj.Parameters.PreProcessing.Spatial_Downsampling &&...
                                    MetaTemp.PreProcessing.Temporal_Downsampling == obj.Parameters.PreProcessing.Temporal_Downsampling &&...
                                    MetaTemp.PreProcessing.Fix_Pixels == obj.Parameters.PreProcessing.Fix_Pixels
                                if isempty(MetaTemp.PreProcessing.Cropping) && isempty(obj.Parameters.PreProcessing.Cropping)
                                    warndlg(['The file was already preprocessed with the same parameters. Please load the corresponding session instead of duplicating.'])
                                    return
                                else
                                    MatchingAndCropping = true;
                                end
                            end
                        end
                    end
                    if MatchingAndCropping
                        Answer = questdlg(['The file was already preprocessed with the same parameters, except for the cropping: is it worth reprocessing/duplicating everything?' newline newline,...
                            'What do you want to do?'],'Please choose...',...
                            'Reprocess again.',...
                            'Go back and load the already preprocessed branch.',...
                            'Go back and load the already preprocessed branch.');

                        if strcmpi(Answer, 'Go back and load the already preprocessed branch.'),
                            return
                        else
                            warndlg(['A new branch will be created. Consider deleting the previous one to save space.'])
                        end
                    end
                end

                obj.Parameters.PreProcessing.ForkNumber = obj.GetFork('PreProcessing')+1;
                % Set the preprocessing files/folders names
                obj.Parameters.PreProcessing.PreProcessingFolder = ...
                    [obj.Parameters.PreProcessing.MainFolder 'PP-' num2str(obj.Parameters.PreProcessing.ForkNumber) filesep];
                mkdir(obj.Parameters.PreProcessing.PreProcessingFolder)
                obj.Parameters.PreProcessing.PreProcessingFile_isxd = ...
                    [obj.Parameters.PreProcessing.PreProcessingFolder obj.Parameters.PreProcessing.Basename '_PP-' num2str(obj.Parameters.PreProcessing.ForkNumber) '.isxd'];

                % Disable all interactions
                obj.DisableAll;


                obj.Parameters.PreProcessing.PreProcessingFile = ...
                    [obj.Parameters.PreProcessing.PreProcessingFolder obj.Parameters.PreProcessing.Basename '_PP-' num2str(obj.Parameters.PreProcessing.ForkNumber) '.h5'];

                % Check which case (inscopix or not)
                if contains(obj.Parameters.PreProcessing.RawFile,'.isxd')
                    % Call the API
                    % We don't know how well the API handles "empty" cropping, so
                    % just in case, let's prepare two different calls
                    if isempty(obj.Parameters.PreProcessing.Cropping)
                        % No cropping
                        isx.preprocess(...
                            obj.Parameters.PreProcessing.RawFile,...
                            obj.Parameters.PreProcessing.PreProcessingFile_isxd,...                    )
                            'temporal_downsample_factor',obj.Parameters.PreProcessing.Temporal_Downsampling,...
                            'spatial_downsample_factor',obj.Parameters.PreProcessing.Spatial_Downsampling,...
                            'fix_defective_pixels', obj.Parameters.PreProcessing.Fix_Pixels,...
                            'trim_early_frames', true);
                    else
                        % Cropping
                        % Get the corners since that's what the API expects
                        % [top, left, bottom, right]
                        Rectangle = floor([...
                            obj.Parameters.PreProcessing.Cropping(2),...
                            obj.Parameters.PreProcessing.Cropping(1),...
                            obj.Parameters.PreProcessing.Cropping(2)+obj.Parameters.PreProcessing.Cropping(4),...
                            obj.Parameters.PreProcessing.Cropping(1)+obj.Parameters.PreProcessing.Cropping(3)]);
                        % Call
                        isx.preprocess(...
                            obj.Parameters.PreProcessing.RawFile,...
                            obj.Parameters.PreProcessing.PreProcessingFile_isxd,...
                            'temporal_downsample_factor',obj.Parameters.PreProcessing.Temporal_Downsampling,...
                            'spatial_downsample_factor',obj.Parameters.PreProcessing.Spatial_Downsampling,...
                            'crop_rect',Rectangle,...
                            'fix_defective_pixels', obj.Parameters.PreProcessing.Fix_Pixels,...
                            'trim_early_frames', true);
                    end
                    % Save as .h5 (for motion correction step) - we can afford to
                    % load everything before saving (at least these days)
                    RawMovie = isx.Movie.read(obj.Parameters.PreProcessing.PreProcessingFile_isxd);
                    % Waiting loop...
                    FileSize = 0;
                    while FileSize == 0
                        pause(0.05)
                        FileProc = dir(obj.Parameters.PreProcessing.PreProcessingFile_isxd);
                        FileSize = FileProc.bytes;
                    end
                    DownSampledFrames = NaN(RawMovie.spacing.num_pixels(1),RawMovie.spacing.num_pixels(2),RawMovie.timing.num_samples);
                    for F = 0 : RawMovie.timing.num_samples-1
                        DownSampledFrames(:,:,F+1) = RawMovie.get_frame_data(F);
                    end
                    DownSampledFrames = (DownSampledFrames); % CANNOT directly convert to uint8 (clipping -might be worth to just upscale it)
                    Times = RawMovie.timing.get_offsets_since_start;
                    DownSampled_Times = [Times.secs_float];
                    h5create(obj.Parameters.PreProcessing.PreProcessingFile,'/mov',size(DownSampledFrames),'Chunksize',[size(DownSampledFrames,[1 2]),10],'Datatype','double');
                    h5write(obj.Parameters.PreProcessing.PreProcessingFile,'/mov',DownSampledFrames);
                else
                    % Save as .h5 (for motion correction step) - we can afford to
                    % load everything before saving (at least these days)
                    DownSampledFrames = h5read(obj.Parameters.PreProcessing.RawFile,'/mov');

                    if ~isempty(obj.Parameters.PreProcessing.Cropping)
                        Rectangle = floor([...
                            obj.Parameters.PreProcessing.Cropping(2),...
                            obj.Parameters.PreProcessing.Cropping(1),...
                            obj.Parameters.PreProcessing.Cropping(2)+obj.Parameters.PreProcessing.Cropping(4),...
                            obj.Parameters.PreProcessing.Cropping(1)+obj.Parameters.PreProcessing.Cropping(3)]);
                        DownSampledFrames = DownSampledFrames(Rectangle(1):Rectangle(3),Rectangle(2):Rectangle(4),:);
                    end


                    if obj.Parameters.PreProcessing.Spatial_Downsampling>1 || obj.Parameters.PreProcessing.Temporal_Downsampling>1
                        % Spatial downsampling
                        if obj.Parameters.PreProcessing.Spatial_Downsampling>1
                            DownSize = floor(size(DownSampledFrames,[1 2])/obj.Parameters.PreProcessing.Spatial_Downsampling);
                            InitialHandled = DownSize*obj.Parameters.PreProcessing.Spatial_Downsampling;
                            Spatial_DownSampled = (NaN([DownSize size(DownSampledFrames,3)]));
                            averagingFunction = @(x) nanmean(x.data(:));
                            Spatial_Downsampling = obj.Parameters.PreProcessing.Spatial_Downsampling;
                            parfor F = 1 : size(Spatial_DownSampled,3)
                                Spatial_DownSampled(:,:,F) = (blockproc(single(DownSampledFrames(1:InitialHandled(1),1:InitialHandled(2),F)), Spatial_Downsampling * [1,1], averagingFunction));
                            end
                        else
                            Spatial_DownSampled = single(DownSampledFrames);
                        end
                        clear DownSampledFrames

                        % Temporal downsampling
                        if obj.Parameters.PreProcessing.Temporal_Downsampling>1
                            DownSize = floor(size(Spatial_DownSampled,3)/obj.Parameters.PreProcessing.Temporal_Downsampling);
                            Temporal_DownSampled = (NaN([size(Spatial_DownSampled,[1 2]),DownSize]));
                            Temporal_DownSampled_Times = NaN(DownSize,1);
                            parfor F = 1 : DownSize
                                IndxF = ((F-1)*obj.Parameters.PreProcessing.Temporal_Downsampling)+(1:obj.Parameters.PreProcessing.Temporal_Downsampling);
                                Temporal_DownSampled(:,:,F) = (nanmean(Spatial_DownSampled(:,:,IndxF)));
                                DownSampled_Times(F,1) = nanmean(obj.Parameters.PreProcessing.TimeStamps(IndxF));
                            end
                            Temporal_DownSampled = uint8(Temporal_DownSampled);
                        else
                            Temporal_DownSampled = uint8(Spatial_DownSampled);
                            DownSampled_Times = obj.Parameters.PreProcessing.TimeStamps;
                        end
                        delete(Spatial_DownSampled)
                    else
                        DownSampled_Times = obj.Parameters.PreProcessing.TimeStamps;
                        Temporal_DownSampled = DownSampledFrames;
                    end
                    h5create(obj.Parameters.PreProcessing.PreProcessingFile,'/mov',size(Temporal_DownSampled),'Chunksize',[size(Temporal_DownSampled,[1 2]),10],'Datatype','uint8');
                    h5write(obj.Parameters.PreProcessing.PreProcessingFile,'/mov',Temporal_DownSampled);
                    clear Temporal_DownSampled
                end
            end
            h5create(obj.Parameters.PreProcessing.PreProcessingFile,'/times',size(DownSampled_Times),'Datatype','double');
            h5write(obj.Parameters.PreProcessing.PreProcessingFile,'/times',DownSampled_Times);
            obj.Parameters.PreProcessing.Date = datetime;

            % Save Meta file
            MetaFile = [obj.Parameters.PreProcessing.PreProcessingFolder obj.Parameters.PreProcessing.Basename '_PP-' num2str(obj.Parameters.PreProcessing.ForkNumber) '_StepMeta.mat'];
            PreProcessing = obj.Parameters.PreProcessing;
            save(MetaFile,'PreProcessing')

            % Enable all interactions
            obj.EnableAll;
            drawnow

            % Switch to second step
            obj.PreFiltering_Layout;

            % Update the tab display
            Tabs = fieldnames(obj.Handles.Tabs);
            for F = 1 : numel(Tabs)
                obj.Handles.Tabs.(Tabs{F}).UserData = '0';
            end
            obj.Handles.Tabs.PreFiltering.UserData = '1';
            obj.CurrentTab = 'PreFiltering';
            obj.TabAppearence;
        end

        % Prefiltering
        function PreFilter(obj)
            % Stop movies
            obj.Playing = false;
            pause(0.02)
            % Check whether we already have a prefiltering session with
            % the same parameters
            MatchingFolders = dir([obj.Parameters.PreProcessing.PreProcessingFolder 'PP-*']);
            if ~isempty(MatchingFolders),
                MatchingNoCropping = false;
                MatchingAndCropping = false;
                MatchingFolders = MatchingFolders([MatchingFolders.isdir]);
                for M = 1 : numel(MatchingFolders)
                    MetaFiles = dir([MatchingFolders.folder '\**\*_StepMeta.mat']);
                    if ~isempty(MetaFiles),
                        MetaTemp = load(fullfile(MetaFiles,folder,MetaFiles.name));
                        if ...
                                MetaTemp.PreFiltering.Sig1 == obj.Parameters.PreFiltering.Sig1 &&...
                                MetaTemp.PreFiltering.Sig2 == obj.Parameters.PreFiltering.Sig2 &&...
                                MetaTemp.PreFiltering.WindowSize == obj.Parameters.PreFiltering.WindowSize,
                            if isempty(MetaTemp.PreFiltering.Cropping) && isempty(obj.Parameters.PreFiltering.Cropping)
                                warndlg(['The file was already prefiltered with the same parameters. Please load the corresponding session instead of duplicating.'])
                                return
                            else
                                MatchingAndCropping = true;
                            end
                        end
                    end
                end
                if MatchingAndCropping
                    Answer = questdlg(['The file was already prefiltered with the same parameters, except for the cropping: is it worth reprocessing/duplicating everything?' newline newline,...
                        'What do you want to do?'],'Please choose...',...
                        'Reprocess again.',...
                        'Go back and load the already prefiltered branch.',...
                        'Go back and load the already prefiltered branch.');

                    if strcmpi(Answer, 'Go back and load the already prefiltered branch.'),
                        return
                    else
                        warndlg(['A new branch will be created. Consider deleting the previous one to save space.'])
                    end
                end
            end
            obj.Parameters.PreFiltering.ForkNumber = obj.GetFork('PreFiltering')+1;
            % Set the preprocessing files/folders names
            obj.Parameters.PreFiltering.PreFilteringFolder = ...
                [obj.Parameters.PreProcessing.PreProcessingFolder 'PF-' num2str(obj.Parameters.PreFiltering.ForkNumber) filesep];
            mkdir(obj.Parameters.PreFiltering.PreFilteringFolder)
            obj.Parameters.PreFiltering.PreFilteringFile = ...
                [obj.Parameters.PreFiltering.PreFilteringFolder obj.Parameters.PreProcessing.Basename '_PP-' num2str(obj.Parameters.PreProcessing.ForkNumber) '_PF-' num2str(obj.Parameters.PreFiltering.ForkNumber) '.h5'];

            % Disable all interactions
            obj.DisableAll;

            % Update status display

            % Prepare cropping mask if needed
            if isfield(obj.Handles.UIElements,'CropArea')
                if ~isempty(obj.Handles.UIElements.CropArea)
                    CroppingMask = obj.Handles.UIElements.CropArea(2).createMask(obj.Parameters.PreFiltering.d1,obj.Parameters.PreFiltering.d2);
                else
                    CroppingMask = true(obj.Parameters.PreFiltering.d1,obj.Parameters.PreFiltering.d2);
                end
            else
                CroppingMask = true(obj.Parameters.PreFiltering.d1,obj.Parameters.PreFiltering.d2);
            end

            % Apply filter and cropping to each frame (we have enough RAM,
            % so everything will be loaded/filtered and then only written)
            Filtered = single(zeros(obj.Parameters.PreFiltering.d1,obj.Parameters.PreFiltering.d2,obj.Parameters.PreFiltering.FramesNum));
            % parpool cannot share the same "reader" from the API, but
            % creating one per loop still saves a bit of time when compared
            % to a normal for loop although some time is spent retrieving
            % data

            % Check if the parallel pool is running, otherwise display a
            % message to inform it is starting (takes a few seconds)

            % Run the loop
            tic
            FileR = obj.Parameters.PreProcessing.PreProcessingFile;
            Kernel = obj.Parameters.PreFiltering.Kernel;
            Threshold =  obj.Parameters.PreFiltering.Threshold;
            PreProcessed = h5read(FileR,'/mov');
            parfor F = 1 : obj.Parameters.PreFiltering.FramesNum
                FilteredTemp = single(PreProcessed(:,:,F));
                FilteredTemp = imfilter(FilteredTemp,Kernel,'symmetric');
                if isempty(Threshold)
                    FilteredTemp(~CroppingMask) = 0;
                else
                    FilteredTemp(~CroppingMask & FilteredTemp<Threshold(1)) = 0;
                    FilteredTemp(FilteredTemp>Threshold(2)) = Threshold(2);
                end
                Filtered(:,:,F) = FilteredTemp;
            end
            toc
            % Prepare .h5 file
            % The optimal chunksize has not been tested for this pipeline
            % and the usual size of our recordings, but we try to take into
            % account the limitations mentioned for .h5: smaller chunksize
            % means that finding one among all the other is slow; bigger
            % chunksize means that we potentially have to read a lot of
            % data to reach the value we want to read, which is also a bit
            % slow. Also, there is a maximum chunksize (4GB)
            WhosFiltered = whos('Filtered');

            %             ChunkFrames = floor(0.1*4e9/(WhosFiltered.bytes/obj.Parameters.PreFiltering.FramesNum));
            %             if ChunkFrames>obj.Parameters.PreFiltering.FramesNum
            %                 ChunkFrames = obj.Parameters.PreFiltering.FramesNum;
            %             end
            h5create(obj.Parameters.PreFiltering.PreFilteringFile,'/mov',[obj.Parameters.PreFiltering.d1,obj.Parameters.PreFiltering.d2,obj.Parameters.PreFiltering.FramesNum],'Chunksize',[obj.Parameters.PreFiltering.d1,obj.Parameters.PreFiltering.d2,10],'Datatype','uint8');
            h5write(obj.Parameters.PreFiltering.PreFilteringFile,'/mov',Filtered);
            obj.Parameters.PreFiltering.Date = datetime;

            % Save Meta file (append just the current step, to prevent
            % any change to the fields already processed if we went
            % back to a previous tab and changed stuff only to watch
            % the movies)
            MetaFilePrevious = [obj.Parameters.PreProcessing.PreProcessingFolder obj.Parameters.PreProcessing.Basename '_PP-' num2str(obj.Parameters.PreProcessing.ForkNumber) '_StepMeta.mat'];
            PreviousMeta = load(MetaFilePrevious);
            MetaFileCurrent = [obj.Parameters.PreFiltering.PreFilteringFolder obj.Parameters.PreProcessing.Basename '_PP-' num2str(obj.Parameters.PreProcessing.ForkNumber) '_PF-' num2str(obj.Parameters.PreFiltering.ForkNumber) '_StepMeta.mat'];
            PreviousMeta.PreFiltering = obj.Parameters.PreFiltering;
            save(MetaFileCurrent,'-struct','PreviousMeta')

            % Enable all interactions
            obj.EnableAll;
            drawnow

            % Switch to third step
            obj.MotionCorrection_Layout;

            % Update the tab display
            Tabs = fieldnames(obj.Handles.Tabs);
            for F = 1 : numel(Tabs),
                obj.Handles.Tabs.(Tabs{F}).UserData = '0';
            end
            obj.Handles.Tabs.MotionCorrection.UserData = '1';
            obj.CurrentTab = 'MotionCorrection';
            obj.TabAppearence;

        end

        % Motion correction
        function MotionCorrect(obj)
            % Disable interactions
            obj.DisableAll;

            obj.Parameters.MotionCorrection.ForkNumber = obj.GetFork('MotionCorrection')+1;

            % Non-rigid motion correction
            options_nonrigid = NoRMCorreSetParms(...
                'd1',obj.Parameters.PreFiltering.d1,...
                'd2',obj.Parameters.PreFiltering.d2,...
                'grid_size', obj.Parameters.MotionCorrection.grid_size,...
                'min_patch_size', ceil([obj.Parameters.MotionCorrection.grid_size([1 2]) * 0.5 1]),...
                'overlap_pre', obj.Parameters.MotionCorrection.overlap_pre,...
                'overlap_post', obj.Parameters.MotionCorrection.overlap_post,...
                'max_dev',obj.Parameters.MotionCorrection.max_dev,...
                'output_type','mat',...
                'correct_bidir',false,...
                'iter',obj.Parameters.MotionCorrection.iter,...
                'max_shift',obj.Parameters.MotionCorrection.max_shift,...
                'init_batch',obj.Parameters.MotionCorrection.init_batch,...
                'boundary','zero',...
                'bin_width',obj.Parameters.MotionCorrection.bin_width);

            % Find shifts on filtered data
            tic; [~,Shifts,~] = normcorre_batch(obj.Parameters.PreFiltering.PreFilteringFile,options_nonrigid); toc


            % The motion corrected filtered movie is not saved -to check
            % the quality of the motion correction, filtering the motion
            % corrected raw movie yields the same with barely any processing
            % cost, saving disk space


            % See notes in PreFiltering regarding .h5's chunksize
            FileDir = dir(obj.Parameters.PreProcessing.RawFile);
            SingleFrameSize = FileDir.bytes / obj.Parameters.PreFiltering.FramesNum;
            Memory_Batch_Size = floor(0.1 * floor(4e9 / SingleFrameSize));

            % Set the motion correction files/folders names
            obj.Parameters.MotionCorrection.MotionCorrectionFolder = ...
                [obj.Parameters.PreFiltering.PreFilteringFolder 'MC-' num2str(obj.Parameters.MotionCorrection.ForkNumber) filesep];
            mkdir( obj.Parameters.MotionCorrection.MotionCorrectionFolder)
            obj.Parameters.MotionCorrection.MotionCorrectionFile = ...
                [ obj.Parameters.MotionCorrection.MotionCorrectionFolder,...
                obj.Parameters.PreProcessing.Basename,...
                '_PP-' num2str(obj.Parameters.PreProcessing.ForkNumber),...
                '_PF-' num2str(obj.Parameters.PreFiltering.ForkNumber),...
                '_MC-' num2str(obj.Parameters.MotionCorrection.ForkNumber),...
                '.h5'];

            %            RawMovieFrames = h5read(obj.Parameters.PreProcessing.PreProcessingFile,'/mov');

            % Apply to the raw movie
            options_nonrigid.output_type = 'hdf5';
            options_nonrigid.mem_batch_size = Memory_Batch_Size;
            options_nonrigid.h5_filename = obj.Parameters.MotionCorrection.MotionCorrectionFile;
            tic; apply_shifts(obj.Parameters.PreProcessing.PreProcessingFile,Shifts,options_nonrigid); toc % apply the shifts to full dataset
            obj.Parameters.MotionCorrection.Date = datetime;

            % Save Meta file (append just the current step, to prevent
            % any change to the fields already processed if we went
            % back to a previous tab and changed stuff only to watch
            % the movies)
            MetaFilePrevious = [obj.Parameters.PreFiltering.PreFilteringFolder obj.Parameters.PreProcessing.Basename '_PP-' num2str(obj.Parameters.PreProcessing.ForkNumber) '_PF-' num2str(obj.Parameters.PreFiltering.ForkNumber) '_StepMeta.mat'];
            PreviousMeta = load(MetaFilePrevious);
            MetaFileCurrent = [...
                obj.Parameters.MotionCorrection.MotionCorrectionFolder,...
                obj.Parameters.PreProcessing.Basename,...
                '_PP-' num2str(obj.Parameters.PreProcessing.ForkNumber),...
                '_PF-' num2str(obj.Parameters.PreFiltering.ForkNumber),...
                '_MC-' num2str(obj.Parameters.MotionCorrection.ForkNumber) '_StepMeta.mat'];
            PreviousMeta.MotionCorrection = obj.Parameters.MotionCorrection;
            save(MetaFileCurrent,'-struct','PreviousMeta')


            % Enable all interactions
            obj.EnableAll;
            drawnow

            % Switch to CNMF-E step
            obj.CNMFE_Layout;

            % Update the tab display
            Tabs = fieldnames(obj.Handles.Tabs);
            for F = 1 : numel(Tabs),
                obj.Handles.Tabs.(Tabs{F}).UserData = '0';
            end
            obj.Handles.Tabs.CNMFE.UserData = '1';
            obj.CurrentTab = 'CNMFE';
            obj.TabAppearence;

        end

        % CNMF-E
        function CNMFE(obj)

            % CNMFE outputs the results to the workspace with an eval
            % somewhere... Couldn't find it, but it takes a huge piece of
            % RAM, so we can compare the variables before/after to clean
            % that
            PreVar = evalin('base','who');

            obj.Parameters.CNMFE.ForkNumber = obj.GetFork('CNMFE')+1;
            obj.Parameters.CNMFE.CNMFEFolder = ...
                [obj.Parameters.MotionCorrection.MotionCorrectionFolder 'CNMFE-' num2str(obj.Parameters.CNMFE.ForkNumber) filesep];
            mkdir(obj.Parameters.CNMFE.CNMFEFolder)

            % Crop data if needed
            if ~isempty(obj.Parameters.CNMFE.Crop)
                disp('Applying cropping before running CNMFE...')
                % Retrieve frames
                MovieRead = h5read(obj.Parameters.MotionCorrection.MotionCorrectionFile,'/mov');

                % For inscopix, convert to uint8
                if ~isinteger(MovieRead)
                    MAX = max(MovieRead,[],'all');
                    MIN = min(MovieRead,[],'all');
                    MovieRead = (MovieRead-MIN)/(MAX-MIN);
                    MovieRead = uint8(round(255*MovieRead));    
                end
                obj.Handles.UIElements.Crop.Position = round(obj.Handles.UIElements.Crop.Position);
                MovieRead = MovieRead(obj.Handles.UIElements.Crop.Position(2):obj.Handles.UIElements.Crop.Position(2)+obj.Handles.UIElements.Crop.Position(4),...
                    obj.Handles.UIElements.Crop.Position(1):obj.Handles.UIElements.Crop.Position(1)+obj.Handles.UIElements.Crop.Position(3),:);


                %                obj.Parameters.CNMFE.DataFile = ...
                %                    [obj.Parameters.CNMFE.CNMFEFolder obj.Parameters.PreProcessing.Basename '_PP-' num2str(obj.Parameters.PreProcessing.ForkNumber) '_PF-' num2str(obj.Parameters.PreFiltering.ForkNumber) '_MC-' num2str(obj.Parameters.MotionCorrection.ForkNumber) '.h5'];

                obj.Parameters.CNMFE.DataFile = ...
                    [obj.Parameters.CNMFE.CNMFEFolder obj.Parameters.PreProcessing.Basename '.h5'];


                % Write as a new .h5
                %                 WhosFiltered = whos('MovieRead');
                %                 ChunkFrames = floor(0.1*4e9/(WhosFiltered.bytes/obj.Parameters.PreFiltering.FramesNum));
                %                 if ChunkFrames > obj.Parameters.PreFiltering.FramesNum,
                %                     ChunkFrames = obj.Parameters.PreFiltering.FramesNum;
                %                 end

%                                 h5create(obj.Parameters.CNMFE.DataFile,'/mov',[size(MovieRead,1),size(MovieRead,2),obj.Parameters.PreFiltering.FramesNum],'Chunksize',[size(MovieRead,1),size(MovieRead,2),10],'Datatype','uint8');


                                
                h5create(obj.Parameters.CNMFE.DataFile,'/mov',[size(MovieRead,1),size(MovieRead,2),obj.Parameters.PreFiltering.FramesNum],'Chunksize',[size(MovieRead,1),size(MovieRead,2),10],'Datatype','uint8');
                h5write(obj.Parameters.CNMFE.DataFile,'/mov',MovieRead);
                disp('Done!')
            else
                obj.Parameters.CNMFE.DataFile = obj.Parameters.MotionCorrection.MotionCorrectionFile;
                MovieRead = h5read(obj.Parameters.CNMFE.DataFile,'/mov');
            end


            H5IC = h5info(obj.Parameters.CNMFE.DataFile);
            CNMFE_MovieSize = H5IC.Datasets.Dataspace.Size;

            obj.Parameters.CNMFE.Processed = false;
            % Some parameters have to be adjusted (depend on others)
            obj.Parameters.CNMFE.Neuron.options.min_pixel = 2*pi*(0.5*obj.Parameters.CNMFE.gSig)^2;  % minimum number of nonzero pixels for each neuron
            obj.Parameters.CNMFE.Neuron.options.dmin = obj.Parameters.CNMFE.gSig*1.25; % minimum distances between two neurons. It is used together with merge_thr
            obj.Parameters.CNMFE.Neuron.options.deconv_options.max_tau = ceil(obj.Parameters.CNMFE.Neuron.options.deconv_options.max_tau*obj.Parameters.PreFiltering.SamplingRate);
            % -------------------------    COMPUTATION    -------------------------  %
            
            % Adjusting the patch size for larger movies/FOV because the
            % noise-estimation step engulfs lots of RAM (handled by Matlab,
            % not adjusted by CNMFE besides patching)

            if 0
            PixelsNum = CNMFE_MovieSize(1) * CNMFE_MovieSize(2) * CNMFE_MovieSize(3);
%             MovieSize_GB_Double = PixelsNum*(64/8)/(2^30);
            WorkablePixelsNum = 0.5*1.5e9;
            if PixelsNum>WorkablePixelsNum
                PatchFactor = round(PixelsNum/WorkablePixelsNum);
                PatchSize = round(CNMFE_MovieSize([1,2])/PatchFactor.^0.5);
            else
                PatchSize = [1 1] * round(128 / obj.Parameters.PreProcessing.Spatial_Downsampling);
            end
            end
            PatchSize = [1 1] * round(256 / obj.Parameters.PreProcessing.Spatial_Downsampling);
            MemoryToUse = obj.Parameters.CNMFE.maxRAM; % GB
            PatchSize_GB_Double = PatchSize(1)*PatchSize(2)* CNMFE_MovieSize(3)*(64/8)/(2^30);

            if PatchSize_GB_Double>MemoryToUse  
                PatchFactor = round(PatchSize_GB_Double/MemoryToUse);
                PatchSize = round(CNMFE_MovieSize([1,2])/(PatchFactor));
                PatchSize_GB_Double = PatchSize(1)*PatchSize(2)* CNMFE_MovieSize(3)*(64/8)/(2^30);
            end

            pars_envs = struct('memory_size_to_use',MemoryToUse, ... % GB, memory space you allow to use in MATLAB
                'memory_size_per_patch', 2, ...             % GB, space for loading data within one patch
                'patch_dims', PatchSize);       % px, patch size
            show_merge = false;
            use_parallel = true;                             % use parallel computation for parallel computing

            % Process the maximum intensity projection with the chosen
            % parameters, so that it is saved for further use
            gSig = obj.Parameters.CNMFE.gSig;
            % Create Kernel
            psf = fspecial('gaussian', ceil(gSig*4+1), gSig);
            ind_nonzero = (psf(:)>=max(psf(:,1)));
            psf = psf-mean(psf(ind_nonzero));
            psf(~ind_nonzero) = 0;

            % Filter with the same kind of filter as CNMFE
            Filtered = single(zeros(CNMFE_MovieSize));

            tic
            parfor FM = 1 : size(Filtered,3)
                FrameFM = single(MovieRead(:,:,FM));
                Filtered(:,:,FM) = imfilter(FrameFM,psf,'replicate');
                %                 Filtered(:,:,FM) = imfilter(single(h5read(DataFile,'/mov',[1,1,FM],[CNMFE_MovieSize(1),CNMFE_MovieSize(2),1])),psf,'replicate');
            end
            toc
if 0
            % Determine the proper number of workers, becomes crucial
            % for huge recordings(otherwise, will engage all cores and 
            % saturate RAM anyway, the way CNMFE implemented it); 
            % CNMFE's code derives several variables from the raw data per 
            % patch, which blows in volume even more: it needs to be taken 
            % into account
            % NB: it's not linearly correlated to movie size 
            % (so trying random stuff to decrease number of workers in a
            % somewhat related manner...)
            MaxNumberWorkers = floor(MemoryToUse/(2*PatchSize_GB_Double));
            try
                if MaxNumberWorkers == 0 || MaxNumberWorkers>64
                    delete(gcp);
                else
                    delete(gcp('nocreate'));
                    parpool(MaxNumberWorkers)
                end
            catch
                delete(gcp);
            end
            %%%%% End of the uggly bit %%%%%
end

            MP = max(Filtered,[],3);
            % Normalize for easier CLim adjustment
            MP = (MP - min(MP,[],'all')) / (max(MP,[],'all') - min(MP,[],'all'));
            obj.Parameters.CNMFE.MaxProjection = MP;

            % Reclaim memory
            clear MovieRead Filtered

            % Distribute data and get ready to run source extraction
            obj.Parameters.CNMFE.Neuron.select_data(obj.Parameters.CNMFE.DataFile);
            obj.Parameters.CNMFE.Neuron.getReady(pars_envs);

            % Initialize neurons from the video data within a selected temporal range
            K = [];                      % maximum number of neurons per patch. when K=[], take as many as possible.
            frame_range = [];            % when [], uses all frames
            save_initialization = false; % save the initialization procedure as a video.
            tic
            try
                use_parallelSpatial = true;
                [center, Cn, PNR] = obj.Parameters.CNMFE.Neuron.initComponents_parallel(K, frame_range, save_initialization, use_parallel);
            catch
                use_parallelSpatial = false;
                [center, Cn, PNR] = obj.Parameters.CNMFE.Neuron.initComponents_parallel(K, frame_range, save_initialization, false);
            end
            toc
            obj.Parameters.CNMFE.Neuron.compactSpatial;
            figure;
            imagesc(Cn.*PNR, [0, 0.98*max(PNR,[],'all')]); colormap bone;
            hold on;
            plot(center(:, 2), center(:, 1), '.r', 'markersize', 10);
            drawnow;
            hold off

            % Estimate the background components
            obj.Parameters.CNMFE.Neuron.update_background_parallel(use_parallel);

            %  merge neurons and update spatial/temporal components
            obj.Parameters.CNMFE.Neuron.merge_neurons_dist_corr(show_merge);
            obj.Parameters.CNMFE.Neuron.merge_high_corr(show_merge, obj.Parameters.CNMFE.merge_thr_spatial);

            % pick neurons from the residual
            if 0
                min_corr_res = obj.Parameters.CNMFE.Neuron.options.min_corr;
                min_pnr_res = obj.Parameters.CNMFE.Neuron.options.min_pnr-0.5;
                seed_method_res = 'auto';  % method for initializing neurons from the residual
                obj.Parameters.CNMFE.Neuron.initComponents_residual_parallel([], save_initialization, use_parallel, min_corr_res, min_pnr_res, seed_method_res);
            end

            % update spatial
            obj.Parameters.CNMFE.Neuron.update_spatial_parallel(use_parallelSpatial, true);

            % merge neurons based on correlations
            obj.Parameters.CNMFE.Neuron.merge_high_corr(show_merge, obj.Parameters.CNMFE.merge_thr_spatial);

            for m=1:2
                % update temporal
                obj.Parameters.CNMFE.Neuron.update_temporal_parallel(use_parallelSpatial);

                % delete bad neurons
                obj.Parameters.CNMFE.Neuron.remove_false_positives();

                % merge neurons based on temporal correlation + distances
                obj.Parameters.CNMFE.Neuron.merge_neurons_dist_corr(show_merge);
            end


            % For our signals, keeping the hals algorithm appareas to work better than
            % switching to using the nonnegative least square method as a solver.
            %     obj.Parameters.CNMFE.Neuron.options.spatial_algorithm = 'nnls';

            % Run more iterations
            obj.Parameters.CNMFE.Neuron.update_background_parallel(use_parallel);
            obj.Parameters.CNMFE.Neuron.update_spatial_parallel(use_parallelSpatial);
            obj.Parameters.CNMFE.Neuron.update_temporal_parallel(use_parallelSpatial);

            K = size(obj.Parameters.CNMFE.Neuron.A,2);
            tags = obj.Parameters.CNMFE.Neuron.tag_neurons_parallel();  % find neurons with fewer nonzero pixels than min_pixel and silent calcium transients
            obj.Parameters.CNMFE.Neuron.remove_false_positives();
            obj.Parameters.CNMFE.Neuron.merge_neurons_dist_corr(show_merge);
            obj.Parameters.CNMFE.Neuron.merge_high_corr(show_merge, obj.Parameters.CNMFE.merge_thr_spatial);

            if K~=size(obj.Parameters.CNMFE.Neuron.A,2)
                obj.Parameters.CNMFE.Neuron.update_spatial_parallel(use_parallelSpatial);
                obj.Parameters.CNMFE.Neuron.update_temporal_parallel(use_parallelSpatial);
                obj.Parameters.CNMFE.Neuron.remove_false_positives();
            end
            obj.Parameters.CNMFE.VanillaNeuron = obj.Parameters.CNMFE.Neuron;
            obj.Parameters.CNMFE.Processed = true;
            MetaFilePrevious = [obj.Parameters.MotionCorrection.MotionCorrectionFolder obj.Parameters.PreProcessing.Basename '_PP-' num2str(obj.Parameters.PreProcessing.ForkNumber) '_PF-' num2str(obj.Parameters.PreFiltering.ForkNumber) '_MC-' num2str(obj.Parameters.MotionCorrection.ForkNumber) '_StepMeta.mat'];
            PreviousMeta = load(MetaFilePrevious);
            MetaFileCurrent = [...
                obj.Parameters.CNMFE.CNMFEFolder,...
                obj.Parameters.PreProcessing.Basename,...
                '_PP-' num2str(obj.Parameters.PreProcessing.ForkNumber),...
                '_PF-' num2str(obj.Parameters.PreFiltering.ForkNumber),...
                '_MC-' num2str(obj.Parameters.MotionCorrection.ForkNumber),...
                '_CNMFE-' num2str(obj.Parameters.CNMFE.ForkNumber) '_StepMeta.mat'];
            PreviousMeta.CNMFE = obj.Parameters.CNMFE;
            save(MetaFileCurrent,'-struct','PreviousMeta')

            drawnow
            % Enable all interactions
            obj.EnableAll;
            drawnow

            % Switch to CNMF-E step
            obj.PostProcessing_Layout(true);

            % Update the tab display
            Tabs = fieldnames(obj.Handles.Tabs);
            for F = 1 : numel(Tabs)
                obj.Handles.Tabs.(Tabs{F}).UserData = '0';
            end
            obj.Handles.Tabs.PostProcessing.UserData = '1';
            obj.CurrentTab = 'PostProcessing';
            obj.TabAppearence;

            PostVar = evalin('base','who');
            NewVar = PostVar(arrayfun(@(x) ~any(strcmpi(PreVar,PostVar{x})),1:numel(PostVar)));
            arrayfun(@(x) evalin('base',['clear ' NewVar{x}]),1:numel(NewVar))

            % Delete custom parpool (could decrease performances if used
            % for other steps as is)
            delete(gcp('nocreate'));

        end

        % Determine forks numbers
        function CurrentFork = GetFork(obj,Level)
            % Retrieve matching folders
            switch Level
                case 'PreProcessing'
                    MatchingFolders = dir([obj.Parameters.PreProcessing.MainFolder 'PP-*']);
                case 'PreFiltering'
                    MatchingFolders = dir([obj.Parameters.PreProcessing.PreProcessingFolder 'PF-*']);
                case 'MotionCorrection'
                    MatchingFolders = dir([obj.Parameters.PreFiltering.PreFilteringFolder 'MC-*']);
                case 'CNMFE'
                    MatchingFolders = dir([obj.Parameters.MotionCorrection.MotionCorrectionFolder 'CNMFE-*']);
            end
            if ~isempty(MatchingFolders),
                MatchingFolders = MatchingFolders([MatchingFolders.isdir]);

                % Retrieve highest indice
                Indices = arrayfun(@(x) strsplit(MatchingFolders(x).name,'-'), 1:numel(MatchingFolders),'UniformOutput',false);
                Indices = arrayfun(@(x) str2double(Indices{x}{2}), 1:numel(MatchingFolders));
                CurrentFork = max(Indices);
            else
                CurrentFork = 0;
            end
        end

        %% Callbacks
        function TabCallback(obj,src,evt)
            if str2double(src.UserData) == 1,
                % Do nothing, already the current tab
                return
            else
                obj.Playing = false;
                drawnow
                % Check that we can switch
                Str = strrep(src.String,'-','');
                if strcmpi(Str,'PostProcessing') && obj.JustLoaded
                    if isfield(obj.Parameters,'CNMFE')
                        if isfield(obj.Parameters.CNMFE,'Neuron')
                            obj.Parameters.CNMFE.Processed = 1;
                        end
                    end
                    Status = obj.([Str '_Layout'])(true);
                    obj.JustLoaded = false;
                else
                    Status = obj.([Str '_Layout']);
                end
                if Status
                    % Update the tab display
                    Tabs = fieldnames(obj.Handles.Tabs);
                    for F = 1 : numel(Tabs),
                        obj.Handles.Tabs.(Tabs{F}).UserData = '0';
                    end
                    src.UserData = '1';
                    obj.CurrentTab = src.String;
                    obj.TabAppearence;
                end
            end
        end

        function TabAppearence(obj)
            Tabs = fieldnames(obj.Handles.Tabs);
            for F = 1 : numel(Tabs),
                JavaObj = findjobj(obj.Handles.Tabs.(Tabs{F}));
                if str2double(obj.Handles.Tabs.(Tabs{F}).UserData) == 1,
                    obj.Handles.Tabs.(Tabs{F}).ForegroundColor = [0.6 0.3 0];
                    obj.Handles.Tabs.(Tabs{F}).BackgroundColor = [0.2 0.2 0.2];
                    JavaObj.Border = [];
                    JavaObj.repaint;
                else
                    obj.Handles.Tabs.(Tabs{F}).ForegroundColor = 'k';
                    obj.Handles.Tabs.(Tabs{F}).BackgroundColor = [0.25 0.25 0.25];
                    BorderColor = java.awt.Color(0.1,0.1,0.1);
                    newBorder = javax.swing.border.LineBorder(BorderColor,3,true);
                    JavaObj.Border = newBorder;
                    JavaObj.repaint;
                end
            end
        end

        % Tabs layouts
        function WipeLayout(obj)
            if isfield(obj.Handles,'UIElements')
                % Retrieve the main categories
                UIElements = fieldnames(obj.Handles.UIElements);
                for U = 1 : numel(UIElements)
                    for S = 1 : numel(obj.Handles.UIElements.(UIElements{U}))
                        if isstruct(obj.Handles.UIElements.(UIElements{U})(S))
                            % Retrieve the subfields
                            UIElements_U = fieldnames(obj.Handles.UIElements.(UIElements{U})(S));
                            Check_U = obj.Handles.UIElements.(UIElements{U})(S);
                            Check_U = Check_U(arrayfun(@(x) isstruct(Check_U(x)),1:numel(Check_U)));
                            % Check which are either uicontrol or axes (not the plots)
                            ToDelete = find(arrayfun(@(x)...
                                strcmpi(class(Check_U.(UIElements_U{x})),'matlab.graphics.axis.Axes')|...
                                strcmpi(class(Check_U.(UIElements_U{x})),'matlab.ui.control.UIControl') ...
                                ,1:numel(UIElements_U)));
                            for TD = 1 : numel(ToDelete)
                                % arrayfun(@(x) delete(obj.Handles.UIElements.(UIElements{U}).(UIElements_U{ToDelete(x)})), 1:numel(ToDelete));
                                Ax = obj.Handles.UIElements.(UIElements{U})(S).(UIElements_U{ToDelete(TD)});
                                delete(Ax)
                            end
                        else
                            if   strcmpi(class(obj.Handles.UIElements.(UIElements{U})),'matlab.graphics.axis.Axes')||...
                                    strcmpi(class(obj.Handles.UIElements.(UIElements{U})),'matlab.ui.control.UIControl'),
                                delete(obj.Handles.UIElements.(UIElements{U}));
                            end
                        end
                    end
                end
                % Reset the fields
                obj.Handles = rmfield(obj.Handles,'UIElements');
            end
        end

        %% Preprocessing tab & functions

        function Status = PreProcessing_Layout(obj)
            Status = false;
            if ~strcmpi(obj.CurrentTab,'PreProcessing'),
                Status = true;
                % Wipe current layout
                obj.WipeLayout;

                % Main display and uielements
                obj.SetMovieDisplay1;

                % Other uielements
                % Movie
                obj.Handles.UIElements.MovieInfos_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.6 0.87 0.25 0.03],...
                    'String','Raw calcium-imaging movie',...
                    'FontSize',16,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.7 0.7 0.7],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.LoadFile_Button = uicontrol('Style','pushbutton',...
                    'Units','Normalized','Position',[0.6 0.84 0.1 0.03],...
                    'String','Load file',...
                    'FontSize',16,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.7 0.7 0.7],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.LoadRaw});

                obj.Handles.UIElements.Basename = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.627 0.8025 0.325 0.03],...
                    'String',obj.Parameters.PreProcessing.Basename,...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.1 0.1 0.1],...
                    'HorizontalAlignment','left');

                obj.Handles.UIElements.Date_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.627 0.765 0.2 0.04],...
                    'String','Date / Time',...
                    'FontSize',14,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.1 0.1 0.1],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.Date_Value = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.73 0.765 0.125 0.04],...
                    'String',datestr(obj.Parameters.PreProcessing.Date),...
                    'FontSize',14,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.1 0.1 0.1],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.Duration_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.627 0.74 0.2 0.04],...
                    'String','Duration',...
                    'FontSize',14,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.3 0.3 0.3],'ForegroundColor',[0.1 0.1 0.1],...
                    'HorizontalAlignment','left');
                if ~isempty(obj.Parameters.PreProcessing.Duration),
                    Str = [num2str(obj.Parameters.PreProcessing.Duration) ' s'];
                else
                    Str = '';
                end
                obj.Handles.UIElements.Duration_Value = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.73 0.74 0.125 0.04],...
                    'String',Str,...
                    'FontSize',14,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.3 0.3 0.3],'ForegroundColor',[0.1 0.1 0.1],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.FrameRate_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.627 0.715 0.2 0.04],...
                    'String','Frame rate',...
                    'FontSize',14,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.1 0.1 0.1],...
                    'HorizontalAlignment','left');
                if ~isempty(obj.Parameters.PreProcessing.Period),
                    Str = [num2str(1/1000*round(1000*1/obj.Parameters.PreProcessing.Period)) ' Hz'];
                else
                    Str = '';
                end
                obj.Handles.UIElements.FrameRate_Value = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.73 0.715 0.125 0.04],...
                    'String',Str,...
                    'FontSize',14,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.1 0.1 0.1],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.FrameNumber_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.627 0.69 0.2 0.04],...
                    'String','Number of frames',...
                    'FontSize',14,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.3 0.3 0.3],'ForegroundColor',[0.1 0.1 0.1],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.FrameNumber_Value = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.73 0.69 0.125 0.04],...
                    'String',num2str(obj.Parameters.PreProcessing.FramesNum),...
                    'FontSize',14,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.3 0.3 0.3],'ForegroundColor',[0.1 0.1 0.1],...
                    'HorizontalAlignment','left');

                obj.Handles.UIElements.LEDPower_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.627 0.665 0.2 0.04],...
                    'String','LED Power',...
                    'FontSize',14,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.1 0.1 0.1],...
                    'HorizontalAlignment','left');
                if ~isempty(obj.Parameters.PreProcessing.Power),
                    Str = [num2str(obj.Parameters.PreProcessing.Power) ' mW'];
                else
                    Str = '';
                end
                obj.Handles.UIElements.LEDPower_Value = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.73 0.665 0.125 0.04],...
                    'String',Str,...
                    'FontSize',14,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.1 0.1 0.1],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.Gain_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.627 0.64 0.2 0.04],...
                    'String','Gain',...
                    'FontSize',14,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.3 0.3 0.3],'ForegroundColor',[0.1 0.1 0.1],...
                    'HorizontalAlignment','left');
                if ~isempty(obj.Parameters.PreProcessing.Duration),
                    Str = ['x' num2str(obj.Parameters.PreProcessing.Gain)];
                else
                    Str = '';
                end
                obj.Handles.UIElements.Gain_Value = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.73 0.64 0.125 0.04],...
                    'String',Str,...
                    'FontSize',14,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.3 0.3 0.3],'ForegroundColor',[0.1 0.1 0.1],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.Focus_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.627 0.615 0.2 0.04],...
                    'String','Focus',...
                    'FontSize',14,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.1 0.1 0.1],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.Focus_Value = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.73 0.615 0.125 0.04],...
                    'String',num2str(obj.Parameters.PreProcessing.Focus),...
                    'FontSize',14,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.1 0.1 0.1],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.LazyFix = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.627 0.59 0.325 0.04],...
                    'BackgroundColor',[0.2 0.2 0.2]);

                % Preprocessing parameters
                obj.Handles.UIElements.Parameters_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.6 0.55 0.25 0.03],...
                    'String','Preprocessing',...
                    'FontSize',16,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.7 0.7 0.7],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.TemporalDownS_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.627 0.5 0.2 0.04],...
                    'String','Temporal downsampling',...
                    'FontSize',15,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.TemporalDownS_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.76 0.505 0.075 0.04],...
                    'String',num2str(obj.Parameters.PreProcessing.Temporal_Downsampling),...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.TemporalDownS_CB});
                obj.Handles.UIElements.SpatialDownS_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.627 0.455 0.2 0.04],...
                    'String','Spatial downsampling',...
                    'FontSize',15,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.SpatialDownS_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.76 0.46 0.075 0.04],...
                    'String',num2str(obj.Parameters.PreProcessing.Spatial_Downsampling),...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.SpatialDownS_CB});
                obj.Handles.UIElements.FixPixels_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.627 0.41 0.2 0.04],...
                    'String','Fix defective pixels',...
                    'FontSize',15,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.FixPixels_CheckBox = uicontrol('Style','checkbox',...
                    'Units','Normalized','Position',[0.76+0.033 0.415 0.075/2 0.04],...
                    'String','',...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'Value',obj.Parameters.PreProcessing.Fix_Pixels,...
                    'Callback',{@(~,~)obj.FixPixels_CB});
                obj.Handles.UIElements.Crop_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.627 0.36 0.2 0.04],...
                    'String','Cropping',...
                    'FontSize',15,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.Crop_Draw = uicontrol('Style','pushbutton',...
                    'Units','Normalized','Position',[0.76 0.365 0.075 0.04],...
                    'String','Draw area',...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.Crop_Draw_CB});
                obj.Handles.UIElements.Crop_Delete = uicontrol('Style','pushbutton',...
                    'Units','Normalized','Position',[0.76+0.075 0.365 0.075 0.04],...
                    'String','Delete',...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Enable','off',...
                    'Callback',{@(~,~)obj.Crop_Delete_CB});

                % Validate
                obj.Handles.UIElements.ValidatePreProc_Button = uicontrol('Style','pushbutton',...
                    'Units','Normalized','Position',[0.627 0.25 0.075 0.04],...
                    'String','PreProcess',...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.15 0.15 0.15],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.PreProcess},...
                    'Enable','off');
            else
                % Update
                obj.Handles.UIElements.Basename.String = obj.Parameters.PreProcessing.Basename;
                obj.Handles.UIElements.Date_Value.String = datestr(obj.Parameters.PreProcessing.Date);
                obj.Handles.UIElements.Duration_Value.String = [num2str(obj.Parameters.PreProcessing.Duration) ' s'];
                obj.Handles.UIElements.FrameRate_Value.String = [num2str(1/1000*round(1000/obj.Parameters.PreProcessing.Period)) ' Hz'];
                obj.Handles.UIElements.FrameNumber_Value.String = num2str(obj.Parameters.PreProcessing.FramesNum);
                obj.Handles.UIElements.LEDPower_Value.String = [num2str(obj.Parameters.PreProcessing.Power) ' mW'];
                obj.Handles.UIElements.Gain_Value.String = ['x' num2str(obj.Parameters.PreProcessing.Gain)];
                obj.Handles.UIElements.Focus_Value.String = num2str(obj.Parameters.PreProcessing.Focus);
                obj.Handles.UIElements.FixPixels_CheckBox.Value = obj.Parameters.PreProcessing.Fix_Pixels;
                obj.Handles.UIElements.TemporalDownS_Edit.String = num2str(obj.Parameters.PreProcessing.Temporal_Downsampling);
                obj.Handles.UIElements.SpatialDownS_Edit.String = num2str(obj.Parameters.PreProcessing.Spatial_Downsampling);
            end

            % If we are changing tabs but have already (pre)processed:
            % update display
            if ~isempty(obj.Parameters.PreProcessing.RawFile),
                % Set movies structure
                obj.CleanMovies;
                obj.CurrentPlayers = struct(...
                    'Axes',obj.Handles.UIElements.MovieDisplay(1).Plot,...
                    'Data','Raw',...
                    'Filter',false,...
                    'Readers',[],...
                    'Plot',[],...
                    'AbsolutePosition',{obj.Handles.UIElements.MovieDisplay(1).AbsolutePosition});

                % Load movie
                obj.LoadMovies;

                % If cropping area is defined, add it
                if ~isempty(obj.Parameters.PreProcessing.Cropping),
                    obj.Handles.UIElements.CropArea = drawrectangle(obj.Handles.UIElements.MovieDisplay(1).Box,...
                        'FaceAlpha',0,...
                        'LineWidth',3,...
                        'Color',[0.6 0.3 0],...
                        'Position',obj.Parameters.PreProcessing.Cropping,...
                        'Deletable',false);
                    obj.Handles.UIElements.Crop_Delete.Enable = 'on';
                    obj.Handles.UIElements.Crop_Draw.Enable = 'off';
                end
                obj.Handles.UIElements.ValidatePreProc_Button.Enable = 'on';
            end
            drawnow
        end

        % Raw file selection
        function LoadRaw(obj)
            % Load file
            if isempty(obj.StartPath),
                Experimenters = DataBase.Lists.GetList('Experimenters');
                CurrentExperimenter = Experimenters(strcmpi(Experimenters(:,2),getenv('username')),1);
                [File, Path] = uigetfile({
                    '*.h5','Raw open scope file';'*.isxd','Raw file'},...
                    'Please choose a file...',['F:\' CurrentExperimenter{1} '\Data\CalciumImaging\']);
            else
                [File, Path] = uigetfile({
                    '*.h5','Raw open scope file';'*.isxd','Raw file'},...
                    'Please choose a file...',[obj.StartPath]);
            end
            if File == 0,
                return
            end

            % Check whether the filename seems correct
            PathLast = strsplit(Path,filesep);
            SessionName = [PathLast{end-2} '_' PathLast{end-1}];
            if ~contains(File,SessionName)
                % Judge and ask what to do
                FigCrow= figure('Units','normalized','Position',[0.025 0.35 0.25 0.35],'ToolBar','none','MenuBar','none');
                Crow = imread('F:\MATLAB\Common\Toolbox\General\CrowOfJudgement.jpg');
                image(Crow)
                GC = gca; GC.Position = [0 0 1 1]; axis(GC,'equal');
                GC.XAxis.Visible = 'off';GC.YAxis.Visible = 'off';
                Answer = questdlg(['The filename doesn''t seem to match the folder structure.' newline,...
                    '(if the naming is wrong, everything will have to be reprocessed)' newline,newline...
                    'What do you want to do?'],'Please choose...',...
                    'Return and double-check file and path names.',...
                    'Continue anyway.',...
                    'Return and double-check file and path names.');
                if strcmpi(Answer, 'Return and double-check file and path names.')
                    close(FigCrow)
                    return
                end
                close(FigCrow)
            end

            % Loosely check whether we have processed it already
            MetaFiles = dir([Path '**\*_StepMeta.mat']);

            % If something was found, warn about it and ask what to do
            if ~isempty(MetaFiles)
                Answer = questdlg(['Some processing was already performed with this file.' newline newline,...
                    'What do you want to do?'],'Please choose...',...
                    'Start processing from the beginning (the other files will not be deleted).',...
                    'Go back and directly load the already processed branch.',...
                    'Start processing from the beginning (the other files will not be deleted).');

                if strcmpi(Answer, 'Go back and directly load the already processed branch.'),
                    return
                end
            end

            % Delete cropping area if present
            if ~isempty(obj.Parameters.PreProcessing.Cropping)
                delete(obj.Handles.UIElements.CropArea);
                obj.Handles.UIElements.Crop_Delete.Enable = 'off';
                obj.Handles.UIElements.Crop_Draw.Enable = 'on';
            end

            % Start with a fresh parameters structure
            obj.Parameters = [];
            obj.Parameters.PreProcessing = obj.DefaultParameters.PreProcessing;
            obj.Parameters.PreProcessing.RawFile = fullfile(Path,File);
            obj.Parameters.PreProcessing.MainFolder = Path;
            if contains(File,'isxd')
                Basename = strsplit(File,'.isxd');
            else
                Basename = strsplit(File,'.h5');
            end
            obj.Parameters.PreProcessing.Basename = Basename{1};
            obj.Parameters.PreProcessing.ForkNumber = obj.GetFork('PreProcessing')+1;


            if contains(File,'isxd')
                RawMovie = isx.Movie.read(obj.Parameters.PreProcessing.RawFile);
                obj.Parameters.PreProcessing.Date =  RawMovie.timing.start.datetime;
                Times = RawMovie.timing.get_offsets_since_start;
                obj.Parameters.PreProcessing.TimeStamps = [Times.secs_float];
                obj.Parameters.PreProcessing.Duration = obj.Parameters.PreProcessing.TimeStamps(end);
                obj.Parameters.PreProcessing.FramesNum = RawMovie.timing.num_samples;
                obj.Parameters.PreProcessing.Period = RawMovie.timing.period.secs_float;
                obj.Parameters.PreProcessing.d1 = RawMovie.spacing.num_pixels(1);
                obj.Parameters.PreProcessing.d2 = RawMovie.spacing.num_pixels(2);
                AcqInfos = RawMovie.get_acquisition_info;

                if isempty(AcqInfos)
                    warning('Old versions of the files do not contain some metadata informations.')
                    obj.Parameters.PreProcessing.Exposure = [];
                    obj.Parameters.PreProcessing.Focus = [];
                    obj.Parameters.PreProcessing.Gain = [];
                    obj.Parameters.PreProcessing.Miniscope = [];
                    obj.Parameters.PreProcessing.Power = [];
                elseif isfield(AcqInfos,'MicroscopeEXLEDPower_mw_mm_2_')
                    obj.Parameters.PreProcessing.Exposure = AcqInfos.ExposureTime_ms_;
                    obj.Parameters.PreProcessing.Focus = AcqInfos.MicroscopeFocus;
                    obj.Parameters.PreProcessing.Gain = AcqInfos.MicroscopeGain;
                    obj.Parameters.PreProcessing.Miniscope = AcqInfos.MicroscopeType;
                    obj.Parameters.PreProcessing.Power = AcqInfos.MicroscopeEXLEDPower_mw_mm_2_;
                else
                    % Newest version, changes in their metadata
                    obj.Parameters.PreProcessing.Exposure = AcqInfos.ExposureTime_ms_;
                    obj.Parameters.PreProcessing.Focus = AcqInfos.MicroscopeFocus;
                    obj.Parameters.PreProcessing.Gain = AcqInfos.MicroscopeGain;
                    obj.Parameters.PreProcessing.Miniscope = AcqInfos.MicroscopeType;
                    obj.Parameters.PreProcessing.Power = AcqInfos.MicroscopeEXLED1Power_mw_mm_2_;
                end
            else
                H5I = h5info(obj.Parameters.PreProcessing.RawFile);
                FieldNames = {H5I.Datasets(:).Name};
                if ~(any(strcmpi(FieldNames,'mov')) && any(strcmpi(FieldNames,'times')))
                    warning('Could not find appropriate fields in the h5 file. Aborting.')
                end

                % Retrieve timestamps
                obj.Parameters.PreProcessing.TimeStamps = h5read(obj.Parameters.PreProcessing.RawFile,'/times');
                obj.Parameters.PreProcessing.Duration = obj.Parameters.PreProcessing.TimeStamps(end);

                % Retrieve metadata
                MetaFile = [obj.Parameters.PreProcessing.MainFolder, obj.Parameters.PreProcessing.Basename, '_CAIMetaData.csv'];
                if exist(MetaFile,'file')~=2
                    % Judge and ask what to do
                    FigCrow= figure('Units','normalized','Position',[0.025 0.35 0.25 0.35],'ToolBar','none','MenuBar','none');
                    Crow = imread('F:\MATLAB\Common\Toolbox\General\CrowOfJudgement.jpg');
                    image(Crow)
                    GC = gca; GC.Position = [0 0 1 1]; axis(GC,'equal');
                    GC.XAxis.Visible = 'off';GC.YAxis.Visible = 'off';
                    Answer = questdlg(['Could not find matching metadata file (MetaData.json) in session''s main folder.' newline,...
                        'What do you want to do?'],'Please choose...',...
                        'Return and double-check where the file is.',...
                        'Continue anyway (ONLY FOR QUICKLY CHECKING THE DATA).',...
                        'Return and double-check where the file is.');
                    close(FigCrow)
                    if strcmpi(Answer, 'Return and double-check file and path names.') || isempty(Answer)
                        return
                    end

                    % Deduce available missing info, set the rest blank
                    IndexMovie = find(strcmpi(FieldNames,'mov'));

                    obj.Parameters.PreProcessing.d1 = H5I.Datasets(IndexMovie).Dataspace.MaxSize(1);
                    obj.Parameters.PreProcessing.d2 = H5I.Datasets(IndexMovie).Dataspace.MaxSize(2);
                    obj.Parameters.PreProcessing.FramesNum = H5I.Datasets(IndexMovie).Dataspace.MaxSize(3);
                    obj.Parameters.PreProcessing.Date = NaN;
                    obj.Parameters.PreProcessing.Miniscope = '';


                    obj.Parameters.PreProcessing.Period = mean(diff(obj.Parameters.PreProcessing.TimeStamps));

                    obj.Parameters.PreProcessing.Exposure = [];
                    obj.Parameters.PreProcessing.Focus = [];
                    obj.Parameters.PreProcessing.Gain = [];
                    obj.Parameters.PreProcessing.Power = [];
                else
                    fid = fopen(MetaFile);
                    raw = fread(fid,inf);
                    str = char(raw');
                    fclose(fid);
                    MetaData = jsondecode(str);

                    obj.Parameters.PreProcessing.Date = datetime(...
                        MetaData.recordingStartTime.year,...
                        MetaData.recordingStartTime.month,...
                        MetaData.recordingStartTime.day,...
                        MetaData.recordingStartTime.hour,...
                        MetaData.recordingStartTime.minute,...
                        MetaData.recordingStartTime.second,...
                        MetaData.recordingStartTime.msec);
                    obj.Parameters.PreProcessing.Miniscope = MetaData.deviceType;

                    if ~isnumeric(MetaData.frameRate)
                        if contains(MetaData.frameRate,'FPS')
                            MetaData.frameRate = strsplit(MetaData.frameRate,'FPS');
                            MetaData.frameRate = str2double(MetaData.frameRate{1});
                        else
                            MetaData.frameRate = 0;
                        end
                    end

                    obj.Parameters.PreProcessing.Period = 1/MetaData.frameRate;

                    obj.Parameters.PreProcessing.d1 = MetaData.ROI.height;
                    obj.Parameters.PreProcessing.d2 = MetaData.ROI.width;
                    obj.Parameters.PreProcessing.FramesNum = numel(obj.Parameters.PreProcessing.TimeStamps);

                    obj.Parameters.PreProcessing.Exposure = [];
                    obj.Parameters.PreProcessing.Focus = MetaData.ewl;
                    obj.Parameters.PreProcessing.Gain = MetaData.gain;
                    obj.Parameters.PreProcessing.Power = MetaData.led0;
                end

            end

            % Update tabs parameters
            obj.CurrentTab = 'PreProcessing';
            obj.PreProcessing_Layout;

            % Set movies structure
            obj.CleanMovies;
            obj.CurrentPlayers = struct(...
                'Axes',obj.Handles.UIElements.MovieDisplay(1).Plot,...
                'Data','Raw',...
                'Filter',false,...
                'Readers',[],...
                'Plot',[],...
                'AbsolutePosition',obj.Handles.UIElements.MovieDisplay(1).AbsolutePosition);

            % Load movie
            obj.LoadMovies;

            % If cropping area is defined, add it
            if ~isempty(obj.Parameters.PreProcessing.Cropping),
                obj.Handles.UIElements.CropArea = drawrectangle(obj.Handles.UIElements.MovieDisplay(1).Box,...
                    'FaceAlpha',0,...
                    'LineWidth',3,...
                    'Color',[0.6 0.3 0],...
                    'Deletable',false,...
                    'Position',obj.Parameters.PreProcessing.Cropping);
                obj.Handles.UIElements.Crop_Delete.Enable = 'on';
                obj.Handles.UIElements.Crop_Draw.Enable = 'off';
            end

            % Update the tree display
            obj.UpdateTree;

            % Update the tab display
            Tabs = fieldnames(obj.Handles.Tabs);
            for F = 1 : numel(Tabs),
                obj.Handles.Tabs.(Tabs{F}).UserData = '0';
            end
            src.UserData = '1';
            obj.CurrentTab = 'PreProcessing';
            obj.TabAppearence;

            obj.Handles.UIElements.ValidatePreProc_Button.Enable = 'on';
            drawnow

        end

        function FixPixels_CB(obj)
            if obj.Handles.UIElements.FixPixels_CheckBox.Value == 1,
                obj.Parameters.PreProcessing.Fix_Pixels = true;
            else
                obj.Parameters.PreProcessing.Fix_Pixels = false;
            end
        end

        function TemporalDownS_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.TemporalDownS_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>=1,
                obj.Parameters.PreProcessing.Temporal_Downsampling = TempInput;
            else
                obj.Handles.UIElements.TemporalDownS_Edit.String = num2str(obj.Parameters.PreProcessing.Temporal_Downsampling);
            end
        end

        function SpatialDownS_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.SpatialDownS_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>=1,
                obj.Parameters.PreProcessing.Spatial_Downsampling = TempInput;
            else
                obj.Handles.UIElements.SpatialDownS_Edit.String = num2str(obj.Parameters.PreProcessing.Spatial_Downsampling);
            end
        end

        function Crop_Draw_CB(obj)
            PlotCrop = false;
            if ~isempty(obj.CurrentPlayers),
                if isfield(obj.Handles.UIElements,'CropArea')
                    if isempty(obj.Handles.UIElements.CropArea)
                        PlotCrop = true;
                    end
                else
                    PlotCrop = true;
                end
            end
            if PlotCrop
                CurrentPlayingStatus = obj.Playing;
                obj.Playing = false;
                obj.Handles.UIElements.CropArea = drawrectangle(obj.Handles.UIElements.MovieDisplay(1).Box,...
                    'FaceAlpha',0,...
                    'LineWidth',3,...
                    'Deletable',false,...
                    'Color',[0.6 0.3 0]);
                obj.Parameters.PreProcessing.Cropping = obj.Handles.UIElements.CropArea.Position;
                addlistener(obj.Handles.UIElements.CropArea,'MovingROI',@(~,~)obj.Crop_Edit_CB);
                drawnow
                obj.Handles.UIElements.Crop_Delete.Enable = 'on';
                obj.Handles.UIElements.Crop_Draw.Enable = 'off';
                obj.Playing = CurrentPlayingStatus;
                obj.PlayMovies;
            end
        end

        function Crop_Edit_CB(obj)
            obj.Parameters.PreProcessing.Cropping = obj.Handles.UIElements.CropArea.Position;
        end

        function Crop_Delete_CB(obj)
            if ~isempty(obj.Parameters.PreProcessing.Cropping),
                delete(obj.Handles.UIElements.CropArea)
                obj.Handles.UIElements.CropArea = [];
                obj.Parameters.PreProcessing.Cropping = [];
                obj.Handles.UIElements.Crop_Delete.Enable = 'off';
                obj.Handles.UIElements.Crop_Draw.Enable = 'on';
                drawnow
            end
        end

        %% Prefiltering tab & functions
        function Status = PreFiltering_Layout(obj)
            Status = false;
            if ~strcmpi(obj.CurrentTab,'PreFiltering'),
                % Check that we can actually switch to that tab...
                if ~isfield(obj.Parameters,'PreProcessing')|| isempty(obj.Parameters.PreProcessing.PreProcessingFile_isxd)
                    % Not preprocessed: we cannot change to the current tab
                    % (normally, the callback should be inactive, just a
                    % double-safety... and before that part is actually implemented)
                    return
                else
                    % At least, the session was preprocessed: we can
                    % change, but we might have to initialize
                    Status = true;
                    if ~isfield(obj.Parameters,'PreFiltering')
                        % Initialize prefiltering
                        obj.Parameters.PreFiltering = obj.DefaultParameters.PreFiltering;
                        if contains(obj.Parameters.PreProcessing.RawFile,'.isxd')
                            % Retrieve movie infos
                            PreProcessed_Movie = isx.Movie.read(obj.Parameters.PreProcessing.PreProcessingFile_isxd);
                            obj.Parameters.PreFiltering.d1 = PreProcessed_Movie.spacing.num_pixels(1);
                            obj.Parameters.PreFiltering.d2 = PreProcessed_Movie.spacing.num_pixels(2);
                            obj.Parameters.PreFiltering.FramesNum = PreProcessed_Movie.timing.num_samples;
                            Times = PreProcessed_Movie.timing.get_offsets_since_start;
                            obj.Parameters.PreFiltering.TimeStamps = ([Times.secs_float])';
                            obj.Parameters.PreFiltering.Duration = obj.Parameters.PreFiltering.TimeStamps(end);
                            obj.Parameters.PreFiltering.SamplingRate = 1/PreProcessed_Movie.timing.period.secs_float;
                        else
                            % Retrieve movie infos
                            H5I = h5info(obj.Parameters.PreProcessing.PreProcessingFile);
                            FieldNames = {H5I.Datasets(:).Name};
                            if ~(any(strcmpi(FieldNames,'mov')) && any(strcmpi(FieldNames,'times')))
                                warning('Could not find appropriate fields in the h5 file. Aborting.')
                            end

                            IndexMovie = find(strcmpi(FieldNames,'mov'));

                            obj.Parameters.PreFiltering.d1 = H5I.Datasets(IndexMovie).Dataspace.MaxSize(1);
                            obj.Parameters.PreFiltering.d2 = H5I.Datasets(IndexMovie).Dataspace.MaxSize(2);
                            obj.Parameters.PreFiltering.FramesNum = H5I.Datasets(IndexMovie).Dataspace.MaxSize(3);

                            obj.Parameters.PreFiltering.TimeStamps = h5read(obj.Parameters.PreProcessing.RawFile,'/times');
                            obj.Parameters.PreFiltering.Duration = obj.Parameters.PreFiltering.TimeStamps(end);
                            obj.Parameters.PreFiltering.SamplingRate = 1/(nanmean(diff(obj.Parameters.PreFiltering.TimeStamps)));
                        end
                        obj.UpdatePFKernel;
                        obj.Parameters.PreFiltering.ForkNumber = obj.GetFork('PreFiltering')+1;
                        % Adjust some default parameters to take into account potential
                        % downsampling
                        obj.Parameters.PreFiltering.Sig1 = obj.Parameters.PreFiltering.Sig1/(obj.Parameters.PreProcessing.Spatial_Downsampling/2);
                        obj.Parameters.PreFiltering.Sig2 = obj.Parameters.PreFiltering.Sig2/(obj.Parameters.PreProcessing.Spatial_Downsampling/2);
                        obj.Parameters.PreFiltering.WindowSize = obj.Parameters.PreFiltering.WindowSize/(obj.Parameters.PreProcessing.Spatial_Downsampling/2);
                    end
                end

                % Wipe current layout
                obj.WipeLayout;

                % Main display and uielements
                obj.SetMovieDisplay1;

                % Filtered movie, display only
                obj.SetMovieDisplay2;

                % Set movies structure
                obj.CleanMovies;

                obj.CurrentPlayers = struct(...
                    'Axes',{obj.Handles.UIElements.MovieDisplay(1).Plot,obj.Handles.UIElements.MovieDisplay(2).Plot },...
                    'Data','PreProcessed',...
                    'Filter',{false,true},...
                    'Kernel','PreFiltering',...
                    'Readers',[],...
                    'Plot',[],...
                    'AbsolutePosition',{obj.Handles.UIElements.MovieDisplay(1).AbsolutePosition,obj.Handles.UIElements.MovieDisplay(2).AbsolutePosition});

                % Load movie
                obj.LoadMovies;

                % Other uielements
                % Prefiltering parameters
                obj.Handles.UIElements.Parameters_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.605 0.375 0.25 0.03],...
                    'String','Filtering (Motion correction)',...
                    'FontSize',16,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.7 0.7 0.7],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.Sig1_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.627 0.3 0.2 0.04],...
                    'String','Sig1',...
                    'FontSize',15,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.Sig1_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.76 0.305 0.075 0.04],...
                    'String',num2str(obj.Parameters.PreFiltering.Sig1),...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.Sig1_CB});
                obj.Handles.UIElements.Sig2_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.627 0.255 0.2 0.04],...
                    'String','Sig2',...
                    'FontSize',15,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.Sig2_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.76 0.26 0.075 0.04],...
                    'String',num2str(num2str(obj.Parameters.PreFiltering.Sig2)),...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.Sig2_CB});
                obj.Handles.UIElements.Window_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.627 0.21 0.2 0.04],...
                    'String','Window size',...
                    'FontSize',15,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.Window_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.76 0.215 0.075 0.04],...
                    'String',num2str(num2str(obj.Parameters.PreFiltering.WindowSize)),...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.Window_CB});
                obj.Handles.UIElements.Crop_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.627 0.16 0.2 0.04],...
                    'String','Cropping',...
                    'FontSize',15,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.Crop_DrawPF = uicontrol('Style','pushbutton',...
                    'Units','Normalized','Position',[0.76 0.165 0.075 0.04],...
                    'String','Draw area',...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.Crop_DrawPF_CB});
                obj.Handles.UIElements.Crop_DeletePF = uicontrol('Style','pushbutton',...
                    'Units','Normalized','Position',[0.76+0.075 0.165 0.075 0.04],...
                    'String','Delete',...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Enable','off',...
                    'Callback',{@(~,~)obj.Crop_DeletePF_CB});
                obj.Handles.UIElements.Threshold_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.627 0.11 0.2 0.04],...
                    'String','Threshold',...
                    'FontSize',15,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.ThresholdPF_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.76 0.115 0.075 0.04],...
                    'String',num2str(num2str(obj.Parameters.PreFiltering.Threshold)),...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.ThresholdPF_CB});
                obj.Handles.UIElements.ThresholdPFHigh_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.76+0.075 0.115 0.075 0.04],...
                    'String','',...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.ThresholdPFHigh_CB});

                % Validate
                obj.Handles.UIElements.ValidatePreProc_Button = uicontrol('Style','pushbutton',...
                    'Units','Normalized','Position',[0.627 0.035 0.075 0.04],...
                    'String','Process',...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.15 0.15 0.15],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.PreFilter},...
                    'Enable','off');
                obj.ThresholdPF_CB;

            else
                % Update (e.g. if we load another session while in the tab)


            end
            % If cropping area is defined, add it
            if ~isempty(obj.Parameters.PreFiltering.Cropping),
                obj.Handles.UIElements.CropArea(2) = drawpolygon(obj.Handles.UIElements.MovieDisplay(2).Box,...
                    'Position',obj.Parameters.PreFiltering.Cropping,...
                    'FaceAlpha',0,...
                    'LineWidth',3,...
                    'Color',[0.6 0.3 0],...
                    'UserData','2');
                obj.Handles.UIElements.CropArea(1) = drawpolygon(obj.Handles.UIElements.MovieDisplay(1).Box,...
                    'Position',obj.Parameters.PreFiltering.Cropping,...
                    'FaceAlpha',0,...
                    'LineWidth',3,...
                    'Color',[0.6 0.3 0],...
                    'UserData','1');
                addlistener(obj.Handles.UIElements.CropArea(:),'MovingROI',@(src,evt)obj.Crop_EditPF_CB(src,evt));
                obj.Handles.UIElements.Crop_DeletePF.Enable = 'on';
                obj.Handles.UIElements.Crop_DrawPF.Enable = 'off';
            end
            obj.Handles.UIElements.ValidatePreProc_Button.Enable = 'on';
        end

        function UpdatePFKernel(obj)
            if obj.Parameters.PreFiltering.Sig1>0 && obj.Parameters.PreFiltering.Sig2>0
                gaussian1 = fspecial('Gaussian', obj.Parameters.PreFiltering.WindowSize, obj.Parameters.PreFiltering.Sig1);
                gaussian2 = fspecial('Gaussian', obj.Parameters.PreFiltering.WindowSize, obj.Parameters.PreFiltering.Sig2);
                if obj.Parameters.PreFiltering.Sig1>obj.Parameters.PreFiltering.Sig2,
                    psf = gaussian1 - gaussian2;
                else
                    psf = gaussian2 - gaussian1;
                end
            elseif obj.Parameters.PreFiltering.Sig1 == 0
                psf = fspecial('Gaussian', obj.Parameters.PreFiltering.WindowSize, obj.Parameters.PreFiltering.Sig2);
            else
                psf = fspecial('Gaussian', obj.Parameters.PreFiltering.WindowSize, obj.Parameters.PreFiltering.Sig1);
            end
            ind_nonzero = (psf(:)>=max(psf(:,1)));
            psf = psf-mean(psf(ind_nonzero));
            psf(~ind_nonzero) = 0;
            obj.Parameters.PreFiltering.Kernel = psf;
        end

        function Sig1_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.Sig1_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>0,
                obj.Parameters.PreFiltering.Sig1 = TempInput;
                obj.UpdatePFKernel;
            elseif TempInput == 0
                if obj.Parameters.PreFiltering.Sig2 == 0,
                    obj.Handles.UIElements.Sig1_Edit.String = num2str(obj.Parameters.PreFiltering.Sig1);
                else
                    obj.Parameters.PreFiltering.Sig1 = TempInput;
                    obj.UpdatePFKernel;
                end
            else
                obj.Handles.UIElements.Sig1_Edit.String = num2str(obj.Parameters.PreFiltering.Sig1);
            end
            obj.UpdateDisplay;
            obj.Handles.UIElements.MovieDisplay(2).Plot.CLim(2) = max(obj.Handles.UIElements.MovieDisplay(2).Plot.Children.CData,[],'all');
        end

        function Sig2_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.Sig2_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>=0,
                obj.Parameters.PreFiltering.Sig2 = TempInput;
                obj.UpdatePFKernel;
            elseif TempInput == 0
                if obj.Parameters.PreFiltering.Sig1 == 0,
                    obj.Handles.UIElements.Sig2_Edit.String = num2str(obj.Parameters.PreFiltering.Sig2);
                else
                    obj.Parameters.PreFiltering.Sig2 = TempInput;
                    obj.UpdatePFKernel;
                end
            else
                obj.Handles.UIElements.Sig2_Edit.String = num2str(obj.Parameters.PreFiltering.Sig2);
            end
            obj.UpdateDisplay;
        end

        function Window_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.Window_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>=1,
                obj.Parameters.PreFiltering.WindowSize = TempInput;
                obj.UpdatePFKernel;
            else
                obj.Handles.UIElements.Window_Edit.String = num2str(obj.Parameters.PreFiltering.WindowSize);
            end
            obj.UpdateDisplay;
        end

        function ThresholdPF_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.ThresholdPF_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput<obj.Handles.UIElements.MovieDisplay(2).Plot.CLim(2),
                obj.Parameters.PreFiltering.Threshold(1) = TempInput;
                obj.Handles.UIElements.MovieDisplay(2).Plot.CLim(1) = TempInput;
            else
                obj.Handles.UIElements.ThresholdPF_Edit.String = num2str(obj.Parameters.PreFiltering.Threshold(1));
            end
            obj.UpdateDisplay;
        end

        function ThresholdPFHigh_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.ThresholdPFHigh_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>obj.Handles.UIElements.MovieDisplay(2).Plot.CLim(1),
                obj.Handles.UIElements.MovieDisplay(2).Plot.CLim(2) = TempInput;
                obj.Parameters.PreFiltering.Threshold(2) = TempInput;
            else
                obj.Handles.UIElements.ThresholdPFHigh_Edit.String = num2str( obj.Parameters.PreFiltering.Threshold(2));
            end
            obj.UpdateDisplay;
        end

        function Crop_DrawPF_CB(obj)
            if ~isempty(obj.CurrentPlayers),
                obj.CurrentPlayingStatus = obj.Playing;
                obj.Playing = false;
                % Disable all callbacks
                obj.DisableAll;

                % Set a callback for the two axes - we want to use the
                % proper axes for the drawpolygon
                obj.Handles.UIElements.MovieDisplay(1).Box.ButtonDownFcn = {@(~,~)obj.MovieDisplay1ButtonDown_PF_CB};
                obj.Handles.UIElements.MovieDisplay(2).Box.ButtonDownFcn = {@(~,~)obj.MovieDisplay2ButtonDown_PF_CB};
            end
        end

        function MovieDisplay1ButtonDown_PF_CB(obj)
            obj.Handles.UIElements.CropArea(1) = drawpolygon(obj.Handles.UIElements.MovieDisplay(1).Box,...
                'FaceAlpha',0,...
                'LineWidth',3,...
                'Color',[0.6 0.3 0],...
                'UserData','1');
            obj.Parameters.PreFiltering.Cropping = obj.Handles.UIElements.CropArea(1).Position;
            obj.Handles.UIElements.CropArea(2) = drawpolygon(obj.Handles.UIElements.MovieDisplay(2).Box,...
                'Position',obj.Parameters.PreFiltering.Cropping,...
                'FaceAlpha',0,...
                'LineWidth',3,...
                'Color',[0.6 0.3 0],...
                'UserData','2');
            addlistener(obj.Handles.UIElements.CropArea(:),'MovingROI',@(src,evt)obj.Crop_EditPF_CB(src,evt));
            drawnow
            obj.EnableAll;
            obj.Handles.UIElements.Crop_DeletePF.Enable = 'on';
            obj.Handles.UIElements.Crop_DrawPF.Enable = 'off';
            obj.Playing = obj.CurrentPlayingStatus;
            obj.PlayMovies;
        end

        function MovieDisplay2ButtonDown_PF_CB(obj)
            obj.Handles.UIElements.CropArea(2) = drawpolygon(obj.Handles.UIElements.MovieDisplay(2).Box,...
                'FaceAlpha',0,...
                'LineWidth',3,...
                'Color',[0.6 0.3 0],...
                'UserData','2');
            obj.Parameters.PreFiltering.Cropping = obj.Handles.UIElements.CropArea(2).Position;
            obj.Handles.UIElements.CropArea(1) = drawpolygon(obj.Handles.UIElements.MovieDisplay(1).Box,...
                'Position',obj.Parameters.PreFiltering.Cropping,...
                'FaceAlpha',0,...
                'LineWidth',3,...
                'Color',[0.6 0.3 0],...
                'UserData','1');
            addlistener(obj.Handles.UIElements.CropArea(:),'MovingROI',@(src,evt)obj.Crop_EditPF_CB(src,evt));
            drawnow
            obj.EnableAll;
            obj.Handles.UIElements.Crop_DeletePF.Enable = 'on';
            obj.Handles.UIElements.Crop_DrawPF.Enable = 'off';
            obj.Playing = obj.CurrentPlayingStatus;
            obj.PlayMovies;
        end

        function Crop_EditPF_CB(obj,src,evt)
            obj.Parameters.PreFiltering.Cropping = obj.Handles.UIElements.CropArea(str2double(src.UserData)).Position;
            if str2double(src.UserData) == 1,
                obj.Handles.UIElements.CropArea(2).Position = obj.Parameters.PreFiltering.Cropping;
            else
                obj.Handles.UIElements.CropArea(1).Position = obj.Parameters.PreFiltering.Cropping;
            end
        end

        function Crop_DeletePF_CB(obj)
            if ~isempty(obj.Parameters.PreFiltering.Cropping),
                delete(obj.Handles.UIElements.CropArea(:))
                obj.Parameters.PreFiltering.Cropping = [];
                obj.Handles.UIElements.Crop_DeletePF.Enable = 'off';
                obj.Handles.UIElements.Crop_DrawPF.Enable = 'on';
                drawnow
            end
        end

        %% Motion correction tab & functions
        function Status = MotionCorrection_Layout(obj)
            Status = false;
            if ~strcmpi(obj.CurrentTab,'MotionCorrection'),
                % Check that we can actually switch to that tab...
                if ~isfield(obj.Parameters,'PreFiltering') || isempty(obj.Parameters.PreFiltering.PreFilteringFile)
                    % Not prefiltered: we cannot change to the current tab
                    % (normally, the callback should be inactive, just a
                    % double-safety... and before that part is actually implemented)
                    return
                else
                    Status = true;
                    % At least, the session was prefiltered: we can
                    % change, but we might have to initialize
                    if ~isfield(obj.Parameters,'MotionCorrection'),
                        % Initialize motion correction step
                        obj.Parameters.MotionCorrection = obj.DefaultParameters.MotionCorrection;

                        % Adjust some default parameters to take into account potential
                        % downsampling
                        obj.Parameters.MotionCorrection.grid_size([1 2]) = obj.Parameters.MotionCorrection.grid_size([1 2])/(obj.Parameters.PreProcessing.Spatial_Downsampling);
                        obj.Parameters.MotionCorrection.max_dev([1 2]) = obj.Parameters.MotionCorrection.max_dev([1 2])/(obj.Parameters.PreProcessing.Spatial_Downsampling);
                        obj.Parameters.MotionCorrection.max_shift([1 2]) = obj.Parameters.MotionCorrection.max_shift([1 2])/(obj.Parameters.PreProcessing.Spatial_Downsampling);
                        obj.Parameters.MotionCorrection.overlap_post([1 2]) = obj.Parameters.MotionCorrection.overlap_post([1 2])/(obj.Parameters.PreProcessing.Spatial_Downsampling);
                        obj.Parameters.MotionCorrection.overlap_pre([1 2]) = obj.Parameters.MotionCorrection.overlap_pre([1 2])/(obj.Parameters.PreProcessing.Spatial_Downsampling);

                        % Create Visualization field
                        obj.Parameters.Visualization.Sig1 = obj.Parameters.PreFiltering.Sig1;
                        obj.Parameters.Visualization.Sig2 = obj.Parameters.PreFiltering.Sig2;
                        obj.Parameters.Visualization.WindowSize = obj.Parameters.PreFiltering.WindowSize;
                        obj.Parameters.Visualization.ThresholdLow = obj.Parameters.PreFiltering.Threshold(1);
                        obj.Parameters.Visualization.ThresholdHigh = [];
                        obj.Parameters.Visualization.Kernel = obj.Parameters.PreFiltering.Kernel;
                    end
                end
                if ~isfield(obj.Parameters,'Visualization'),
                    % Create Visualization field
                    obj.Parameters.Visualization.Sig1 = obj.Parameters.PreFiltering.Sig1;
                    obj.Parameters.Visualization.Sig2 = obj.Parameters.PreFiltering.Sig2;
                    obj.Parameters.Visualization.WindowSize = obj.Parameters.PreFiltering.WindowSize;
                    obj.Parameters.Visualization.ThresholdLow = obj.Parameters.PreFiltering.Threshold(1);
                    obj.Parameters.Visualization.ThresholdHigh = [];
                    obj.Parameters.Visualization.Kernel = obj.Parameters.PreFiltering.Kernel;
                end


                % Wipe current layout
                obj.WipeLayout;

                % Main display and uielements
                obj.SetMovieDisplay1;

                % Motion corrected movie, display only
                obj.SetMovieDisplay2;

                % Set movies structure
                obj.CleanMovies;

                if isempty(obj.Parameters.MotionCorrection.MotionCorrectionFile)
                    obj.CurrentPlayers = struct(...
                        'Axes',{obj.Handles.UIElements.MovieDisplay(1).Plot},...
                        'Data','PreProcessed',...
                        'Filter',{true},...
                        'Kernel','Visualization',...
                        'Readers',[],...
                        'Plot',[],...
                        'AbsolutePosition',{obj.Handles.UIElements.MovieDisplay(1).AbsolutePosition});
                else
                    obj.CurrentPlayers = struct(...
                        'Axes',{obj.Handles.UIElements.MovieDisplay(1).Plot,obj.Handles.UIElements.MovieDisplay(2).Plot},...
                        'Data',{'PreProcessed','MotionCorrected'},...
                        'Filter',true,...
                        'Kernel','Visualization',...
                        'Readers',[],...
                        'Plot',[],...
                        'AbsolutePosition',{obj.Handles.UIElements.MovieDisplay(1).AbsolutePosition,obj.Handles.UIElements.MovieDisplay(2).AbsolutePosition});
                end

                % Load movie
                obj.LoadMovies;

                % Other uielements
                % Display options
                obj.Handles.UIElements.Parameters_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.1950 0.255 0.25 0.03],...
                    'String','Display options (not for processing)',...
                    'FontSize',16,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.7 0.7 0.7],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.Sig1_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.217 0.210 0.2 0.035],...
                    'String','Sig1',...
                    'FontSize',15,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.Sig1_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.350 0.215 0.075 0.035],...
                    'String',num2str(obj.Parameters.Visualization.Sig1),...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.Sig1_Visu_CB});
                obj.Handles.UIElements.Sig2_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.217 0.17 0.2 0.035],...
                    'String','Sig2',...
                    'FontSize',15,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.Sig2_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.350 0.175 0.075 0.035],...
                    'String',num2str(num2str(obj.Parameters.Visualization.Sig2)),...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.Sig2_Visu_CB});
                obj.Handles.UIElements.Window_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.217 0.13 0.2 0.035],...
                    'String','Window size',...
                    'FontSize',15,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.Window_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.350 0.135 0.075 0.035],...
                    'String',num2str(num2str(obj.Parameters.Visualization.WindowSize)),...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.Window_Visu_CB});
                obj.Handles.UIElements.Threshold_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.217 0.09 0.2 0.035],...
                    'String','Color limit',...
                    'FontSize',15,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.ThresholdLow_Visu_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.350 0.095 0.075 0.035],...
                    'String',num2str(num2str(obj.Parameters.Visualization.ThresholdLow)),...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.ThresholdLow_Visu_CB});
                obj.Handles.UIElements.ThresholdHigh_Visu_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.350+0.075 0.095 0.075 0.035],...
                    'String',num2str(num2str(obj.Parameters.Visualization.ThresholdHigh)),...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.ThresholdHigh_Visu_CB});
                obj.Handles.UIElements.GridDisplay_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.217 0.05 0.2 0.035],...
                    'String','Display representative tile',...
                    'FontSize',15,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.GridDisplay_Check = uicontrol('Style','checkbox',...
                    'Units','Normalized','Position',[0.350+0.0325 0.055 0.075 0.035],...
                    'String','',...
                    'Value',1,...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.GridDisplay_CB});

                % Motion Correction parameters
                obj.Handles.UIElements.Parameters_LegendP = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.605 0.375 0.25 0.03],...
                    'String','Motion correction: Parameters',...
                    'FontSize',16,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.7 0.7 0.7],...
                    'HorizontalAlignment','left');
                % bin_width
                obj.Handles.UIElements.bin_width_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.627 0.33 0.2 0.035],...
                    'String','Bin width',...
                    'Tooltip','Width of each bin',...
                    'FontSize',15,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.bin_width_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.76 0.335 0.075 0.035],...
                    'String',num2str(obj.Parameters.MotionCorrection.bin_width),...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.bin_width_CB});
                % grid_size
                obj.Handles.UIElements.grid_size_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.627 0.29 0.2 0.035],...
                    'String','Grid size',...
                    'Tooltip','Size of non-overlapping regions',...
                    'FontSize',15,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.grid_size_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.76 0.295 0.075 0.035],...
                    'String',num2str(num2str(obj.Parameters.MotionCorrection.grid_size(1))),...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.grid_size_CB});
                % init_batch
                obj.Handles.UIElements.init_batch_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.627 0.25 0.2 0.035],...
                    'String','Initial batch',...
                    'Tooltip','Length of initial batch',...
                    'FontSize',15,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.init_batch_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.76 0.255 0.075 0.035],...
                    'String',num2str(num2str(obj.Parameters.MotionCorrection.init_batch)),...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.init_batch_CB});
                % iter
                obj.Handles.UIElements.iter_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.627 0.21 0.2 0.035],...
                    'String','Iteration',...
                    'Tooltip','Number of data passes',...
                    'FontSize',15,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.iter_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.76 0.215 0.075 0.035],...
                    'String',num2str(num2str(obj.Parameters.MotionCorrection.iter)),...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.iter_CB});
                % max_dev
                obj.Handles.UIElements.max_dev_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.627 0.17 0.2 0.035],...
                    'String','Maximum deviation',...
                    'Tooltip','Maximum deviation of patch shift from rigid shift',...
                    'FontSize',15,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.max_dev_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.76 0.175 0.075 0.035],...
                    'String',num2str(num2str(obj.Parameters.MotionCorrection.max_dev(1))),...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.max_dev_CB});
                % max_shift
                obj.Handles.UIElements.max_shift_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.627 0.13 0.2 0.035],...
                    'String','Maximum rigid',...
                    'Tooltip','Maximum rigid shift in each direction',...
                    'FontSize',15,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.max_shift_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.76 0.135 0.075 0.035],...
                    'String',num2str(num2str(obj.Parameters.MotionCorrection.max_shift(1))),...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.max_shift_CB});
                % overlap_pre
                obj.Handles.UIElements.overlap_pre_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.627 0.09 0.2 0.035],...
                    'String','Overlap (pre)',...
                    'Tooltip','Length of initial batch',...
                    'FontSize',15,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.overlap_pre_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.76 0.095 0.075 0.035],...
                    'String',num2str(num2str(obj.Parameters.MotionCorrection.overlap_pre(1))),...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.overlap_pre_CB});
                % overlap_post
                obj.Handles.UIElements.overlap_post_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.627 0.05 0.2 0.035],...
                    'String','Overlap (post)',...
                    'Tooltip','Length of initial batch',...
                    'FontSize',15,'FontName','Arial','FontWeight','n',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.overlap_post_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.76 0.055 0.075 0.035],...
                    'String',num2str(num2str(obj.Parameters.MotionCorrection.overlap_post(1))),...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.overlap_post_CB});

                % Validate
                obj.Handles.UIElements.ValidateMotionCorrection_Button = uicontrol('Style','pushbutton',...
                    'Units','Normalized','Position',[0.965-0.075 0.055 0.075 0.035],...
                    'String','Process',...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.15 0.15 0.15],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.MotionCorrect},...
                    'Enable','off');
                obj.ThresholdLow_Visu_CB;
                obj.ThresholdHigh_Visu_CB;
                % Plot representative tile (just centered on the FOV; it is
                % probably possible to retrieve the exact tiling used by
                % NormCorre but that's not a priority right now)
                Center = 0.5 * [obj.Parameters.PreFiltering.d2,obj.Parameters.PreFiltering.d1];
                obj.Handles.UIElements.TileOverlap = fill(...
                    [Center(1) - 0.5*obj.Parameters.MotionCorrection.grid_size(1) + obj.Parameters.MotionCorrection.overlap_pre(1),...
                    Center(1) + 0.5*obj.Parameters.MotionCorrection.grid_size(1) - obj.Parameters.MotionCorrection.overlap_pre(1),...
                    Center(1) + 0.5*obj.Parameters.MotionCorrection.grid_size(1) - obj.Parameters.MotionCorrection.overlap_pre(1),...
                    Center(1) - 0.5*obj.Parameters.MotionCorrection.grid_size(1) + obj.Parameters.MotionCorrection.overlap_pre(1)],...
                    [Center(2) - 0.5*obj.Parameters.MotionCorrection.grid_size(2) + obj.Parameters.MotionCorrection.overlap_pre(2),...
                    Center(2) - 0.5*obj.Parameters.MotionCorrection.grid_size(2) + obj.Parameters.MotionCorrection.overlap_pre(2),...
                    Center(2) + 0.5*obj.Parameters.MotionCorrection.grid_size(2) - obj.Parameters.MotionCorrection.overlap_pre(2),...
                    Center(2) + 0.5*obj.Parameters.MotionCorrection.grid_size(2) - obj.Parameters.MotionCorrection.overlap_pre(2)],...
                    [0.6 0.3 0],'EdgeColor',[0.6 0.3 0],'LineWidth',0.5,'FaceAlpha',0.1,...
                    'Parent',obj.Handles.UIElements.MovieDisplay(1).Box);
                obj.Handles.UIElements.GridTile = fill(...
                    [Center(1) - 0.5*obj.Parameters.MotionCorrection.grid_size(1),...
                    Center(1) + 0.5*obj.Parameters.MotionCorrection.grid_size(1),...
                    Center(1) + 0.5*obj.Parameters.MotionCorrection.grid_size(1),...
                    Center(1) - 0.5*obj.Parameters.MotionCorrection.grid_size(1)],...
                    [Center(2) - 0.5*obj.Parameters.MotionCorrection.grid_size(2),...
                    Center(2) - 0.5*obj.Parameters.MotionCorrection.grid_size(2),...
                    Center(2) + 0.5*obj.Parameters.MotionCorrection.grid_size(2),...
                    Center(2) + 0.5*obj.Parameters.MotionCorrection.grid_size(2)],...
                    [0.6 0.3 0],'EdgeColor',[0.6 0.3 0],'LineWidth',1.5,'FaceAlpha',0.1,...
                    'Parent',obj.Handles.UIElements.MovieDisplay(1).Box);

                obj.Handles.UIElements.ValidateMotionCorrection_Button.Enable = 'on';
            else
                % Update (e.g. if we load another session while in the tab)


            end

        end

        function GridUpdate(obj)
            Center = 0.5 * [obj.Parameters.PreFiltering.d2,obj.Parameters.PreFiltering.d1];
            obj.Handles.UIElements.TileOverlap.XData = [...
                Center(1) - 0.5*obj.Parameters.MotionCorrection.grid_size(1) + obj.Parameters.MotionCorrection.overlap_pre(1),...
                Center(1) + 0.5*obj.Parameters.MotionCorrection.grid_size(1) - obj.Parameters.MotionCorrection.overlap_pre(1),...
                Center(1) + 0.5*obj.Parameters.MotionCorrection.grid_size(1) - obj.Parameters.MotionCorrection.overlap_pre(1),...
                Center(1) - 0.5*obj.Parameters.MotionCorrection.grid_size(1) + obj.Parameters.MotionCorrection.overlap_pre(1)];
            obj.Handles.UIElements.TileOverlap.YData = [...
                Center(2) - 0.5*obj.Parameters.MotionCorrection.grid_size(2) + obj.Parameters.MotionCorrection.overlap_pre(2),...
                Center(2) - 0.5*obj.Parameters.MotionCorrection.grid_size(2) + obj.Parameters.MotionCorrection.overlap_pre(2),...
                Center(2) + 0.5*obj.Parameters.MotionCorrection.grid_size(2) - obj.Parameters.MotionCorrection.overlap_pre(2),...
                Center(2) + 0.5*obj.Parameters.MotionCorrection.grid_size(2) - obj.Parameters.MotionCorrection.overlap_pre(2)];
            obj.Handles.UIElements.GridTile.XData = [...
                Center(1) - 0.5*obj.Parameters.MotionCorrection.grid_size(1),...
                Center(1) + 0.5*obj.Parameters.MotionCorrection.grid_size(1),...
                Center(1) + 0.5*obj.Parameters.MotionCorrection.grid_size(1),...
                Center(1) - 0.5*obj.Parameters.MotionCorrection.grid_size(1)];
            obj.Handles.UIElements.GridTile.YData = [...
                Center(2) - 0.5*obj.Parameters.MotionCorrection.grid_size(2),...
                Center(2) - 0.5*obj.Parameters.MotionCorrection.grid_size(2),...
                Center(2) + 0.5*obj.Parameters.MotionCorrection.grid_size(2),...
                Center(2) + 0.5*obj.Parameters.MotionCorrection.grid_size(2)];
            drawnow
        end

        function GridDisplay_CB(obj)
            if  obj.Handles.UIElements.GridDisplay_Check.Value,
                obj.Handles.UIElements.TileOverlap.Visible = 'on';
                obj.Handles.UIElements.GridTile.Visible = 'on';
            else
                obj.Handles.UIElements.TileOverlap.Visible = 'off';
                obj.Handles.UIElements.GridTile.Visible = 'off';
            end
        end

        function bin_width_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.bin_width_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>=1,
                obj.Parameters.MotionCorrection.bin_width = TempInput;
            else
                obj.Handles.UIElements.bin_width_Edit.String = num2str(obj.Parameters.MotionCorrection.bin_width);
            end
        end

        function grid_size_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.grid_size_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>=1,
                obj.Parameters.MotionCorrection.grid_size = [TempInput TempInput 1];
                obj.GridUpdate;
            else
                obj.Handles.UIElements.grid_size_Edit.String = num2str(obj.Parameters.MotionCorrection.grid_size(1));
            end
        end

        function init_batch_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.init_batch_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>=1,
                obj.Parameters.MotionCorrection.init_batch = TempInput;
            else
                obj.Handles.UIElements.init_batch_Edit.String = num2str(obj.Parameters.MotionCorrection.init_batch);
            end
        end

        function iter_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.iter_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>=1,
                obj.Parameters.MotionCorrection.iter = TempInput;
            else
                obj.Handles.UIElements.iter_Edit.String = num2str(obj.Parameters.MotionCorrection.iter);
            end
        end

        function max_dev_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.max_dev_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>=1,
                obj.Parameters.MotionCorrection.max_dev = [TempInput TempInput 1];
            else
                obj.Handles.UIElements.max_dev_Edit.String = num2str(obj.Parameters.MotionCorrection.max_dev(1));
            end
        end

        function max_shift_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.max_shift_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>=1,
                obj.Parameters.MotionCorrection.max_shift = [TempInput TempInput 1];
            else
                obj.Handles.UIElements.max_shift_Edit.String = num2str(obj.Parameters.MotionCorrection.max_shift(1));
            end
        end

        function overlap_pre_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.overlap_pre_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>=1,
                obj.Parameters.MotionCorrection.overlap_pre = [TempInput TempInput 1];
                obj.GridUpdate;
            else
                obj.Handles.UIElements.overlap_pre_Edit.String = num2str(obj.Parameters.MotionCorrection.overlap_pre(1));
            end
        end

        function overlap_post_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.overlap_post_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>=1,
                obj.Parameters.MotionCorrection.overlap_post = [TempInput TempInput 1];
            else
                obj.Handles.UIElements.overlap_post_Edit.String = num2str(obj.Parameters.MotionCorrection.overlap_post(1));
            end
        end

        %% CNMFE tab & functions

        function Status = CNMFE_Layout(obj)
            Status = false;
            if ~strcmpi(obj.CurrentTab,'CNMFE'),
                % Check that we can actually switch to that tab...
                if ~isfield(obj.Parameters,'MotionCorrection') || isempty(obj.Parameters.MotionCorrection.MotionCorrectionFile)
                    % Not motion corrected: we cannot change to the current tab
                    % (normally, the callback should be inactive, just a
                    % double-safety... and before that part is actually implemented)
                    return
                else
                    Status = true;
                    % At least, the session was prefiltered: we can
                    % change, but we might have to initialize

                    if ~isfield(obj.Parameters,'CNMFE'),

                        % Initialize CNMF-E step
                        obj.Parameters.CNMFE = obj.DefaultParameters.CNMFE;

                        % Create a Source2D object (CNMFE main class)
                        Neuron = Sources2D;

                        % Set the tag for transfer and indicating the
                        % analysis is fully performed to false
                        obj.Parameters.CNMFE.Extracted = false;

                        % Set some of the core parameters' values from the
                        % previous steps
                        obj.Parameters.CNMFE.gSig = max([obj.Parameters.PreFiltering.Sig1 obj.Parameters.PreFiltering.Sig2]);
                        obj.Parameters.CNMFE.gSiz = obj.Parameters.CNMFE.gSig*4+1; % same as CNMFE is approximating it for center_psf true

                        % CNMFE parameters are *a bit* messy to set, and the
                        % function setting the default values is also a bit
                        % hard to handle as the values are in a list on
                        % their own -to play it safe, let's just replace
                        % the values as intended, after the instance is
                        % created; the following lines are actually adapted
                        % from the demo code


                        % -------------------------      SPATIAL      -------------------------  %
                        % NB: Formula CNMFE: fspecial('gaussian', round(gSiz), gSig) AND center_psf
                        % false; otherwise, uses ceil(gSig*4+1), gSig
                        %   Formula Motion Correction: psf = fspecial('gaussian', round(2*gSiz),
                        %   gSig); with e.g. for SubS 2x gSig = 2; gSiz = 10;
                        Neuron.options.gSig = obj.Parameters.CNMFE.gSig;               % pixel, gaussian width of a gaussian kernel for filtering the data. 0 means no filtering
                        Neuron.options.gSiz =  obj.Parameters.CNMFE.gSiz;              % pixel, neuron diameter
                        Neuron.options.ssub = 1;           % spatial downsampling factor
                        % determine the search locations by selecting a round area
                        Neuron.options.search_method = 'ellipse';
                        Neuron.options.dist = 5;
                        Neuron.options.bSiz = ceil(Neuron.options.dist/obj.Parameters.PreProcessing.Spatial_Downsampling);
                        Neuron.options.spatial_constraints = struct('connected', true, 'circular', false);  % you can include following constraints: 'circular'
                        Neuron.options.spatial_algorithm = 'hals_thresh';
                        % -------------------------      TEMPORAL     -------------------------  %
                        Neuron.options.tsub = 1;                                       % temporal downsampling factor
                        Neuron.options.deconv_flag = obj.Parameters.CNMFE.Constrained; % run deconvolution or not; set to false for unconstrained analysis
                        Neuron.options.deconv_options = struct('type', 'ar1', ...      % model of the calcium traces. {'ar1', 'ar2'}
                            'method', 'foopsi', ...                     % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
                            'smin', -5, ...                             % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level
                            'optimize_pars', true, ...                  % optimize AR coefficients
                            'optimize_b', true, ...                     % optimize the baseline);
                            'max_tau', obj.Parameters.CNMFE.max_tau);   % maximum decay time (unit: frame);
                        Neuron.Fs = obj.Parameters.PreFiltering.SamplingRate;
                        Neuron.options.nk = obj.Parameters.CNMFE.detrend_nk;             % detrending the slow fluctuation. usually 1 is fine (no detrending)
                        Neuron.options.detrend_method = 'spline';  % compute the local minimum as an estimation of trend.

                        % -------------------------     BACKGROUND    -------------------------  %
                        Neuron.options.background_model = 'ring';  % model of the background {'ring', 'svd'(default), 'nmf'}
                        Neuron.options.nb = 1;             % number of background sources for each patch (only be used in SVD and NMF model)
                        Neuron.options.ring_radius = ceil(0.5*1.05 * obj.Parameters.CNMFE.gSiz);  % when the ring model used, it is the radius of the ring used in the background model.
                        Neuron.options.num_neighbors = []; % number of neighbors for each neuron
                        Neuron.options.bg_ssub = 2;        % downsample background for a faster speed... BUT LEADS TO ERROR in update_background_parallel: either debug or leave to 1

                        % -------------------------      MERGING      -------------------------  %
                        Neuron.options.merge_thr = obj.Parameters.CNMFE.merge_thr;  % thresholds for merging neurons; [spatial overlap ratio, temporal correlation of calcium traces, spike correlation]
                        Neuron.options.method_dist = 'mean';                        % method for computing neuron distances {'mean', 'max'}

                        % -------------------------  INITIALIZATION   -------------------------  %
                        Neuron.options.bd = 1;                                  % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
                        Neuron.options.center_psf = true;                       % set the value as true when the background fluctuation is large (usually 1p data)
                        Neuron.options.min_corr = obj.Parameters.CNMFE.min_corr;% minimum local correlation for a seeding pixel
                        Neuron.options.min_pnr = obj.Parameters.CNMFE.min_pnr;  % minimum peak-to-noise ratio for a seeding pixel
                        % Place the Sources2D object in our own object
                        obj.Parameters.CNMFE.Neuron = Neuron;

                        % MP parameters: inherited from the main parameters
                        % by default
                        obj.Parameters.CNMFE.MaxProjection_Sig = obj.Parameters.CNMFE.gSig;
                        obj.Parameters.CNMFE.MaxProjection_Sub = obj.Parameters.CNMFE.CorrSubsampling;
                    end
                end

                % Wipe current layout
                obj.WipeLayout;

                % The display here is a bit different: we only need the
                % correlation and PNR plots (and the product), and the UI
                % to set the different parameters' values

                % Initialize empty plots (long processing time, we'd better
                % be sure of the parameters we want to use)

                obj.Handles.UIElements.Cn_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.195 0.9025 0.2 0.02],...
                    'String','Correlation (Cn)',...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center');
                obj.Handles.UIElements.PNR_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.195+0.205 0.9025 0.2 0.02],...
                    'String','Peak-to-noise (PNR)',...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center');
                obj.Handles.UIElements.Product_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.195+0.205*2 0.9025 0.2 0.02],...
                    'String','Product',...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center');
                obj.Handles.UIElements.Parameters_1 = axes(...
                    'Units','normalized',...
                    'Position',[0.195+0.205*3 0.65 0.155 0.25],...
                    'Color',[0.2 0.2 0.2],'XColor','k','YColor','k','Box','on','LineWidth',1.5);
                obj.Handles.UIElements.Parameters_1.XTick = [];
                obj.Handles.UIElements.Parameters_1.YTick = [];
                obj.Handles.UIElements.Parameters_2 = axes(...
                    'Units','normalized',...
                    'Position',[0.195+0.205*3 0.39 0.155 0.25],...
                    'Color',[0.2 0.2 0.2],'XColor','k','YColor','k','Box','on','LineWidth',1.5);
                obj.Handles.UIElements.Parameters_2.XTick = [];
                obj.Handles.UIElements.Parameters_2.YTick = [];

                % UI to process the Correlation/PNR
                Shifts = [0,0.26];
                if ~isfield(obj.Handles,'Values') || ~isfield(obj.Handles.Values,'Sig')
                    obj.Handles.Values.Sig([1 2]) = obj.Parameters.CNMFE.gSig;
                    obj.Handles.Values.CorrSubsampling([1 2]) = obj.Parameters.CNMFE.CorrSubsampling;
                    obj.Handles.Values.min_corr([1 2]) = obj.Parameters.CNMFE.min_corr;
                    obj.Handles.Values.min_pnr([1 2]) = obj.Parameters.CNMFE.min_pnr;
                    obj.Handles.Values.File{1} = 'PreProcessed';
                    obj.Handles.Values.File{2} = 'MotionCorrected';
                    obj.Handles.Values.PNR = cell(2,1);
                    obj.Handles.Values.Cn = cell(2,1);
                    obj.Handles.Values.MaxProjection = [];
                    obj.Handles.Values.ProjectionSig = obj.Parameters.CNMFE.gSig;
                    obj.Handles.Values.ProjectionSubsampling = obj.Parameters.CNMFE.CorrSubsampling;
                    obj.Handles.Values.ProjectionCLim = [0 1];
                    obj.Handles.Values.Binarize([1 2]) = false;
                    obj.Handles.Values.SetOnce = false;
                end
                for CS = 1 : 2,
                    obj.Handles.UIElements.Cn(CS) = axes(...
                        'Units','normalized',...
                        'Position',[0.195 0.65-Shifts(CS) 0.2 0.25],...
                        'Color','k','XColor','none','YColor','none');hold on
                    if isfield(obj.Handles,'Values') && isfield(obj.Handles.Values,'Cn')
                        imagesc(obj.Handles.Values.Cn{CS},'Parent',obj.Handles.UIElements.Cn(CS))
                    end
                    obj.Handles.UIElements.PNR(CS) = axes(...
                        'Units','normalized',...
                        'Position',[0.195+0.205 0.65-Shifts(CS) 0.2 0.25],...
                        'Color','k','XColor','none','YColor','none');hold on
                    if isfield(obj.Handles,'Values') && isfield(obj.Handles.Values,'PNR')
                        imagesc(obj.Handles.Values.PNR{CS},'Parent',obj.Handles.UIElements.PNR(CS))
                    end
                    obj.Handles.UIElements.Product(CS) = axes(...
                        'Units','normalized',...
                        'Position',[0.195+0.205*2 0.65-Shifts(CS) 0.2 0.25],...
                        'Color','k','XColor','none','YColor','none');hold on
                    %Greedy initialization as follow in CNMFE code:
                    % v_search = Cn.*PNR;
                    % v_search(or(Cn<min_corr, PNR<min_pnr)) = 0;
                    Prod = [];
                    if isfield(obj.Handles,'Values') && isfield(obj.Handles.Values,'Cn')&& isfield(obj.Handles.Values,'PNR')
                        Prod = obj.Handles.Values.Cn{CS}.*obj.Handles.Values.PNR{CS};
                        Prod(obj.Handles.Values.Cn{CS}<obj.Handles.Values.min_corr(CS) | obj.Handles.Values.PNR{CS}<obj.Handles.Values.min_pnr(CS)) = 0;
                        imagesc(Prod,'Parent',obj.Handles.UIElements.Product(CS))
                    end
                    if isfield(obj.Handles,'Values') && isfield(obj.Handles.Values,'Cn') && ~isempty(obj.Handles.Values.Cn{CS})
                        obj.Handles.UIElements.Cn(CS).CLim = [obj.Handles.Values.min_corr(CS) 1];
                        if obj.Handles.Values.min_pnr(CS) < 0.98*max(obj.Handles.Values.PNR{CS},[],'all')
                            obj.Handles.UIElements.PNR(CS).CLim = [obj.Handles.Values.min_pnr(CS) 0.98*max(obj.Handles.Values.PNR{CS},[],'all')];
                        else
                            obj.Handles.UIElements.PNR(CS).CLim = [0.98*max(obj.Handles.Values.PNR{CS},[],'all')-0.01 0.98*max(obj.Handles.Values.PNR{CS},[],'all')];
                            obj.Handles.Values.min_pnr(CS) = 0.98*max(obj.Handles.Values.PNR{CS},[],'all')-0.01;
                        end
                        obj.Handles.UIElements.Product(CS).CLim = [obj.Handles.Values.min_corr(CS)*obj.Handles.Values.min_pnr(CS) 0.98*max(obj.Handles.Values.PNR{CS},[],'all')];
                        if obj.Handles.Values.Binarize(CS),
                            obj.Handles.UIElements.PNR(CS).CLim = obj.Handles.UIElements.PNR(CS).CLim(1)*[1 1+eps];
                            obj.Handles.UIElements.Cn(CS).CLim = obj.Handles.UIElements.Cn(CS).CLim(1)*[1 1+eps];
                            obj.Handles.UIElements.Product(CS).CLim = obj.Handles.UIElements.Product(CS).CLim(1)*[1 1+eps];
                        end

                        % Adjust limits to preserve ratio
                        [D1,D2] = size(obj.Handles.Values.Cn{CS});
                        obj.Handles.UIElements.Cn(CS).Units = 'pixels';
                        AbsolutePosition = obj.Handles.UIElements.Cn(CS).Position;
                        obj.Handles.UIElements.Cn(CS).Units = 'normalized';
                        if D2/D1 >= AbsolutePosition(3)/AbsolutePosition(4),
                            % Y needs to be adjusted
                            YDelta = D2*AbsolutePosition(4)/AbsolutePosition(3) - D1;
                            obj.Handles.UIElements.Cn(CS).YLim = [1-0.5*YDelta D1+0.5*YDelta];
                            obj.Handles.UIElements.Cn(CS).XLim = [1 D2];
                            obj.Handles.UIElements.PNR(CS).YLim = [1-0.5*YDelta D1+0.5*YDelta];
                            obj.Handles.UIElements.PNR(CS).XLim = [1 D2];
                            obj.Handles.UIElements.Product(CS).YLim = [1-0.5*YDelta D1+0.5*YDelta];
                            obj.Handles.UIElements.Product(CS).XLim = [1 D2];
                        else
                            % X needs to be adjusted
                            XDelta = D1*AbsolutePosition(3)/AbsolutePosition(4) - D2;
                            obj.Handles.UIElements.Cn(CS).XLim = [1-0.5*XDelta D2+0.5*XDelta];
                            obj.Handles.UIElements.Cn(CS).YLim = [1 D1];
                            obj.Handles.UIElements.PNR(CS).XLim = [1-0.5*XDelta D2+0.5*XDelta];
                            obj.Handles.UIElements.PNR(CS).YLim = [1 D1];
                            obj.Handles.UIElements.Product(CS).XLim = [1-0.5*XDelta D2+0.5*XDelta];
                            obj.Handles.UIElements.Product(CS).YLim = [1 D1];
                        end
                    end

                    obj.Handles.UIElements.Config_Legend(CS) = uicontrol('Style','text',...
                        'Units','Normalized','Position',[0.195+0.205*3 0.88-Shifts(CS) 0.1555 0.02],...
                        'String',['Config Test #' num2str(CS)],...
                        'FontSize',13,'FontName','Arial','FontWeight','b',...
                        'BackgroundColor','k','ForegroundColor',[0.6 0.6 0.6],...
                        'HorizontalAlignment','center');
                    % Data selection
                    obj.Handles.UIElements.File_Legend(CS) = uicontrol('Style','text',...
                        'Units','Normalized','Position',[0.195+0.205*3+0.005 0.83-Shifts(CS) 0.05 0.03],...
                        'String','Data',...
                        'FontSize',13,'FontName','Arial','FontWeight','b',...
                        'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                        'HorizontalAlignment','left');
                    obj.Handles.UIElements.File_Preprocessed_Button(CS) = uicontrol('Style','togglebutton',...
                        'Units','Normalized','Position',[0.195+0.205*3+0.005+0.025 0.835-Shifts(CS) 0.035 0.03],...
                        'String','Raw','Value',0,...
                        'FontSize',13,'FontName','Arial','FontWeight','b',...
                        'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                        'HorizontalAlignment','center',...
                        'UserData',num2str(CS),...
                        'Value',0,...
                        'Callback',@(src,evt)obj.CNMFE_Preprocessed_CB(src,evt));
                    obj.Handles.UIElements.File_MotionCorrected_Button(CS) = uicontrol('Style','togglebutton',...
                        'Units','Normalized','Position',[0.195+0.205*3+0.005+0.025+0.035 0.835-Shifts(CS) 0.0825 0.03],...
                        'String','Motion corrected','Value',0,...
                        'FontSize',13,'FontName','Arial','FontWeight','b',...
                        'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                        'HorizontalAlignment','center',...
                        'UserData',num2str(CS),...
                        'Value',0,...
                        'Callback',@(src,evt)obj.CNMFE_MotionCorrected_CB(src,evt));
                    if isfield(obj.Handles,'Values') && isfield(obj.Handles.Values,'File')
                        if strcmpi(obj.Handles.Values.File{CS},'MotionCorrected'),
                            obj.Handles.UIElements.File_Preprocessed_Button(CS).Value = 0;
                            obj.Handles.UIElements.File_MotionCorrected_Button(CS).Value = 1;
                        else
                            obj.Handles.UIElements.File_Preprocessed_Button(CS).Value = 1;
                            obj.Handles.UIElements.File_MotionCorrected_Button(CS).Value = 0;
                        end
                    end
                    obj.Handles.UIElements.Sig_Legend(CS) = uicontrol('Style','text',...
                        'Units','Normalized','Position',[0.195+0.205*3+0.005 0.8-Shifts(CS) 0.0575 0.025],...
                        'String','Sigma',...
                        'FontSize',13,'FontName','Arial','FontWeight','b',...
                        'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                        'HorizontalAlignment','right');
                    obj.Handles.UIElements.Sig_Edit(CS) = uicontrol('Style','edit',...
                        'Units','Normalized','Position',[0.195+0.205*3+0.005+0.025+0.035 0.8025-Shifts(CS) 0.0825/2 0.025],...
                        'String', num2str(obj.Handles.Values.Sig(CS)),'Value',0,...
                        'FontSize',13,'FontName','Arial','FontWeight','b',...
                        'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                        'HorizontalAlignment','center',...
                        'UserData',num2str(CS),...
                        'Callback',@(src,evt)obj.CNMFE_Sig_Edit_CB(src,evt));

                    obj.Handles.UIElements.Sub_Legend(CS) = uicontrol('Style','text',...
                        'Units','Normalized','Position',[0.195+0.205*3+0.005 0.77-Shifts(CS) 0.0575 0.025],...
                        'String','Subsampling',...
                        'FontSize',13,'FontName','Arial','FontWeight','b',...
                        'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                        'HorizontalAlignment','right');
                    obj.Handles.UIElements.Sub_Edit(CS) = uicontrol('Style','edit',...
                        'Units','Normalized','Position',[0.195+0.205*3+0.005+0.025+0.035 0.7725-Shifts(CS) 0.0825/2 0.025],...
                        'String',num2str(obj.Handles.Values.CorrSubsampling(CS)),'Value',0,...
                        'FontSize',13,'FontName','Arial','FontWeight','b',...
                        'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                        'HorizontalAlignment','center',...
                        'UserData',num2str(CS),...
                        'Callback',@(src,evt)obj.CNMFE_Sub_Edit_CB(src,evt));

                    obj.Handles.UIElements.Binarize_Checkbox(CS) = uicontrol('Style','checkbox',...
                        'Units','Normalized','Position',[0.195+0.205*3+0.005 0.735-Shifts(CS) 0.0575 0.03],...
                        'String',' Binarize',...
                        'FontSize',13,'FontName','Arial','FontWeight','b',...
                        'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                        'HorizontalAlignment','left',...
                        'UserData',num2str(CS),...
                        'Value',obj.Handles.Values.Binarize(CS),...
                        'Callback',@(src,evt)obj.Binarize_CB(src,evt));

                    obj.Handles.UIElements.Apply_Button(CS) = uicontrol('Style','pushbutton',...
                        'Units','Normalized','Position',[0.195+0.205*3+0.005+0.025+0.035 0.735-Shifts(CS) 0.0825/2 0.03],...
                        'String','Apply',...
                        'FontSize',13,'FontName','Arial','FontWeight','b',...
                        'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                        'HorizontalAlignment','right',...
                        'UserData',num2str(CS),...
                        'Callback',@(src,evt)obj.Process_Correlation(src,evt));
                    obj.Handles.UIElements.Set_Button(CS) = uicontrol('Style','pushbutton',...
                        'Units','Normalized','Position',[0.195+0.205*3+0.005+0.025+0.035+0.0825/2 0.735-Shifts(CS) 0.0825/2 0.03],...
                        'String','Set',...
                        'FontSize',13,'FontName','Arial','FontWeight','b',...
                        'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                        'UserData',num2str(CS),...
                        'Enable','off',...
                        'Callback',@(src,evt)obj.CNMFE_Set_CB(src,evt));

                    obj.Handles.UIElements.CnSlider_Legend(CS) = uicontrol('Style','text',...
                        'Units','Normalized','Position',[0.195+0.205*3+0.005 0.695-Shifts(CS) 0.05 0.025],...
                        'String','Cn',...
                        'FontSize',13,'FontName','Arial','FontWeight','b',...
                        'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                        'HorizontalAlignment','left');
                    obj.Handles.UIElements.CnSlider_Edit(CS) = uicontrol('Style','edit',...
                        'Units','Normalized','Position',[0.195+0.205*3+0.005+0.025 0.6975-Shifts(CS) 0.035 0.025],...
                        'String',num2str(obj.Handles.Values.min_corr(CS),'%.2f'),...
                        'FontSize',13,'FontName','Arial','FontWeight','b',...
                        'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                        'HorizontalAlignment','center',...
                        'UserData',num2str(CS),...
                        'Callback',@(src,evt)obj.CnSlider_Edit_CB(src,evt));
                    obj.Handles.UIElements.PNRSlider_Legend(CS) = uicontrol('Style','text',...
                        'Units','Normalized','Position',[0.195+0.205*3+0.005 0.665-Shifts(CS) 0.05 0.025],...
                        'String','PNR',...
                        'FontSize',13,'FontName','Arial','FontWeight','b',...
                        'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                        'HorizontalAlignment','left');
                    obj.Handles.UIElements.PNRSlider_Edit(CS) = uicontrol('Style','edit',...
                        'Units','Normalized','Position',[0.195+0.205*3+0.005+0.025 0.6675-Shifts(CS) 0.035 0.025],...
                        'String',num2str(obj.Handles.Values.min_pnr(CS),'%.2f'),...
                        'FontSize',13,'FontName','Arial','FontWeight','b',...
                        'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                        'HorizontalAlignment','center',...
                        'UserData',num2str(CS),...
                        'Callback',@(src,evt)obj.PNRSlider_Edit_CB(src,evt));

                    obj.Handles.UIElements.CnSlider_Axes(CS) = axes('Units','normalized',...
                        'Position',[0.195+0.205*3+0.005+0.025+0.035+0.005 0.695-Shifts(CS) 0.0775 0.03],...
                        'Color',[0.2,0.2,0.2],'LineWidth',1); hold on;
                    obj.Handles.UIElements.PNRSlider_Axes(CS) = axes('Units','normalized',...
                        'Position',[0.195+0.205*3+0.005+0.025+0.035+0.005 0.665-Shifts(CS) 0.0775 0.03],...
                        'Color',[0.2,0.2,0.2],'LineWidth',1); hold on

                    obj.Handles.UIElements.PNRSlider_BaseLine(CS) = plot([0 20],[0 0],...
                        'Parent',obj.Handles.UIElements.PNRSlider_Axes(CS),...
                        'LineWidth',1.5,'Color','k');
                    obj.Handles.UIElements.PNRSlider_Slider(CS) = plot(obj.Handles.Values.min_pnr(CS)*[1 1],[-1 1],...
                        'Parent',obj.Handles.UIElements.PNRSlider_Axes(CS),...
                        'LineWidth',3,'Color','k',...
                        'UserData',num2str(CS),...
                        'ButtonDownFcn',{@(src,evt)obj.PNRSliderCB(src,evt)});
                    obj.Handles.UIElements.PNRSlider_Axes(CS).XAxis.Visible = 'off';
                    obj.Handles.UIElements.PNRSlider_Axes(CS).YAxis.Visible = 'off';

                    if ~isempty(Prod)
                        obj.Handles.UIElements.PNRSlider_Axes(CS).XLim = [0 max(obj.Handles.Values.PNR{CS},[],'all')];
                        obj.Handles.UIElements.PNRSlider_BaseLine(CS).XData(2) = max(obj.Handles.Values.PNR{CS},[],'all');
                    else
                        obj.Handles.UIElements.PNRSlider_Axes(CS).XLim = [0 20];
                    end
                    obj.Handles.UIElements.PNRSlider_Axes(CS).YLim = [-1.5 1.5];

                    obj.Handles.UIElements.CnSlider_BaseLine(CS) = plot([0 1],[0 0],...
                        'Parent',obj.Handles.UIElements.CnSlider_Axes(CS),...
                        'LineWidth',1.5,'Color','k');
                    obj.Handles.UIElements.CnSlider_Slider(CS) = plot(obj.Handles.Values.min_corr(CS)*[1 1],[-1 1],...
                        'Parent',obj.Handles.UIElements.CnSlider_Axes(CS),...
                        'LineWidth',3,'Color','k',...
                        'UserData',num2str(CS),...
                        'ButtonDownFcn',{@(src,evt)obj.CnSliderCB(src,evt)});
                    obj.Handles.UIElements.CnSlider_Axes(CS).XAxis.Visible = 'off';
                    obj.Handles.UIElements.CnSlider_Axes(CS).YAxis.Visible = 'off';
                    obj.Handles.UIElements.CnSlider_Axes(CS).XLim = [0 1];
                    obj.Handles.UIElements.CnSlider_Axes(CS).YLim = [-1.5 1.5];
                end

                % Maximum projection to set neuron size and cropping
                obj.Handles.UIElements.MovieDisplay(3).Plot = axes(...
                    'Units','normalized',...
                    'Position',[0.195 0.04 0.256 0.32],...
                    'Color','k','XColor','none','YColor','none',...
                    'Colormap',bone);hold on
                obj.Handles.UIElements.MovieDisplay(3).MaxProjection = imagesc(obj.Handles.Values.MaxProjection,'Parent',obj.Handles.UIElements.MovieDisplay(3).Plot);
                obj.Handles.UIElements.MovieDisplay(3).Plot.CLim = obj.Handles.Values.ProjectionCLim;

                if ~isempty(obj.Handles.Values.MaxProjection),
                    % Adjust limits to preserve ratio
                    [D1,D2] = size(obj.Handles.Values.MaxProjection);
                    obj.Handles.UIElements.MovieDisplay(3).Plot.Units = 'pixels';
                    AbsolutePosition = obj.Handles.UIElements.MovieDisplay(3).Plot.Position;
                    obj.Handles.UIElements.MovieDisplay(3).Plot.Units = 'normalized';
                    if D2/D1 >= AbsolutePosition(3)/AbsolutePosition(4),
                        % Y needs to be adjusted
                        YDelta = D2*AbsolutePosition(4)/AbsolutePosition(3) - D1;
                        obj.Handles.UIElements.MovieDisplay(3).YLim = [1-0.5*YDelta D1+0.5*YDelta];
                        obj.Handles.UIElements.MovieDisplay(3).XLim = [1 D2];
                    else
                        % X needs to be adjusted
                        XDelta = D1*AbsolutePosition(3)/AbsolutePosition(4) - D2;
                        obj.Handles.UIElements.MovieDisplay(3).Plot.XLim = [1-0.5*XDelta D2+0.5*XDelta];
                        obj.Handles.UIElements.MovieDisplay(3).Plot.YLim = [1 D1];
                    end
                end
                if ~isempty(obj.Parameters.CNMFE.Cell),
                    Center = obj.Parameters.CNMFE.Cell;
                else
                    Center = [0.5*diff(obj.Handles.UIElements.MovieDisplay(3).Plot.XLim) + obj.Handles.UIElements.MovieDisplay(3).Plot.XLim(1),...
                        0.5*diff(obj.Handles.UIElements.MovieDisplay(3).Plot.YLim) + obj.Handles.UIElements.MovieDisplay(3).Plot.YLim(1)];
                    obj.Parameters.CNMFE.Cell = Center;
                end
                % Circle equation for the ROI
                Theta = 0:pi/50:2*pi;
                XCoor = obj.Parameters.CNMFE.gSiz/2 * cos(Theta) + Center(1);
                YCoor = obj.Parameters.CNMFE.gSiz/2 * sin(Theta) + Center(2);
                obj.Handles.UIElements.NeuronROI = plot(XCoor,YCoor,'LineWidth',2,'Color','g',...
                    'Parent',obj.Handles.UIElements.MovieDisplay(3).Plot,...
                    'ButtonDownFcn', @(~,~)obj.CNMFE_ProjectionMoveNeuron_CB);

                obj.Handles.UIElements.MovieDisplay(3).ParametersBox = axes(...
                    'Units','normalized',...
                    'Position',[0.195+obj.Handles.UIElements.MovieDisplay(3).Plot.Position(3)+0.005 0.04 0.144 0.32],...
                    'Color',[0.2 0.2 0.2],'XColor','k','YColor','k','Box','on','LineWidth',1.5);
                obj.Handles.UIElements.MovieDisplay(3).ParametersBox.XTick = [];
                obj.Handles.UIElements.MovieDisplay(3).ParametersBox.YTick = [];
                obj.Handles.UIElements.MaxProjection_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(1) obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(2)+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(4)-0.02 obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(3)+0.0005 0.02],...
                    'String','Spatial parameters',...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor','k','ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center');

                obj.Handles.UIElements.ProjectionSig_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.0075+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(1) obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(2)+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(4)-0.08 0.055 0.025],...
                    'String','Sigma',...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','right');
                obj.Handles.UIElements.ProjectionSig_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(1)+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(3)/2 obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(2)+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(4)-0.08 0.0825/2 0.025],...
                    'String', num2str(obj.Handles.Values.ProjectionSig),'Value',0,...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',@(~,~)obj.CNMFE_ProjectionSig_Edit_CB);

                obj.Handles.UIElements.ProjectionSub_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.0075+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(1) obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(2)+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(4)-0.11 0.055 0.025],...
                    'String','Subsampling',...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','right');
                obj.Handles.UIElements.ProjectionSub_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(1)+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(3)/2 obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(2)+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(4)-0.11 0.0825/2 0.025],...
                    'String',num2str(obj.Handles.Values.ProjectionSubsampling),'Value',0,...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',@(~,~)obj.CNMFE_ProjectionSub_Edit_CB);

                obj.Handles.UIElements.ProcessProjection_Button = uicontrol('Style','pushbutton',...
                    'Units','Normalized','Position',[obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(1)+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(3)/2 obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(2)+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(4)-0.155 0.0825/1.5 0.04],...
                    'String','Process',...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','right',...
                    'Callback',@(~,~)obj.Process_MaxProjection);

                obj.Handles.UIElements.ProjectionCLim_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.0075+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(1) obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(2)+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(4)-0.205 0.055/2 0.025],...
                    'String','CLim',...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','right');

                obj.Handles.UIElements.CLimSlider_Axes = axes('Units','normalized',...
                    'Position',[0.0075 + 0.005 + obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(1) + 0.055/2,...
                    obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(2)+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(4)-0.205,...
                    0.0870 0.03],...
                    'Color',[0.2,0.2,0.2],'LineWidth',1); hold on;
                obj.Handles.UIElements.CLimSlider_BaseLine = plot([0 1],[0 0],...
                    'Parent',obj.Handles.UIElements.CLimSlider_Axes,...
                    'LineWidth',1.5,'Color','k');
                obj.Handles.UIElements.CLimSlider_SliderLow = plot(obj.Handles.Values.ProjectionCLim(1)*[1 1],[-1 1],...
                    'Parent',obj.Handles.UIElements.CLimSlider_Axes,...
                    'LineWidth',3,'Color','k',...
                    'UserData',num2str(CS),...
                    'ButtonDownFcn',{@(~,~)obj.CLimSliderLowCB});
                obj.Handles.UIElements.CLimSlider_SliderHigh = plot(obj.Handles.Values.ProjectionCLim(2)*[1 1],[-1 1],...
                    'Parent',obj.Handles.UIElements.CLimSlider_Axes,...
                    'LineWidth',3,'Color','k',...
                    'UserData',num2str(CS),...
                    'ButtonDownFcn',{@(~,~)obj.CLimSliderHighCB});

                obj.Handles.UIElements.CLimSlider_Axes.XAxis.Visible = 'off';
                obj.Handles.UIElements.CLimSlider_Axes.YAxis.Visible = 'off';
                obj.Handles.UIElements.CLimSlider_Axes.XLim = [0 1];
                obj.Handles.UIElements.CLimSlider_Axes.YLim = [-1.5 1.5];

                obj.Handles.UIElements.ProjectionCell_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.0075+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(1) obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(2)+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(4)-0.255 0.035 0.025],...
                    'String','Cell size',...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','right');
                obj.Handles.UIElements.ProjectionCellSize_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.0075+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(1)+0.035+0.005 obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(2)+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(4)-0.255+0.0025 0.0795/2 0.025],...
                    'String', num2str(obj.Parameters.CNMFE.gSiz),...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',@(~,~)obj.ProjectionCellSize_Edit_CB);
                obj.Handles.UIElements.ProjectionCellSize_Set = uicontrol('Style','pushbutton',...
                    'Units','Normalized','Position',[0.0795/2+0.0075+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(1)+0.035+0.005 obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(2)+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(4)-0.255+0.0025 0.0795/2 0.025],...
                    'String', 'Set',...
                    'Enable','on',...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center'); % Just to trigger the callback from the edit box, no need to do anything else
                obj.Handles.UIElements.ProjectionCrop_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.0075+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(1) obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(2)+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(4)-0.295 0.035 0.025],...
                    'String','Crop',...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','right');
                obj.Handles.UIElements.ProjectionCrop_Draw = uicontrol('Style','pushbutton',...
                    'Units','Normalized','Position',[0.0075+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(1)+0.035+0.005 obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(2)+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(4)-0.295+0.0025 0.0795/2 0.025],...
                    'String', 'Draw','Value',0,...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',@(~,~)obj.CNMFE_ProjectionDrawCrop_CB);
                obj.Handles.UIElements.ProjectionCrop_Delete = uicontrol('Style','pushbutton',...
                    'Units','Normalized','Position',[0.0795/2+0.0075+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(1)+0.035+0.005 obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(2)+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(4)-0.295+0.0025 0.0795/2 0.025],...
                    'String', 'Delete','Value',0,...
                    'Enable','off',...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',@(~,~)obj.CNMFE_ProjectionDeleteCrop_CB);

                if ~isempty(obj.Parameters.CNMFE.Crop)
                    if numel(obj.Parameters.CNMFE.Crop)==4
                        obj.Handles.UIElements.Crop = drawrectangle(obj.Handles.UIElements.MovieDisplay(3).Plot,...
                            'Position',obj.Parameters.CNMFE.Crop,...
                            'FaceAlpha',0,...
                            'LineWidth',3,...
                            'Color',[0.6 0.3 0],...
                            'Deletable',false,...
                            'FaceSelectable',false);
                        addlistener(obj.Handles.UIElements.Crop(:),'MovingROI',@(~,~)obj.CNMFE_ProjectionEditCrop_CB);
                        obj.Handles.UIElements.ProjectionCrop_Delete.Enable = 'on';
                        obj.Handles.UIElements.ProjectionCrop_Draw.Enable = 'off';
                        drawnow
                    else
                        obj.Parameters.CNMFE.Crop = [];
                    end
                end


                % All parameters
                obj.Handles.UIElements.AllParametersBox = axes(...
                    'Units','normalized',...
                    'Position',[obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(1)+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(3)+0.005 0.04 0.36 0.32],...
                    'Color',[0.2 0.2 0.2],'XColor','k','YColor','k','Box','on','LineWidth',1.5);
                obj.Handles.UIElements.AllParametersBox.XTick = [];
                obj.Handles.UIElements.AllParametersBox.YTick = [];
                obj.Handles.UIElements.AllParameters_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[obj.Handles.UIElements.AllParametersBox.Position(1) obj.Handles.UIElements.AllParametersBox.Position(2)+obj.Handles.UIElements.AllParametersBox.Position(4)-0.02 obj.Handles.UIElements.AllParametersBox.Position(3)+0.0005 0.02],...
                    'String','CNMF-E parameters',...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor','k','ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center');

                % Spatial
                obj.Handles.UIElements.SpatialCNMFE_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.61+0.015 0.31 0.075 0.025],...
                    'String','Spatial',...
                    'FontSize',15,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.SigmaCNMFE_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.6105+0.015, 0.2750, 0.0450, 0.0250],...
                    'String','Sigma',...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','right');
                obj.Handles.UIElements.SigmaCNMFE_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.6105+0.015+0.005+0.0450, 0.2750, 0.0413, 0.0250],...
                    'String', num2str(obj.Parameters.CNMFE.gSig),...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Enable','inactive');

                obj.Handles.UIElements.NeuronCNMFE_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.6105+0.015, 0.2450, 0.0450, 0.0250],...
                    'String','Cell size',...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','right');
                obj.Handles.UIElements.NeuronCNMFE_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.6105+0.015+0.005+0.0450, 0.2450, 0.0413, 0.0250],...
                    'String', num2str(obj.Parameters.CNMFE.gSiz),...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Enable','inactive');

                obj.Handles.UIElements.MergingCNMFE_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.6105+0.015, 0.2000, 0.0450, 0.0250],...
                    'String','Merging',...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','right');
                obj.Handles.UIElements.MergingThreshSpatialCNMFE_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.6105+0.015+0.005+0.0450, 0.2150, 0.0413, 0.0250],...
                    'String', num2str(obj.Parameters.CNMFE.merge_thr_spatial(1)),...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Tooltip','Spatial correlation threshold (0-1)',...
                    'Callback',@(~,~)obj.MergingThreshSpatialCNMFE_Edit_CB);
                obj.Handles.UIElements.MergingThreshTemporalCNMFE_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.6105+0.015+0.005+0.0450+0.0413, 0.2150, 0.0413, 0.0250],...
                    'String', num2str(obj.Parameters.CNMFE.merge_thr_spatial(2)),...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Tooltip','Temporal correlation threshold (0-1)',...
                    'Callback',@(~,~)obj.MergingThreshTemporalCNMFE_Edit_CB);
                obj.Handles.UIElements.MergingValueCNMFE_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.6105+0.015+0.005+0.0450, 0.185, 0.0413, 0.0250],...
                    'String', num2str(obj.Parameters.CNMFE.merge_thr),...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Tooltip',sprintf('thresholds for merging neurons.\n[spatial overlap ratio, temporal correlation of calcium traces, spike correlation]'),...
                    'Callback',@(~,~)obj.MergingValueCNMFE_Edit_CB);
                obj.Handles.UIElements.RingRadiusCNMFE_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.610, 0.1550, 0.0450+0.015, 0.0250],...
                    'String','Ring radius',...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','right');
                obj.Handles.UIElements.RingRadiusCNMFE_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.6105+0.015+0.005+0.0450, 0.1550, 0.0413, 0.0250],...
                    'String', num2str(obj.Parameters.CNMFE.Neuron.options.ring_radius),...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Tooltip','Ring radius for background estimation; typically slightly larger than cell radius.',...
                    'Callback',@(~,~)obj.RingRadiusCNMFE_Edit_CB);


                % Temporal
                obj.Handles.UIElements.TemporalCNMFE_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.61+0.015 0.115 0.075 0.025],...
                    'String','Temporal',...
                    'FontSize',15,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','left');

                obj.Handles.UIElements.DecayCNMFE_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.6105+0.015, 0.08, 0.0450, 0.0250],...
                    'String','Max decay',...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','right');
                obj.Handles.UIElements.DecayCNMFE_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.6105+0.015+0.005+0.0450, 0.08, 0.0413, 0.0250],...
                    'String', num2str(obj.Parameters.CNMFE.max_tau),...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Tooltip','Maximum decay time in seconds.',...
                    'Callback',@(~,~)obj.DecayCNMFE_Edit_CB);

                obj.Handles.UIElements.DetrendCNMFE_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.6105+0.015, 0.05, 0.0450, 0.0250],...
                    'String','Detrend',...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','right');

                obj.Handles.UIElements.DetrendCNMFE_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.6105+0.015+0.005+0.0450, 0.05, 0.0413, 0.0250],...
                    'String', num2str(obj.Parameters.CNMFE.Neuron.options.nk),...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Tooltip','Factor to detrend the signal (<1/30 total time).',...
                    'Callback',@(~,~)obj.DetrendCNMFE_Edit_CB);

                % Initialization
                obj.Handles.UIElements.InitializationCNMFE_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.79+0.015 0.31 0.075 0.025],...
                    'String','Initialization',...
                    'FontSize',15,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','left');
                obj.Handles.UIElements.mincorrCNMFE_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.6105+0.015+0.18, 0.2750, 0.0450*2, 0.0250],...
                    'String','Min correlation',...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','right');
                obj.Handles.UIElements.mincorrCNMFE_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.6105+0.015+0.18+0.005+0.0450*2, 0.2750, 0.0413, 0.0250],...
                    'String', num2str(obj.Parameters.CNMFE.Neuron.options.min_corr),...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Enable','inactive');

                obj.Handles.UIElements.minPNRCNMFE_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.6105+0.015+0.18, 0.2450, 0.0450*2, 0.0250],...
                    'String','Min peak-to-noise',...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','right');
                obj.Handles.UIElements.minPNRCNMFE_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.6105+0.015+0.18+0.005+0.0450*2, 0.2450, 0.0413, 0.0250],...
                    'String', num2str(obj.Parameters.CNMFE.Neuron.options.min_pnr),...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Enable','inactive');

                obj.Handles.UIElements.maxRAMCNMFE_Legend = uicontrol('Style','text',...
                    'Units','Normalized','Position',[0.6105+0.015+0.18, 0.2150, 0.0450*2, 0.0250],...
                    'String','Max. RAM',...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','right');
                % Legacy
                if ~isfield(obj.Parameters.CNMFE,'maxRAM')
                    obj.Parameters.CNMFE.maxRAM = 64;
                end
                obj.Handles.UIElements.maxRAMCNMFE_Edit = uicontrol('Style','edit',...
                    'Units','Normalized','Position',[0.6105+0.015+0.18+0.005+0.0450*2, 0.2150, 0.0413, 0.0250],...
                    'String', num2str(obj.Parameters.CNMFE.maxRAM),...
                    'FontSize',13,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Tooltip','Maximum RAM that CNMFE is allowed to recruit.',...
                    'Callback',@(~,~)obj.maxRAMCNMFE_Edit_CB);


                % Validate
                if obj.Handles.Values.SetOnce, % Loose way of disabling direct processing without checking parameters
                    EnV = 'on';
                else
                    EnV = 'off';
                end
                obj.Handles.UIElements.ValidateCNMFE_Button = uicontrol('Style','pushbutton',...
                    'Units','Normalized','Position',[0.9418-0.075 0.07 0.075 0.04],...
                    'String','Process',...
                    'FontSize',14,'FontName','Arial','FontWeight','b',...
                    'BackgroundColor',[0.15 0.15 0.15],'ForegroundColor',[0.6 0.6 0.6],...
                    'HorizontalAlignment','center',...
                    'Callback',{@(~,~)obj.CNMFE},...
                    'Enable',EnV);

            else
                % Update (e.g. if we load another session while in the tab)


            end

        end

        function MergingThreshSpatialCNMFE_Edit_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.MergingThreshSpatialCNMFE_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>0 && TempInput<1,
                obj.Parameters.CNMFE.merge_thr_spatial(1) = TempInput;
            else
                obj.Handles.UIElements.MergingThreshSpatialCNMFE_Edit.String = num2str(obj.Parameters.CNMFE.merge_thr_spatial(1));
            end
        end

        function MergingThreshTemporalCNMFE_Edit_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.MergingThreshTemporalCNMFE_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>0 && TempInput<1,
                obj.Parameters.CNMFE.merge_thr_spatial(2) = TempInput;
            else
                obj.Handles.UIElements.MergingThreshTemporalCNMFE_Edit.String = num2str(obj.Parameters.CNMFE.merge_thr_spatial(2));
            end
        end

        function MergingValueCNMFE_Edit_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.MergingValueCNMFE_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>0 && TempInput<1,
                obj.Parameters.CNMFE.merge_thr = TempInput;
                obj.Parameters.CNMFE.Neuron.options.merge_thr = TempInput;
            else
                obj.Handles.UIElements.MergingValueCNMFE_Edit.String = num2str(obj.Parameters.CNMFE.merge_thr);
            end
        end

        function CNMFE_ProjectionSig_Edit_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.ProjectionSig_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>0,
                obj.Handles.Values.ProjectionSig = TempInput;
            else
                obj.Handles.UIElements.ProjectionSig_Edit.String = num2str(obj.Handles.Values.ProjectionSig);
            end
        end

        function CNMFE_ProjectionSub_Edit_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.ProjectionSub_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>0,
                obj.Handles.Values.ProjectionSubsampling = TempInput;
            else
                obj.Handles.UIElements.ProjectionSig_Edit.String = num2str(obj.Handles.Values.ProjectionSig);
            end
        end

        function DecayCNMFE_Edit_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.DecayCNMFE_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>0 && rem(TempInput,1) == 0,
                obj.Parameters.CNMFE.max_tau = TempInput;
                obj.Parameters.CNMFE.Neuron.options.deconv_options.max_tau = TempInput;
            else
                obj.Handles.UIElements.DecayCNMFE_Edit.String = num2str(obj.Parameters.CNMFE.max_tau);
            end
        end

        function DetrendCNMFE_Edit_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.DetrendCNMFE_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>0 && rem(TempInput,1) == 0,
                obj.Parameters.CNMFE.Neuron.options.nk = TempInput;
            else
                obj.Handles.UIElements.DetrendCNMFE_Edit.String = num2str(obj.Parameters.Neuron.options.nk);
            end
        end


        function maxRAMCNMFE_Edit_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.maxRAMCNMFE_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>0 && rem(TempInput,1) == 0,
                obj.Parameters.CNMFE.maxRAM = TempInput;
            else
                obj.Handles.UIElements.maxRAMCNMFE_Edit.String = num2str(obj.Parameters.CNMFE.maxRAM);
            end
        end
        

        function ProjectionCellSize_Edit_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.ProjectionCellSize_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>0 && rem(TempInput,1) == 0,
                Center = obj.Parameters.CNMFE.Cell;
                obj.Parameters.CNMFE.gSiz = TempInput;
                Theta = 0:pi/50:2*pi;
                XCoor = obj.Parameters.CNMFE.gSiz/2 * cos(Theta) + Center(1);
                YCoor = obj.Parameters.CNMFE.gSiz/2 * sin(Theta) + Center(2);
                obj.Handles.UIElements.NeuronROI.XData  = XCoor;
                obj.Handles.UIElements.NeuronROI.YData  = YCoor;
                obj.Handles.UIElements.NeuronCNMFE_Edit.String = num2str(TempInput);
                obj.Parameters.CNMFE.gSiz = TempInput;
                obj.Parameters.CNMFE.Neuron.options.gSiz = TempInput;
                obj.Parameters.CNMFE.Neuron.options.ring_radius = ceil(0.5*1.05*TempInput);
                obj.Handles.UIElements.RingRadiusCNMFE_Edit.String = num2str(ceil(0.5*1.05*TempInput));
            else
                obj.Handles.UIElements.ProjectionCellSize_Edit.String = num2str(obj.Parameters.CNMFE.gSiz);
            end
        end

        function RingRadiusCNMFE_Edit_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.RingRadiusCNMFE_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>0 && rem(TempInput,1) == 0
                obj.Parameters.CNMFE.Neuron.options.ring_radius = TempInput;
            else
                obj.Handles.UIElements.RingRadiusCNMFE_Edit.String = num2str(obj.Parameters.CNMFE.Neuron.options.ring_radius);
            end
        end


        function CNMFE_ProjectionMoveNeuron_CB(obj)
            if ~obj.Dragging
                obj.Handles.Values.LastCoordinates = obj.Handles.UIElements.MovieDisplay(3).Plot.CurrentPoint(1,[1 2]);
                obj.Dragging = true;
                obj.Handles.MainFigure.WindowButtonMotionFcn = @(~,~)obj.MovingNeuronROI;
                obj.Handles.MainFigure.WindowButtonUpFcn = @(~,~)obj.CNMFE_ProjectionMoveNeuron_CB;
            else
                obj.Dragging = false;
                obj.Handles.MainFigure.WindowButtonMotionFcn = [];
                obj.Handles.MainFigure.WindowButtonUpFcn = [];
            end
        end

        function MovingNeuronROI(obj)
            DeltaCoordinates = obj.Handles.UIElements.MovieDisplay(3).Plot.CurrentPoint(1,[1 2]) - obj.Handles.Values.LastCoordinates;
            obj.Handles.UIElements.NeuronROI.XData = obj.Handles.UIElements.NeuronROI.XData + DeltaCoordinates(1);
            obj.Handles.UIElements.NeuronROI.YData = obj.Handles.UIElements.NeuronROI.YData + DeltaCoordinates(2);
            obj.Handles.Values.LastCoordinates = obj.Handles.UIElements.MovieDisplay(3).Plot.CurrentPoint(1,[1 2]);
            obj.Parameters.CNMFE.Cell = obj.Parameters.CNMFE.Cell + DeltaCoordinates;
            drawnow
        end


        function CNMFE_ProjectionDrawCrop_CB(obj)
            obj.DisableAll; drawnow
            obj.Handles.UIElements.Crop = drawrectangle(obj.Handles.UIElements.MovieDisplay(3).Plot,...
                'FaceAlpha',0,...
                'LineWidth',3,...
                'Color',[0.6 0.3 0],...
                'Deletable',false,...
                'FaceSelectable',false);
            obj.Parameters.CNMFE.Crop = obj.Handles.UIElements.Crop.Position;
            obj.EnableAll; drawnow
            addlistener(obj.Handles.UIElements.Crop(:),'MovingROI',@(~,~)obj.CNMFE_ProjectionEditCrop_CB);
            obj.Handles.UIElements.ProjectionCrop_Delete.Enable = 'on';
            obj.Handles.UIElements.ProjectionCrop_Draw.Enable = 'off';
        end

        function CNMFE_ProjectionEditCrop_CB(obj)
            obj.Parameters.CNMFE.Crop = obj.Handles.UIElements.Crop.Position;

        end

        function CNMFE_ProjectionDeleteCrop_CB(obj)
            delete(obj.Handles.UIElements.Crop)
            obj.Parameters.CNMFE.Crop = [];

            obj.Handles.UIElements.ProjectionCrop_Delete.Enable = 'off';
            obj.Handles.UIElements.ProjectionCrop_Draw.Enable = 'on';
            drawnow
        end

        function Process_MaxProjection(obj)
            obj.DisableAll;
            drawnow
            R = obj.Handles.Values.ProjectionSubsampling;
            gSig = obj.Handles.Values.ProjectionSig;
            % Create Kernel
            psf = fspecial('gaussian', ceil(gSig*4+1), gSig);
            ind_nonzero = (psf(:)>=max(psf(:,1)));
            psf = psf-mean(psf(ind_nonzero));
            psf(~ind_nonzero) = 0;
            % Retrieve frames
            MovieRead = h5read(obj.Parameters.MotionCorrection.MotionCorrectionFile,'/mov');
            if isinteger(MovieRead)
                MovieRead = single(MovieRead);
            end
            MovieRead = MovieRead(:,:,1:R:end);
            % Filter with the same kind of filter as CNMFE
            Filtered = zeros(size(MovieRead));
            parfor FM = 1 : size(Filtered,3)
                Filtered(:,:,FM) = imfilter(MovieRead(:,:,FM),psf,'replicate');
            end
            MP = max(Filtered,[],3);
            % Normalize for easier CLim adjustment
            MP = (MP - min(MP,[],'all')) / (max(MP,[],'all') - min(MP,[],'all'));

            obj.Handles.UIElements.MovieDisplay(3).MaxProjection.CData = MP;
            obj.Handles.Values.MaxProjection = MP;
            % Adjust limits to preserve ratio
            [D1,D2] = size(obj.Handles.Values.MaxProjection);
            obj.Handles.UIElements.MovieDisplay(3).Plot.Units = 'pixels';
            AbsolutePosition = obj.Handles.UIElements.MovieDisplay(3).Plot.Position;
            obj.Handles.UIElements.MovieDisplay(3).Plot.Units = 'normalized';
            if D2/D1 >= AbsolutePosition(3)/AbsolutePosition(4),
                % Y needs to be adjusted
                YDelta = D2*AbsolutePosition(4)/AbsolutePosition(3) - D1;
                obj.Handles.UIElements.MovieDisplay(3).Plot.YLim = [1-0.5*YDelta D1+0.5*YDelta];
                obj.Handles.UIElements.MovieDisplay(3).Plot.XLim = [1 D2];
            else
                % X needs to be adjusted
                XDelta = D1*AbsolutePosition(3)/AbsolutePosition(4) - D2;
                obj.Handles.UIElements.MovieDisplay(3).Plot.XLim = [1-0.5*XDelta D2+0.5*XDelta];
                obj.Handles.UIElements.MovieDisplay(3).Plot.YLim = [1 D1];
            end
            obj.EnableAll;
            drawnow
        end

        function Process_Correlation(obj,src,~)
            obj.DisableAll;drawnow
            CS = str2double(src.UserData);
            Movie = obj.Handles.Values.File{CS};
            Sig = obj.Handles.Values.Sig(CS);
            R = obj.Handles.Values.CorrSubsampling(CS);

            options.d1 = obj.Parameters.PreFiltering.d1;        % image height
            options.d2 = obj.Parameters.PreFiltering.d2;        % image width
            options.gSig = Sig;    % width of the gaussian kernel approximating one neuron
            options.gSiz = 1;    % average size of neurons
            options.center_psf = true;
            switch Movie
                case 'PreProcessed'
                    MovieRead = h5read(obj.Parameters.PreProcessing.PreProcessingFile,'/mov');
                case 'MotionCorrected'
                    MovieRead = h5read(obj.Parameters.MotionCorrection.MotionCorrectionFile,'/mov');
            end

            tic;[Cn, PNR] = correlation_image_endoscope(MovieRead(:,:,1:R:end), options);toc
            obj.Handles.Values.PNR{CS} = PNR;
            obj.Handles.Values.Cn{CS} = Cn;
            obj.Handles.UIElements.Cn(CS).Children.CData = Cn;
            obj.Handles.UIElements.Cn(CS).CLim = [obj.Handles.Values.min_corr(CS) 1];
            obj.Handles.UIElements.PNR(CS).Children.CData = PNR;
            if obj.Handles.Values.min_pnr(CS) < 0.98*max(PNR,[],'all')
                obj.Handles.UIElements.PNR(CS).CLim = [obj.Handles.Values.min_pnr(CS) 0.98*max(PNR,[],'all')];
            else
                obj.Handles.UIElements.PNR(CS).CLim = [0.98*max(PNR,[],'all')-0.01 0.98*max(PNR,[],'all')];
                obj.Handles.Values.min_pnr(CS) = 0.98*max(PNR,[],'all')-0.01;
                obj.Handles.UIElements.PNRSlider_Edit(CS).String = num2str(0.98*max(PNR,[],'all')-0.01,'%.2f');
                obj.Handles.UIElements.PNRSlider_Slider(CS).XData = (0.98*max(PNR,[],'all')-0.01)*[1 1];
            end
            Prod = Cn.*PNR;
            Prod(Prod<obj.Handles.Values.min_pnr(CS) | Prod<obj.Handles.Values.min_corr(CS)) = 0;
            obj.Handles.UIElements.Product(CS).Children.CData = Prod;
            obj.Handles.UIElements.Product(CS).CLim = [obj.Handles.Values.min_corr(CS)*obj.Handles.Values.min_pnr(CS) 0.98*max(PNR,[],'all')];

            obj.Handles.UIElements.PNRSlider_Axes(CS).XLim = [0 max(obj.Handles.Values.PNR{CS},[],'all')];
            obj.Handles.UIElements.PNRSlider_BaseLine(CS).XData(2) = max(obj.Handles.Values.PNR{CS},[],'all');
            % Adjust limits to preserve ratio
            [D1,D2] = size(Cn);
            obj.Handles.UIElements.Cn(CS).Units = 'pixels';
            AbsolutePosition = obj.Handles.UIElements.Cn(CS).Position;
            obj.Handles.UIElements.Cn(CS).Units = 'normalized';
            if D2/D1 >= AbsolutePosition(3)/AbsolutePosition(4)
                % Y needs to be adjusted
                YDelta = D2*AbsolutePosition(4)/AbsolutePosition(3) - D1;
                obj.Handles.UIElements.Cn(CS).YLim = [1-0.5*YDelta D1+0.5*YDelta];
                obj.Handles.UIElements.Cn(CS).XLim = [1 D2];
                obj.Handles.UIElements.PNR(CS).YLim = [1-0.5*YDelta D1+0.5*YDelta];
                obj.Handles.UIElements.PNR(CS).XLim = [1 D2];
                obj.Handles.UIElements.Product(CS).YLim = [1-0.5*YDelta D1+0.5*YDelta];
                obj.Handles.UIElements.Product(CS).XLim = [1 D2];
            else
                % X needs to be adjusted
                XDelta = D1*AbsolutePosition(3)/AbsolutePosition(4) - D2;
                obj.Handles.UIElements.Cn(CS).XLim = [1-0.5*XDelta D2+0.5*XDelta];
                obj.Handles.UIElements.Cn(CS).YLim = [1 D1];
                obj.Handles.UIElements.PNR(CS).XLim = [1-0.5*XDelta D2+0.5*XDelta];
                obj.Handles.UIElements.PNR(CS).YLim = [1 D1];
                obj.Handles.UIElements.Product(CS).XLim = [1-0.5*XDelta D2+0.5*XDelta];
                obj.Handles.UIElements.Product(CS).YLim = [1 D1];
            end
            obj.EnableAll;
            if strcmpi(Movie,'MotionCorrected'),
                obj.Handles.UIElements.Set_Button(str2double(src.UserData)).Enable = 'on';
            end
            drawnow
            clear MovieRead
        end


        
        function CNMFE_Preprocessed_CB(obj,src,~)
            CS = str2double(src.UserData);
            obj.Handles.Values.File{CS} = 'PreProcessed';
            obj.Handles.UIElements.File_Preprocessed_Button(CS).Value = 1;
            obj.Handles.UIElements.File_MotionCorrected_Button(CS).Value = 0;
        end

        function CNMFE_MotionCorrected_CB(obj,src,~)
            CS = str2double(src.UserData);
            obj.Handles.Values.File{CS} = 'MotionCorrected';
            obj.Handles.UIElements.File_Preprocessed_Button(CS).Value = 0;
            obj.Handles.UIElements.File_MotionCorrected_Button(CS).Value = 1;
        end

        function CNMFE_Sig_Edit_CB(obj,src,~)
            TempInput = str2double(src.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>0 && obj.Handles.Values.Sig(str2double(src.UserData)) ~= TempInput,
                obj.Handles.Values.Sig(str2double(src.UserData)) = TempInput;
                obj.Handles.UIElements.Set_Button(str2double(src.UserData)).Enable = 'off';
            else
                src.String = num2str(obj.Handles.Values.Sig(str2double(src.UserData)));
            end
        end

        function CNMFE_Sub_Edit_CB(obj,src,~)
            TempInput = str2double(src.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>0 && rem(TempInput,1) == 0 && obj.Handles.Values.CorrSubsampling(str2double(src.UserData)) ~= TempInput,
                obj.Handles.Values.CorrSubsampling(str2double(src.UserData)) = TempInput;
                obj.Handles.UIElements.Set_Button(str2double(src.UserData)).Enable = 'off';
            else
                src.String = num2str(obj.Handles.Values.CorrSubsampling(str2double(src.UserData)));
            end
        end


        function CNMFE_Set_CB(obj,src,~)
            CS = str2double(src.UserData);
            obj.Parameters.CNMFE.Neuron.options.min_corr = obj.Handles.Values.min_corr(CS);
            obj.Parameters.CNMFE.min_corr = obj.Handles.Values.min_corr(CS);
            obj.Handles.UIElements.mincorrCNMFE_Edit.String = num2str(obj.Handles.Values.min_corr(CS));
            obj.Parameters.CNMFE.Neuron.options.min_pnr = obj.Handles.Values.min_pnr(CS);
            obj.Parameters.CNMFE.min_pnr = obj.Handles.Values.min_pnr(CS);
            obj.Handles.UIElements.minPNRCNMFE_Edit.String = num2str(obj.Handles.Values.min_pnr(CS));
            obj.Handles.UIElements.ValidateCNMFE_Button.Enable = 'on';
            obj.Handles.Values.SetOnce = true;
            obj.Parameters.CNMFE.Neuron.options.gSig = obj.Handles.Values.Sig(CS);
            obj.Parameters.CNMFE.gSig = obj.Handles.Values.Sig(CS);
            obj.Handles.UIElements.SigmaCNMFE_Edit.String = obj.Handles.Values.Sig(CS);
        end

        function CnSlider_Edit_CB(obj,src,~)
            TempInput = str2double(src.String);
            CS = str2double(src.UserData);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>=0 && TempInput<1,
                obj.Handles.Values.min_corr(CS) = TempInput;
                obj.Handles.UIElements.CnSlider_Slider(CS).XData = TempInput * [1 1];
                if ~isempty(obj.Handles.Values.PNR{CS})
                    obj.Handles.UIElements.Cn(CS).CLim = [obj.Handles.Values.min_corr(CS) 1];
                    Prod = obj.Handles.Values.Cn{CS}.*obj.Handles.Values.PNR{CS};
                    Prod(obj.Handles.Values.Cn{CS}<obj.Handles.Values.min_corr(CS) | obj.Handles.Values.PNR{CS}<obj.Handles.Values.min_pnr(CS)) = 0;
                    obj.Handles.UIElements.Product(CS).Children.CData = Prod;
                    obj.Handles.UIElements.Product(CS).CLim = [obj.Handles.Values.min_corr(CS)*obj.Handles.Values.min_pnr(CS) 0.98*max(obj.Handles.Values.PNR{CS},[],'all')];
                    if obj.Handles.Values.Binarize(CS),
                        obj.Handles.UIElements.Cn(CS).CLim = obj.Handles.UIElements.Cn(CS).CLim(1)*[1 1+eps];
                        obj.Handles.UIElements.Product(CS).CLim = obj.Handles.UIElements.Product(CS).CLim(1)*[1 1+eps];
                    end
                end
                src.String = num2str(obj.Handles.Values.min_corr(CS),'%.2f');
            else
                src.String = num2str(obj.Handles.Values.min_corr(CS),'%.2f');
            end
        end

        function PNRSlider_Edit_CB(obj,src,~)
            TempInput = str2double(src.String);
            CS = str2double(src.UserData);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>=0 && TempInput<=15,
                if ~isempty(obj.Handles.Values.PNR{CS}),
                    if TempInput>max(obj.Handles.Values.PNR{CS},[],'all'),
                        src.String = num2str(obj.Handles.Values.min_pnr(CS),'%.2f');
                        return
                    end
                end
                obj.Handles.Values.min_pnr(CS) = TempInput;
                obj.Handles.UIElements.PNRSlider_Slider(CS).XData = TempInput * [1 1];
                if ~isempty(obj.Handles.Values.PNR{CS})
                    obj.Handles.UIElements.PNR(CS).CLim = [obj.Handles.Values.min_pnr(CS) 0.98*max(obj.Handles.Values.PNR{CS},[],'all')];
                    Prod = obj.Handles.Values.Cn{CS}.*obj.Handles.Values.PNR{CS};
                    Prod(obj.Handles.Values.Cn{CS}<obj.Handles.Values.min_corr(CS) | obj.Handles.Values.PNR{CS}<obj.Handles.Values.min_pnr(CS)) = 0;
                    obj.Handles.UIElements.Product(CS).Children.CData = Prod;
                    obj.Handles.UIElements.Product(CS).CLim = [obj.Handles.Values.min_corr(CS)*obj.Handles.Values.min_pnr(CS) 0.98*max(obj.Handles.Values.PNR{CS},[],'all')];
                    if obj.Handles.Values.Binarize(CS),
                        obj.Handles.UIElements.PNR(CS).CLim = obj.Handles.UIElements.PNR(CS).CLim(1)*[1 1+eps];
                        obj.Handles.UIElements.Product(CS).CLim = obj.Handles.UIElements.Product(CS).CLim(1)*[1 1+eps];
                    end
                end
                src.String = num2str(obj.Handles.Values.min_pnr(CS),'%.2f');
            else
                src.String = num2str(obj.Handles.Values.min_pnr(CS),'%.2f');
            end
        end


        function PNRSliderCB(obj,src,~)
            if ~obj.Dragging,
                obj.Handles.Values.Moved = str2double(src.UserData);
                obj.Dragging = true;
                obj.Handles.MainFigure.WindowButtonMotionFcn = @(~,~)obj.MovingSliderPNR;
                obj.Handles.MainFigure.WindowButtonUpFcn = @(~,~)obj.PNRSliderCB;
            else
                obj.Dragging = false;
                obj.Handles.MainFigure.WindowButtonMotionFcn = [];
                obj.Handles.MainFigure.WindowButtonUpFcn = [];
            end
        end

        function MovingSliderPNR(obj)
            CS = obj.Handles.Values.Moved;
            CurrentCursor = obj.Handles.UIElements.PNRSlider_Axes(CS).CurrentPoint;
            if CurrentCursor(1)>=obj.Handles.UIElements.PNRSlider_Axes(CS).XLim(1) &&  CurrentCursor(1)<=obj.Handles.UIElements.PNRSlider_Axes(CS).XLim(2),
                if ~isempty(obj.Handles.Values.PNR{CS}),
                    if CurrentCursor(1)>0.98*max(obj.Handles.Values.PNR{CS},[],'all'),
                        return
                    end
                end
                obj.Handles.Values.min_pnr(CS) = CurrentCursor(1);
                obj.Handles.UIElements.PNRSlider_Edit(CS).String = num2str(CurrentCursor(1),'%.2f');
                obj.Handles.UIElements.PNRSlider_Slider(CS).XData = CurrentCursor(1) * [1 1];
                if ~isempty(obj.Handles.Values.PNR{CS}),
                    obj.Handles.UIElements.PNR(CS).CLim = [obj.Handles.Values.min_pnr(CS) 0.98*max(obj.Handles.Values.PNR{CS},[],'all')];
                    Prod = obj.Handles.Values.Cn{CS}.*obj.Handles.Values.PNR{CS};
                    Prod(obj.Handles.Values.Cn{CS}<obj.Handles.Values.min_corr(CS) | obj.Handles.Values.PNR{CS}<obj.Handles.Values.min_pnr(CS)) = 0;
                    obj.Handles.UIElements.Product(CS).Children.CData = Prod;
                    obj.Handles.UIElements.Product(CS).CLim = [obj.Handles.Values.min_corr(CS)*obj.Handles.Values.min_pnr(CS) 0.98*max(obj.Handles.Values.PNR{CS},[],'all')];
                    if obj.Handles.Values.Binarize(CS),
                        obj.Handles.UIElements.PNR(CS).CLim = obj.Handles.UIElements.PNR(CS).CLim(1)*[1 1+eps];
                        obj.Handles.UIElements.Product(CS).CLim = obj.Handles.UIElements.Product(CS).CLim(1)*[1 1+eps];
                    end
                end
            end
        end

        function CnSliderCB(obj,src,~)
            if ~obj.Dragging,
                obj.Handles.Values.Moved = str2double(src.UserData);
                obj.Dragging = true;
                obj.Handles.MainFigure.WindowButtonMotionFcn = @(~,~)obj.MovingSliderCn;
                obj.Handles.MainFigure.WindowButtonUpFcn = @(~,~)obj.CnSliderCB;
            else
                obj.Dragging = false;
                obj.Handles.MainFigure.WindowButtonMotionFcn = [];
                obj.Handles.MainFigure.WindowButtonUpFcn = [];
            end
        end

        function MovingSliderCn(obj)
            CS = obj.Handles.Values.Moved;
            CurrentCursor = obj.Handles.UIElements.CnSlider_Axes(CS).CurrentPoint;
            if CurrentCursor(1)>=obj.Handles.UIElements.CnSlider_Axes(CS).XLim(1) &&  CurrentCursor(1)<=obj.Handles.UIElements.CnSlider_Axes(CS).XLim(2),
                obj.Handles.Values.min_corr(CS) = CurrentCursor(1);
                obj.Handles.UIElements.CnSlider_Edit(CS).String = num2str(CurrentCursor(1),'%.2f');
                obj.Handles.UIElements.CnSlider_Slider(CS).XData = CurrentCursor(1) * [1 1];
                if ~isempty(obj.Handles.Values.Cn{CS}),
                    obj.Handles.UIElements.Cn(CS).CLim = [obj.Handles.Values.min_corr(CS) 1];
                    Prod = obj.Handles.Values.Cn{CS}.*obj.Handles.Values.PNR{CS};
                    Prod(obj.Handles.Values.Cn{CS}<obj.Handles.Values.min_corr(CS) | obj.Handles.Values.PNR{CS}<obj.Handles.Values.min_pnr(CS)) = 0;
                    obj.Handles.UIElements.Product(CS).Children.CData = Prod;
                    obj.Handles.UIElements.Product(CS).CLim = [obj.Handles.Values.min_corr(CS)*obj.Handles.Values.min_pnr(CS) 0.98*max(obj.Handles.Values.PNR{CS},[],'all')];
                    if obj.Handles.Values.Binarize(CS),
                        obj.Handles.UIElements.Cn(CS).CLim = obj.Handles.UIElements.Cn(CS).CLim(1)*[1 1+eps];
                        obj.Handles.UIElements.Product(CS).CLim = obj.Handles.UIElements.Product(CS).CLim(1)*[1 1+eps];
                    end
                end
            end
        end

        function CLimSliderLowCB(obj)
            if ~obj.Dragging,
                obj.Dragging = true;
                obj.Handles.MainFigure.WindowButtonMotionFcn = @(~,~)obj.MovingCLimSliderLow;
                obj.Handles.MainFigure.WindowButtonUpFcn = @(~,~)obj.CLimSliderLowCB;
            else
                obj.Dragging = false;
                obj.Handles.MainFigure.WindowButtonMotionFcn = [];
                obj.Handles.MainFigure.WindowButtonUpFcn = [];
            end
        end

        function MovingCLimSliderLow(obj)
            CurrentCursor = obj.Handles.UIElements.CLimSlider_Axes.CurrentPoint;
            if CurrentCursor(1)>=obj.Handles.UIElements.CLimSlider_Axes.XLim(1) &&  CurrentCursor(1)<=obj.Handles.UIElements.CLimSlider_Axes.XLim(2) && CurrentCursor(1)<obj.Handles.UIElements.MovieDisplay(3).Plot.CLim(2),
                obj.Handles.Values.ProjectionCLim(1) = CurrentCursor(1);
                obj.Handles.UIElements.CLimSlider_SliderLow.XData = CurrentCursor(1) * [1 1];
                obj.Handles.UIElements.MovieDisplay(3).Plot.CLim = obj.Handles.Values.ProjectionCLim;
            end
        end

        function CLimSliderHighCB(obj)
            if ~obj.Dragging,
                obj.Dragging = true;
                obj.Handles.MainFigure.WindowButtonMotionFcn = @(~,~)obj.MovingCLimSliderHigh;
                obj.Handles.MainFigure.WindowButtonUpFcn = @(~,~)obj.CLimSliderHighCB;
            else
                obj.Dragging = false;
                obj.Handles.MainFigure.WindowButtonMotionFcn = [];
                obj.Handles.MainFigure.WindowButtonUpFcn = [];
            end
        end

        function MovingCLimSliderHigh(obj)
            CurrentCursor = obj.Handles.UIElements.CLimSlider_Axes.CurrentPoint;
            if CurrentCursor(1)>=obj.Handles.UIElements.CLimSlider_Axes.XLim(1) &&  CurrentCursor(1)<=obj.Handles.UIElements.CLimSlider_Axes.XLim(2) &&  CurrentCursor(1)>obj.Handles.UIElements.MovieDisplay(3).Plot.CLim(1),
                obj.Handles.Values.ProjectionCLim(2) = CurrentCursor(1);
                obj.Handles.UIElements.CLimSlider_SliderHigh.XData = CurrentCursor(1) * [1 1];
                obj.Handles.UIElements.MovieDisplay(3).Plot.CLim = obj.Handles.Values.ProjectionCLim;
            end
        end

        function Binarize_CB(obj,src,~)
            CS = str2double(src.UserData);
            obj.Handles.Values.Binarize(str2double(src.UserData)) = logical(src.Value);
            if ~isempty(obj.Handles.Values.Cn{CS})
                obj.Handles.UIElements.Cn(CS).CLim = [obj.Handles.Values.min_corr(CS) 1];
                if obj.Handles.Values.min_pnr(CS) < 0.98*max(obj.Handles.Values.PNR{CS},[],'all')
                    obj.Handles.UIElements.PNR(CS).CLim = [obj.Handles.Values.min_pnr(CS) 0.98*max(obj.Handles.Values.PNR{CS},[],'all')];
                else
                    obj.Handles.UIElements.PNR(CS).CLim = [0.98*max(obj.Handles.Values.PNR{CS},[],'all')-0.01 0.98*max(obj.Handles.Values.PNR{CS},[],'all')];
                    obj.Handles.Values.min_pnr(CS) = 0.98*max(obj.Handles.Values.PNR{CS},[],'all')-0.01;
                end
                obj.Handles.UIElements.Product(CS).CLim = [obj.Handles.Values.min_corr(CS)*obj.Handles.Values.min_pnr(CS) 0.98*max(obj.Handles.Values.PNR{CS},[],'all')];
                if obj.Handles.Values.Binarize(CS),
                    obj.Handles.UIElements.PNR(CS).CLim = obj.Handles.UIElements.PNR(CS).CLim(1)*[1 1+eps];
                    obj.Handles.UIElements.Cn(CS).CLim = obj.Handles.UIElements.Cn(CS).CLim(1)*[1 1+eps];
                    obj.Handles.UIElements.Product(CS).CLim = obj.Handles.UIElements.Product(CS).CLim(1)*[1 1+eps];
                end
            end
        end

        %% Postprocessing tab & functions
        function Status = PostProcessing_Layout(obj,varargin)
            if numel(varargin)>0
                Force = varargin{1};
            else
                Force = false;
            end
            Status = false;
            if (~strcmpi(obj.CurrentTab,'PostProcessing')) || Force
                if ~(isfield(obj.Parameters,'CNMFE')) || ~obj.Parameters.CNMFE.Processed
                    return
                else
                    Status = true;

                    % Wipe current layout
                    obj.WipeLayout;

                    obj.Handles.Values.CaMovieInit = false;
                    obj.Handles.Values.CaMovieFilteredInit = false;


                    % Maximum projection (to see the ROIs better)
                    if (~isfield(obj.Handles,'Values') || ~isfield(obj.Handles.Values,'MaxProjectionPP'))
                        % We need to initialize
                        obj.Handles.Values.PostProcessing.Filtering.Sig1 = 0.01/obj.Parameters.PreProcessing.Spatial_Downsampling;
                        obj.Handles.Values.PostProcessing.Filtering.Sig2 = 2/obj.Parameters.PreProcessing.Spatial_Downsampling;
                        obj.Handles.Values.PostProcessing.Filtering.WindowSize = 60/obj.Parameters.PreProcessing.Spatial_Downsampling;
                        obj.Handles.Values.PostProcessing.AutoCLim = false;
                        obj.Handles.Values.PostProcessing.Movie1Filtered = true;
                        obj.Handles.Values.PostProcessing.Movie2Mode = 'none'; % 'none', 'Thermal', 'RGB', 'PS3', 'ELP' -if available, otherwise nothing happens
                        obj.Handles.Values.ProjectionCLim = [0 1];
                        obj.Handles.Values.CaMovieInit = false;
                        obj.Handles.Values.CaMovieFilteredInit = false;
                        obj.UpdatePFKernelPP;

                        % Check which behaviour movies are available
                        obj.Handles.Values.PostProcessing.AvailableMovieModes = {'none'};
                        Folder = fileparts(obj.Parameters.PreProcessing.RawFile);
                        RGBFile = [Folder filesep obj.Parameters.PreProcessing.Basename '.avi'];
                        ThermalFile = [Folder filesep obj.Parameters.PreProcessing.Basename '.mj2'];
                        PS3File = [Folder filesep obj.Parameters.PreProcessing.Basename '_PS3.avi'];
                        ELPFile = [Folder filesep obj.Parameters.PreProcessing.Basename '_ELP.avi'];

                        if exist(RGBFile,'file') == 2
                            obj.Handles.Values.PostProcessing.AvailableMovieModes = [obj.Handles.Values.PostProcessing.AvailableMovieModes,'RGB'];
                        end
                        if exist(ThermalFile,'file') == 2
                            obj.Handles.Values.PostProcessing.AvailableMovieModes = [obj.Handles.Values.PostProcessing.AvailableMovieModes,'Thermal'];
                        end
                        if exist(PS3File,'file') == 2
                            obj.Handles.Values.PostProcessing.AvailableMovieModes = [obj.Handles.Values.PostProcessing.AvailableMovieModes,'PS3'];
                        end
                        if exist(ELPFile,'file') == 2
                            obj.Handles.Values.PostProcessing.AvailableMovieModes = [obj.Handles.Values.PostProcessing.AvailableMovieModes,'ELP'];
                        end

                        if ~isfield(obj.Parameters.CNMFE,'MaxProjection')
                            gSig = obj.Parameters.CNMFE.gSig;
                            % Create Kernel
                            psf = fspecial('gaussian', ceil(gSig*4+1), gSig);
                            ind_nonzero = (psf(:)>=max(psf(:,1)));
                            psf = psf-mean(psf(ind_nonzero));
                            psf(~ind_nonzero) = 0;
                            % Retrieve frames
                            MovieRead = h5read(obj.Parameters.MotionCorrection.MotionCorrectionFile,'/mov');
                            if isinteger(MovieRead)
                                MovieRead = single(MovieRead);
                            end
                            MovieRead = MovieRead(:,:,1:end);
                            % Filter with the same kind of filter as CNMFE
                            Filtered = zeros(size(MovieRead));
                            parfor FM = 1 : size(Filtered,3),
                                Filtered(:,:,FM) = imfilter(MovieRead(:,:,FM),psf,'replicate');
                            end
                            MP = max(Filtered,[],3);
                            % Normalize for easier CLim adjustment
                            MP = (MP - min(MP,[],'all')) / (max(MP,[],'all') - min(MP,[],'all'));
                            obj.Parameters.CNMFE.MaxProjection = MP;
                        end

                        % Set movies structure
                        obj.CleanMovies;

                        % Main display and uielements
                        obj.SetMovieDisplay1;

                        if strcmpi(obj.Handles.Values.PostProcessing.Movie2Mode,'none')
                            obj.CurrentPlayers = struct(...
                                'Axes',{obj.Handles.UIElements.MovieDisplay(1).Plot},...
                                'Data','MotionCorrected',...
                                'Filter',{obj.Handles.Values.PostProcessing.Movie1Filtered},...
                                'Kernel','PostProcessing',...
                                'Readers',[],...
                                'Plot',[],...
                                'AbsolutePosition',{obj.Handles.UIElements.MovieDisplay(1).AbsolutePosition});
                            obj.Handles.UIElements.Traces = axes(...
                                'Units','normalized',...
                                'Position',[0.61 0.1 0.3500 0.8],...
                                'Color','w','Parent',obj.Handles.MainFigure);hold on
                        else
                            obj.SetMovieDisplay2;
                            switch obj.Handles.Values.PostProcessing.Movie2Mode
                                case 'RGB'
                                    MFile = [Folder filesep obj.Parameters.PreProcessing.Basename '.avi'];
                                case 'Thermal'
                                    MFile = [Folder filesep obj.Parameters.PreProcessing.Basename '.mj2'];
                                case 'PS3'
                                    MFile = [Folder filesep obj.Parameters.PreProcessing.Basename '_PS3.avi'];
                                case 'ELP'
                                    MFile = [Folder filesep obj.Parameters.PreProcessing.Basename '_ELP.avi'];
                            end
                            obj.CurrentPlayers = struct(...
                                'Axes',{obj.Handles.UIElements.MovieDisplay(1).Plot,obj.Handles.UIElements.MovieDisplay(2).Plot},...
                                'Data',{'MotionCorrected',MFile},...
                                'Filter',{obj.Handles.Values.PostProcessing.Movie1Filtered,false},...
                                'Kernel','PostProcessing',...
                                'Readers',[],...
                                'Plot',[],...
                                'AbsolutePosition',{obj.Handles.UIElements.MovieDisplay(1).AbsolutePosition,obj.Handles.UIElements.MovieDisplay(2).AbsolutePosition});

                            obj.Handles.UIElements.Traces = axes(...
                                'Units','normalized',...
                                'Position',[0.6050 0.1 0.3600 0.4],...
                                'Color','w','Parent',obj.Handles.MainFigure);hold on
                        end
                        obj.Handles.UIElements.Traces.YTick = [];
                        drawnow
                        obj.Handles.UIElements.Traces.LineWidth = 2.5;
                        obj.Handles.UIElements.Traces.XColor = [0.8 0.8 0.8];
                        obj.Handles.UIElements.Traces.YColor = [0.8 0.8 0.8];
                        obj.Handles.UIElements.Traces.TickDir = 'out';
                        obj.Handles.UIElements.Traces.Box = 'off';
                        obj.Handles.UIElements.Traces.XLim = obj.Handles.UIElements.MovieDisplay(1).Plot.XLim;

                        obj.Handles.Values.PostProcessing.Filtering.Sig1 = 0.01/obj.Parameters.PreProcessing.Spatial_Downsampling;
                        obj.Handles.Values.PostProcessing.Filtering.Sig2 = 2/obj.Parameters.PreProcessing.Spatial_Downsampling;
                        obj.Handles.Values.PostProcessing.Filtering.WindowSize = 60/obj.Parameters.PreProcessing.Spatial_Downsampling;
                        obj.Handles.Values.PostProcessing.AutoCLim = false;
                        obj.Handles.Values.PostProcessing.DisplayNumbers = true;

                        % Load movie
                        obj.LoadMovies;
                        obj.Handles.Values.CaMovieInit = false;

                        if obj.Handles.Values.PostProcessing.AutoCLim
                            obj.Handles.UIElements.MovieDisplay(1).Plot.CLimMode = 'auto';
                            if obj.Handles.Values.PostProcessing.Movie1Filtered
                                obj.Handles.Values.CaMovieFilteredCLim = [-10 0.75*obj.Handles.UIElements.MovieDisplay(1).Plot.CLim(2)];
                                obj.Handles.UIElements.MovieDisplay(1).Plot.CLim =  obj.Handles.Values.CaMovieFilteredCLim;
                                obj.Handles.Values.CaMovieFilteredInit = true;
                            else
                                obj.Handles.Values.CaMovieCLim = obj.Handles.UIElements.MovieDisplay(1).Plot.CLim;
                                obj.Handles.Values.CaMovieInit = true;
                            end
                        else
                            if obj.Handles.Values.PostProcessing.Movie1Filtered
                                if obj.Handles.Values.CaMovieFilteredInit
                                    obj.Handles.UIElements.MovieDisplay(1).Plot.CLim = obj.Handles.Values.CaMovieFilteredCLim;
                                else
                                    obj.Handles.UIElements.MovieDisplay(1).Plot.CLimMode = 'auto';
                                    obj.Handles.Values.CaMovieFilteredCLim = [-10 0.75*obj.Handles.UIElements.MovieDisplay(1).Plot.CLim(2)];
                                    obj.Handles.UIElements.MovieDisplay(1).Plot.CLim =  obj.Handles.Values.CaMovieFilteredCLim;
                                    obj.Handles.Values.CaMovieFilteredInit = true;
                                    obj.Handles.UIElements.MovieDisplay(1).Plot.CLimMode = 'manual';
                                end
                            else
                                if obj.Handles.Values.CaMovieInit
                                    obj.Handles.UIElements.MovieDisplay(1).Plot.CLim = obj.Handles.Values.CaMovieCLim;
                                else
                                    obj.Handles.UIElements.MovieDisplay(1).Plot.CLimMode = 'auto';
                                    obj.Handles.Values.CaMovieCLim = obj.Handles.UIElements.MovieDisplay(1).Plot.CLim;
                                    obj.Handles.Values.CaMovieInit = true;
                                    obj.Handles.UIElements.MovieDisplay(1).Plot.CLimMode = 'manual';
                                end
                            end
                        end


                        % Max projection axes
                        obj.Handles.UIElements.MovieDisplay(3).Plot = axes(...
                            'Units','normalized',...
                            'Position',[0.195 0.04  0.2048 0.2560],...
                            'Color','k','XColor','none','YColor','none',...
                            'Colormap',bone,'UserData','3','Parent',obj.Handles.MainFigure);hold on;
                        obj.Handles.UIElements.MovieDisplay(3).Box = axes(...
                            'Units','normalized',...
                            'Position',[0.195 0.04  0.2048 0.2560],...
                            'Color','none','XColor','none','YColor','none',...
                            'Colormap',bone,'UserData','3','Parent',obj.Handles.MainFigure);hold on
                        obj.Handles.UIElements.MovieDisplay(3).Plot.YDir = 'normal';
                        obj.Handles.UIElements.MovieDisplay(3).Box.YDir = 'normal';
                        obj.Handles.UIElements.MovieDisplay(3).Plot.Toolbar.Visible = 'off';
                        obj.Handles.UIElements.MovieDisplay(3).Plot.Interactions = [];
                        obj.Handles.UIElements.MovieDisplay(3).Box.Toolbar.Visible = 'on';
                        obj.Handles.UIElements.MovieDisplay(3).Tb = axtoolbar(obj.Handles.UIElements.MovieDisplay(3).Box, {'zoomin', 'zoomout', 'pan'});
                        Button = axtoolbarbtn(obj.Handles.UIElements.MovieDisplay(3).Tb, 'push');
                        Button.Icon = 'restoreview';
                        Button.ButtonPushedFcn = @obj.Restoreview;

                        imagesc(obj.Parameters.CNMFE.MaxProjection,'Parent',obj.Handles.UIElements.MovieDisplay(3).Plot)
                        obj.Handles.UIElements.MovieDisplay(3).Plot.CLim = [0 1];

                        % Adjust limits to preserve ratio
                        [D1,D2] = size(obj.Parameters.CNMFE.MaxProjection);
                        obj.Handles.UIElements.MovieDisplay(3).Plot.Units = 'pixels';
                        AbsolutePosition = obj.Handles.UIElements.MovieDisplay(3).Plot.Position;
                        obj.Handles.UIElements.MovieDisplay(3).AbsolutePosition = AbsolutePosition;
                        obj.Handles.UIElements.MovieDisplay(3).Plot.Units = 'normalized';
                        if D2/D1 >= AbsolutePosition(3)/AbsolutePosition(4)
                            % Y needs to be adjusted
                            YDelta = D2*AbsolutePosition(4)/AbsolutePosition(3) - D1;
                            obj.Handles.UIElements.MovieDisplay(3).Plot.YLim = [1-0.5*YDelta D1+0.5*YDelta];
                            obj.Handles.UIElements.MovieDisplay(3).Plot.XLim = [1 D2];
                        else
                            % X needs to be adjusted
                            XDelta = D1*AbsolutePosition(3)/AbsolutePosition(4) - D2;
                            obj.Handles.UIElements.MovieDisplay(3).Plot.XLim = [1-0.5*XDelta D2+0.5*XDelta];
                            obj.Handles.UIElements.MovieDisplay(3).Plot.YLim = [1 D1];
                        end
                        obj.Handles.UIElements.MovieDisplay(3).SourceXLim = obj.Handles.UIElements.MovieDisplay(3).Plot.XLim;
                        obj.Handles.UIElements.MovieDisplay(3).SourceYLim = obj.Handles.UIElements.MovieDisplay(3).Plot.YLim;
                        obj.Handles.UIElements.MovieDisplay(3).Box.XLim =  obj.Handles.UIElements.MovieDisplay(3).Plot.XLim;
                        obj.Handles.UIElements.MovieDisplay(3).Box.YLim =  obj.Handles.UIElements.MovieDisplay(3).Plot.YLim;


                        obj.Handles.UIElements.MovieDisplay(3).ParametersBox = axes(...
                            'Units','normalized','Parent',obj.Handles.MainFigure,...
                            'Position',[0.405 0.04 0.15 0.2560],...
                            'Color',[0.2 0.2 0.2],'XColor','k','YColor','k','Box','on','LineWidth',1.5);
                        obj.Handles.UIElements.MovieDisplay(3).ParametersBox.XTick = [];
                        obj.Handles.UIElements.MovieDisplay(3).ParametersBox.YTick = [];
                        obj.Handles.UIElements.MaxProjection_Legend = uicontrol('Style','text',...
                            'Units','Normalized','Position',[obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(1) obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(2)+obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(4)-0.02 obj.Handles.UIElements.MovieDisplay(3).ParametersBox.Position(3)+0.0005 0.02],...
                            'String','ROIs operations',...
                            'FontSize',13,'FontName','Arial','FontWeight','b',...
                            'BackgroundColor','k','ForegroundColor',[0.6 0.6 0.6],...
                            'HorizontalAlignment','center','Parent',obj.Handles.MainFigure);


                        obj.Handles.UIElements.MergeROIs = uicontrol('Style','pushbutton',...
                            'Units','Normalized','Position',[0.4125 0.23 0.065 0.035],...
                            'String', 'Merge selected',...
                            'Tooltip', ['Merge the different ROIs selected, and run a fast update (spatial/temporal).' newline 'It is advised to "Run iteration" at the end of the manual operations and before exporting.'] ,...
                            'FontSize',12,'FontName','Arial','FontWeight','b',...
                            'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                            'HorizontalAlignment','center',...
                            'Callback',@(~,~)obj.PostProc_MergeROIs_CB,'Parent',obj.Handles.MainFigure);
                        %                      obj.Handles.UIElements.CancelMergeROIs = uicontrol('Style','pushbutton',...
                        %                          'Units','Normalized','Position',[0.485 0.23 0.065 0.035],...
                        %                          'String', 'Cancel last',...
                        %                          'FontSize',12,'FontName','Arial','FontWeight','b',...
                        %                          'BackgroundColor',[0.3 0.3 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                        %                          'HorizontalAlignment','center',...
                        %                          'Callback',@(~,~)obj.PostProc_CancelMergeROIs_CB);

                        obj.Handles.UIElements.DeleteROIs = uicontrol('Style','pushbutton',...
                            'Units','Normalized','Position',[0.4125 0.19 0.065 0.035],...
                            'String', 'Delete selected',...,...
                            'Tooltip', 'Delete the selected ROIs.',...
                            'FontSize',12,'FontName','Arial','FontWeight','b',...
                            'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                            'HorizontalAlignment','center',...
                            'Callback',@(~,~)obj.PostProc_DeleteROIs_CB,'Parent',obj.Handles.MainFigure);
                        %                      obj.Handles.UIElements.CancelDeleteROIs = uicontrol('Style','pushbutton',...
                        %                          'Units','Normalized','Position',[0.485 0.19 0.065 0.035],...
                        %                          'String', 'Cancel last',...
                        %                          'FontSize',12,'FontName','Arial','FontWeight','b',...
                        %                          'BackgroundColor',[0.3 0.2 0.3],'ForegroundColor',[0.6 0.6 0.6],...
                        %                          'HorizontalAlignment','center',...
                        %                          'Callback',@(~,~)obj.PostProc_CancelDeleteROIs_CB);

                        obj.Handles.UIElements.RunIteration = uicontrol('Style','pushbutton',...
                            'Units','Normalized','Position',[0.4125 0.15 0.065 0.035],...
                            'String', 'Run iteration',...
                            'Tooltip', ['Run a few loops of spatial/temporal updating, and background estimation.' newline 'It is "allowed" to do it several times, and it can improve the quality of the ROIs/traces. But usually, past a certain point, it can also deteriorate them.'],...
                            'FontSize',12,'FontName','Arial','FontWeight','b',...
                            'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                            'HorizontalAlignment','center',...
                            'Callback',@(~,~)obj.PostProc_RunIteration_CB,'Parent',obj.Handles.MainFigure);

                        %
                        %                      obj.Handles.UIElements.CancelRunIteration = uicontrol('Style','pushbutton',...
                        %                          'Units','Normalized','Position',[0.485 0.15 0.065 0.035],...
                        %                          'String', 'Cancel last',...
                        %                          'FontSize',12,'FontName','Arial','FontWeight','b',...
                        %                          'BackgroundColor',[0.2 0.3 0.3],'ForegroundColor',[0.6 0.6 0.6],...
                        %                          'HorizontalAlignment','center',...
                        %                          'Callback',@(~,~)obj.PostProc_CancelRunIteration_CB);

                        obj.Handles.UIElements.CancelLast = uicontrol('Style','pushbutton',...
                            'Units','Normalized','Position',[0.485 0.19 0.065 0.035],...
                            'String', ['Cancel last' newline 'action'],...,...
                            'Tooltip','Cancel the last operation (merging, deleting, or running iteration).',...
                            'FontSize',12,'FontName','Arial','FontWeight','b',...
                            'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                            'HorizontalAlignment','center',...
                            'Callback',@(~,~)obj.PostProc_CancelLast_CB,'Parent',obj.Handles.MainFigure);

                        obj.Handles.UIElements.Restore = uicontrol('Style','pushbutton',...
                            'Units','Normalized','Position',[0.485 0.15 0.065 0.035],...
                            'String', ['Restore'],...,...
                            'Tooltip', 'Restore everything as it was just when CNMFE was performed. All manual postprocessing will be lost.',...
                            'FontSize',12,'FontName','Arial','FontWeight','b',...
                            'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                            'HorizontalAlignment','center',...
                            'Callback',@(~,~)obj.PostProc_Restore,'Parent',obj.Handles.MainFigure);

                        obj.Handles.UIElements.CancelSelection = uicontrol('Style','pushbutton',...
                            'Units','Normalized','Position',[0.485 0.23 0.065 0.035],...
                            'String', 'Unselect all',...,...
                            'FontSize',12,'FontName','Arial','FontWeight','b',...
                            'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                            'HorizontalAlignment','center',...
                            'Callback',@(~,~)obj.PostProc_CancelSelection_CB,'Parent',obj.Handles.MainFigure);

                        obj.Handles.UIElements.Save = uicontrol('Style','pushbutton',...
                            'Units','Normalized','Position',[0.4125 0.105 0.065 0.035],...
                            'String', 'Save',...,...
                            'Tooltip', 'Save the changes made manually, as a safepoint or to continue later.',...
                            'FontSize',13,'FontName','Arial','FontWeight','b',...
                            'BackgroundColor',[0.3 0.3 0.45],'ForegroundColor',[0.6 0.6 0.6],...
                            'HorizontalAlignment','center',...
                            'Callback',@(~,~)obj.PostProc_Save_CB,'Parent',obj.Handles.MainFigure);


                        obj.Handles.UIElements.Extract = uicontrol('Style','pushbutton',...
                            'Units','Normalized','Position',[0.485 0.105 0.065 0.035],...
                            'String', 'Extract',...,...
                            'Tooltip', 'Export the results to a file that will be used for analyses, and tag the session to move all the files generated during the processing to the longterm storage.',...
                            'FontSize',13,'FontName','Arial','FontWeight','b',...
                            'BackgroundColor',[0.3 0.3 0.45],'ForegroundColor',[0.6 0.6 0.6],...
                            'HorizontalAlignment','center',...
                            'Callback',@(~,~)obj.PostProc_Extract_CB,'Parent',obj.Handles.MainFigure);


                        obj.Handles.UIElements.ProjectionCLim_Legend = uicontrol('Style','text',...
                            'Units','Normalized','Position',[0.4125 0.05 0.055/2 0.03],...
                            'String','CLim',...
                            'FontSize',13,'FontName','Arial','FontWeight','b',...
                            'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                            'HorizontalAlignment','right','Parent',obj.Handles.MainFigure);
                        obj.Handles.UIElements.CLimSlider_Axes = axes('Units','normalized',...
                            'Position',[0.4125+0.005+0.055/2 0.05 0.1025 0.03],...
                            'Color',[0.2,0.2,0.2],'LineWidth',1,'Parent',obj.Handles.MainFigure); hold on;
                        obj.Handles.UIElements.CLimSlider_BaseLine = plot([0 1],[0 0],...
                            'Parent',obj.Handles.UIElements.CLimSlider_Axes,...
                            'LineWidth',1.5,'Color','k','Parent',obj.Handles.UIElements.CLimSlider_Axes);
                        obj.Handles.UIElements.CLimSlider_SliderLow = plot(obj.Handles.Values.ProjectionCLim(1)*[1 1],[-1 1],...
                            'Parent',obj.Handles.UIElements.CLimSlider_Axes,...
                            'LineWidth',3,'Color','k',...
                            'ButtonDownFcn',{@(~,~)obj.CLimSliderLowCB},'Parent',obj.Handles.UIElements.CLimSlider_Axes);
                        obj.Handles.UIElements.CLimSlider_SliderHigh = plot(obj.Handles.Values.ProjectionCLim(2)*[1 1],[-1 1],...
                            'Parent',obj.Handles.UIElements.CLimSlider_Axes,...
                            'LineWidth',3,'Color','k',...
                            'ButtonDownFcn',{@(~,~)obj.CLimSliderHighCB},'Parent',obj.Handles.UIElements.CLimSlider_Axes);

                        obj.Handles.UIElements.CLimSlider_Axes.XAxis.Visible = 'off';
                        obj.Handles.UIElements.CLimSlider_Axes.YAxis.Visible = 'off';
                        obj.Handles.UIElements.CLimSlider_Axes.XLim = [0 1];
                        obj.Handles.UIElements.CLimSlider_Axes.YLim = [-1.5 1.5];



                        % Filtering
                        obj.Handles.UIElements.Filtering = uicontrol('Style','checkbox',...
                            'Units','Normalized','Position',[0.1950 0.9075 0.05 0.025],...
                            'String', 'Filtered',...
                            'Value', obj.Handles.Values.PostProcessing.Movie1Filtered,...
                            'FontSize',12,'FontName','Arial','FontWeight','b',...
                            'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                            'HorizontalAlignment','center',...
                            'Callback',@(~,~)obj.PostProc_Filtering_CB,'Parent',obj.Handles.MainFigure);

                        obj.Handles.UIElements.Sig1PP_Legend = uicontrol('Style','text',...
                            'Units','Normalized','Position',[0.24 0.903 0.02 0.025],...
                            'String', 'Sig1',...
                            'FontSize',12,'FontName','Arial','FontWeight','b',...
                            'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                            'HorizontalAlignment','left','Parent',obj.Handles.MainFigure);
                        obj.Handles.UIElements.Sig1PP_Edit = uicontrol('Style','edit',...
                            'Units','Normalized','Position',[0.26 0.9075 0.035 0.025],...
                            'String', num2str(obj.Handles.Values.PostProcessing.Filtering.Sig1),...
                            'FontSize',12,'FontName','Arial','FontWeight','b',...
                            'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                            'HorizontalAlignment','center',...
                            'Callback',@(~,~)obj.Sig1PP_CB,'Parent',obj.Handles.MainFigure);


                        obj.Handles.UIElements.Sig2PP_Legend = uicontrol('Style','text',...
                            'Units','Normalized','Position',[0.3 0.903 0.02 0.025],...
                            'String', 'Sig2',...
                            'FontSize',12,'FontName','Arial','FontWeight','b',...
                            'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                            'HorizontalAlignment','left','Parent',obj.Handles.MainFigure);
                        obj.Handles.UIElements.Sig2PP_Edit = uicontrol('Style','edit',...
                            'Units','Normalized','Position',[0.32 0.9075 0.035 0.025],...
                            'String', num2str(obj.Handles.Values.PostProcessing.Filtering.Sig2),...
                            'FontSize',12,'FontName','Arial','FontWeight','b',...
                            'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                            'HorizontalAlignment','center',...
                            'Callback',@(~,~)obj.Sig2PP_CB,'Parent',obj.Handles.MainFigure);

                        obj.Handles.UIElements.WindowSizePP_Legend = uicontrol('Style','text',...
                            'Units','Normalized','Position',[0.3625 0.903 0.02 0.025],...
                            'String', 'Win',...
                            'FontSize',12,'FontName','Arial','FontWeight','b',...
                            'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                            'HorizontalAlignment','left','Parent',obj.Handles.MainFigure);
                        obj.Handles.UIElements.WindowPP_Edit = uicontrol('Style','edit',...
                            'Units','Normalized','Position',[0.38 0.9075 0.035 0.025],...
                            'String', num2str(obj.Handles.Values.PostProcessing.Filtering.WindowSize),...
                            'FontSize',12,'FontName','Arial','FontWeight','b',...
                            'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                            'HorizontalAlignment','center',...
                            'Callback',@(~,~)obj.WindowPP_CB,'Parent',obj.Handles.MainFigure);

                        if ~obj.Handles.UIElements.Filtering.Value
                            obj.Handles.UIElements.Sig1PP_Edit.Enable = 'off';
                            obj.Handles.UIElements.Sig2PP_Edit.Enable = 'off';
                            obj.Handles.UIElements.WindowSizePP_Edit.Enable = 'off';
                        end
                        obj.Handles.UIElements.CLimPP_Legend = uicontrol('Style','text',...
                            'Units','Normalized','Position',[0.425 0.903 0.025 0.025],...
                            'String', 'CLim',...
                            'FontSize',12,'FontName','Arial','FontWeight','b',...
                            'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                            'HorizontalAlignment','left','Parent',obj.Handles.MainFigure);
                        delete(obj.Handles.UIElements.BasenameLegend)

                        obj.Handles.UIElements.AutoCLim = uicontrol('Style','checkbox',...
                            'Units','Normalized','Position',[0.45 0.9075 0.03 0.025],...
                            'String', 'Auto',...
                            'Value', obj.Handles.Values.PostProcessing.AutoCLim,...
                            'FontSize',12,'FontName','Arial','FontWeight','b',...
                            'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                            'HorizontalAlignment','center',...
                            'Callback',@(~,~)obj.PostProc_AutoCLim_CB,'Parent',obj.Handles.MainFigure);

                        if obj.Handles.Values.PostProcessing.Movie1Filtered
                            PPCLim = obj.Handles.Values.CaMovieFilteredCLim;
                        else
                            PPCLim = obj.Handles.Values.CaMovieCLim;
                        end

                        obj.Handles.UIElements.ThresholdLowPostProc_Edit = uicontrol('Style','edit',...
                            'Units','Normalized','Position',[0.485 0.9075 0.035 0.025],...
                            'String', num2str(PPCLim(1)),...
                            'FontSize',12,'FontName','Arial','FontWeight','b',...
                            'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                            'HorizontalAlignment','center',...
                            'Callback',@(~,~)obj.ThresholdLowPostProc_CB,'Parent',obj.Handles.MainFigure);

                        obj.Handles.UIElements.ThresholdHighPostProc_Edit = uicontrol('Style','edit',...
                            'Units','Normalized','Position',[0.525 0.9075 0.035 0.025],...
                            'String', num2str(PPCLim(2)),...
                            'FontSize',12,'FontName','Arial','FontWeight','b',...
                            'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                            'HorizontalAlignment','center',...
                            'Callback',@(~,~)obj.ThresholdHighPostProc_CB,'Parent',obj.Handles.MainFigure);
                        obj.PostProc_AutoCLim_CB;

                        % Second movie choice
                        IndxMovieChoice = find(strcmpi(obj.Handles.Values.PostProcessing.AvailableMovieModes,obj.Handles.Values.PostProcessing.Movie2Mode));
                        obj.Handles.UIElements.MovieLegend2 = uicontrol('Style','popupmenu',...
                            'Units','Normalized','Position',[0.6050 0.903 0.075 0.025],...
                            'String',obj.Handles.Values.PostProcessing.AvailableMovieModes,...
                            'Value',IndxMovieChoice,...
                            'FontSize',12,'FontName','Arial','FontWeight','b',...
                            'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                            'HorizontalAlignment','left','Parent',obj.Handles.MainFigure);


                        obj.Handles.Values.Contours = obj.Parameters.CNMFE.Neuron.get_contours(0.8);
                        Traces = obj.Parameters.CNMFE.Neuron.C_raw;
                        Mins = repmat(min(Traces,[],2),1,length(Traces(1,:)));
                        Maxs = repmat(max(Traces,[],2),1,length(Traces(1,:)));
                        obj.Handles.Values.Traces = smoothdata((Traces - Mins)./(Maxs - Mins),2,'gaussian',5);
                        obj.Handles.TracesPlot = [];
                        hold on
                        % Plot contours
                        if isempty(obj.Parameters.CNMFE.Crop)
                            obj.Parameters.CNMFE.Crop = [0 0];
                        end
                        obj.Handles.Contours_HandlesMP = arrayfun(@(x) fill(obj.Handles.Values.Contours{x}(1,:),obj.Handles.Values.Contours{x}(2,:),[0.6 0.6 0.6],'FaceColor','none','EdgeColor',[0.6 0.6 0.6],'LineWidth',1,'ButtonDownFcn',{@(src,evt)obj.ContourClick(src,evt)},'UserData',x,'Parent',obj.Handles.UIElements.MovieDisplay(3).Box),1:numel(obj.Handles.Values.Contours));
                        obj.Handles.Contours_HandlesMovie = arrayfun(@(x) fill(obj.Parameters.CNMFE.Crop(1)+obj.Handles.Values.Contours{x}(1,:),obj.Parameters.CNMFE.Crop(2)+obj.Handles.Values.Contours{x}(2,:),[0.6 0.6 0.6],'FaceColor','none','EdgeColor',[0.6 0.6 0.6],'LineWidth',1,'ButtonDownFcn',{@(src,evt)obj.ContourClick(src,evt)},'UserData',x,'Parent',obj.Handles.UIElements.MovieDisplay(1).Box),1:numel(obj.Handles.Values.Contours));
                        % Retrieve default plot colors to use when clicking contours
                        DefC = DefColors;
                        obj.Handles.Values.ColorsList = repmat(DefC,[1e4 1]);
                        obj.Handles.Values.SelectedROIs = false(1e4,1);
                        obj.Handles.Values.MergeIndex = false(1e4,1);

                        if obj.Handles.Values.PostProcessing.DisplayNumbers
                            DisplayOpt = 'on';
                        else
                            DisplayOpt = 'off';
                        end

                        % Additional UI elements
                        obj.Handles.UIElements.DisplayNumbers = uicontrol('Style','checkbox',...
                            'Units','Normalized','Position',[ ...
                            obj.Handles.UIElements.MovieDisplay(1).Box.Position(1)+obj.Handles.UIElements.MovieDisplay(1).Box.Position(3)+0.0025...
                            obj.Handles.UIElements.MovieDisplay(1).Box.Position(2)+0.8*obj.Handles.UIElements.MovieDisplay(1).Box.Position(4)...
                            0.05 0.025],...
                            'String', 'Numbers',...
                            'Value', obj.Handles.Values.PostProcessing.DisplayNumbers,...
                            'FontSize',12,'FontName','Arial','FontWeight','b',...
                            'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                            'HorizontalAlignment','center',...
                            'Callback',@(~,~)obj.PostProc_DisplayNumbers_CB,'Parent',obj.Handles.MainFigure);


                        % Check whether we have a longitudinal alignment
                        CatFile = [obj.Parameters.CNMFE.CNMFEFolder filesep obj.Parameters.PreProcessing.Basename '_CellCategories.mat'];
                        if exist(CatFile,'file')==2
                            CatTemp = load(CatFile);
                            obj.Parameters.CNMFE.Categories = CatTemp.Categories;
                            obj.Handles.Values.PostProcessing.DisplayGoodROIs = true;
                            obj.Handles.Values.PostProcessing.DisplayBadROIs = true;

                            obj.Handles.UIElements.DisplayGoodROIs = uicontrol('Style','checkbox',...
                                'Units','Normalized','Position',[ ...
                                obj.Handles.UIElements.MovieDisplay(1).Box.Position(1)+obj.Handles.UIElements.MovieDisplay(1).Box.Position(3)+0.0025...
                                obj.Handles.UIElements.MovieDisplay(1).Box.Position(2)+0.7*obj.Handles.UIElements.MovieDisplay(1).Box.Position(4)...
                                0.05 0.025],...
                                'String', 'Good',...
                                'Value', obj.Handles.Values.PostProcessing.DisplayGoodROIs,...
                                'FontSize',12,'FontName','Arial','FontWeight','b',...
                                'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                                'HorizontalAlignment','center',...
                                'Callback',@(~,~)obj.PostProc_DisplayGoodROIs_CB,'Parent',obj.Handles.MainFigure);
                            obj.Handles.UIElements.DisplayBadROIs = uicontrol('Style','checkbox',...
                                'Units','Normalized','Position',[ ...
                                obj.Handles.UIElements.MovieDisplay(1).Box.Position(1)+obj.Handles.UIElements.MovieDisplay(1).Box.Position(3)+0.0025...
                                obj.Handles.UIElements.MovieDisplay(1).Box.Position(2)+0.6*obj.Handles.UIElements.MovieDisplay(1).Box.Position(4)...
                                0.05 0.025],...
                                'String', 'Bad',...
                                'Value',obj.Handles.Values.PostProcessing.DisplayBadROIs,...
                                'FontSize',12,'FontName','Arial','FontWeight','b',...
                                'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                                'HorizontalAlignment','center',...
                                'Callback',@(~,~)obj.PostProc_DisplayBadROIs_CB,'Parent',obj.Handles.MainFigure);

                            if obj.Handles.Values.PostProcessing.DisplayGoodROIs
                                obj.Handles.Values.PostProcessing.DisplayGoodROIs = 'on';
                            else
                                obj.Handles.Values.PostProcessing.DisplayGoodROIs = 'off';
                            end
                            if obj.Handles.Values.PostProcessing.DisplayBadROIs
                                obj.Handles.Values.PostProcessing.DisplayBadROIs = 'on';
                            else
                                obj.Handles.Values.PostProcessing.DisplayBadROIs = 'off';
                            end
                            IndxGood = find(obj.Parameters.CNMFE.Categories==3);
                            IndxBad = find(obj.Parameters.CNMFE.Categories==1);
                            obj.Handles.GoodROIs_HandlesMP = arrayfun(@(x) fill(obj.Handles.Values.Contours{x}(1,:),obj.Handles.Values.Contours{x}(2,:),[0.6 0.6 0.6],'FaceColor','g','FaceAlpha',0.5,'Visible',obj.Handles.Values.PostProcessing.DisplayGoodROIs,'PickableParts','none','EdgeColor','none','LineWidth',1,'Parent',obj.Handles.UIElements.MovieDisplay(3).Box),IndxGood);
                            obj.Handles.GoodROIs_HandlesMovie = arrayfun(@(x) fill(obj.Parameters.CNMFE.Crop(1)+obj.Handles.Values.Contours{x}(1,:),obj.Parameters.CNMFE.Crop(2)+obj.Handles.Values.Contours{x}(2,:),[0.6 0.6 0.6],'FaceColor','g','FaceAlpha',0.5,'Visible',obj.Handles.Values.PostProcessing.DisplayGoodROIs,'PickableParts','none','EdgeColor','none','LineWidth',1,'Parent',obj.Handles.UIElements.MovieDisplay(1).Box),IndxGood);

                            obj.Handles.BadROIs_HandlesMP = arrayfun(@(x) fill(obj.Handles.Values.Contours{x}(1,:),obj.Handles.Values.Contours{x}(2,:),[0.6 0.6 0.6],'FaceColor','r','FaceAlpha',0.5,'Visible',obj.Handles.Values.PostProcessing.DisplayGoodROIs,'PickableParts','none','EdgeColor','none','LineWidth',1,'Parent',obj.Handles.UIElements.MovieDisplay(3).Box),IndxBad);
                            obj.Handles.BadROIs_HandlesMovie = arrayfun(@(x) fill(obj.Parameters.CNMFE.Crop(1)+obj.Handles.Values.Contours{x}(1,:),obj.Parameters.CNMFE.Crop(2)+obj.Handles.Values.Contours{x}(2,:),[0.6 0.6 0.6],'FaceColor','r','FaceAlpha',0.5,'Visible',obj.Handles.Values.PostProcessing.DisplayGoodROIs,'PickableParts','none','EdgeColor','none','LineWidth',1,'Parent',obj.Handles.UIElements.MovieDisplay(1).Box),IndxBad);
                        end

                        obj.Handles.Numbers_HandlesMovie =  arrayfun(@(x) text(obj.Parameters.CNMFE.Crop(1)+nanmean(obj.Handles.Values.Contours{x}(1,:)),obj.Parameters.CNMFE.Crop(2)+nanmean(obj.Handles.Values.Contours{x}(2,:)),num2str(x),'FontName','Arial','FontSize', 10,'Color','g','FontWeight','bold','Parent',obj.Handles.UIElements.MovieDisplay(1).Box,'HorizontalAlignment','center','Clipping','on','Visible',DisplayOpt,'Clipping','on','HitTest','on','ButtonDownFcn',{@(src,evt)obj.ContourClick(src,evt)}),1:numel(obj.Handles.Values.Contours));
                        obj.Handles.Numbers_HandlesMP =  arrayfun(@(x) text(nanmean(obj.Handles.Values.Contours{x}(1,:)),nanmean(obj.Handles.Values.Contours{x}(2,:)),num2str(x),'FontName','Arial','FontSize',8,'Color','g','FontWeight','bold','Parent',obj.Handles.UIElements.MovieDisplay(3).Box,'HorizontalAlignment','center','Clipping','on','Visible',DisplayOpt,'Clipping','on','HitTest','on','ButtonDownFcn',{@(src,evt)obj.ContourClick(src,evt)}),1:numel(obj.Handles.Values.Contours));
                    end
                end
            end
        end



        function PostProc_DisplayGoodROIs_CB(obj)
            if obj.Handles.UIElements.DisplayGoodROIs.Value
                set([obj.Handles.GoodROIs_HandlesMP(:)], 'Visible','on');
                set([obj.Handles.GoodROIs_HandlesMovie(:)], 'Visible','on');
                obj.Handles.Values.PostProcessing.DisplayGoodROIs = 'on';
            else
                set([obj.Handles.GoodROIs_HandlesMP(:)], 'Visible','off');
                set([obj.Handles.GoodROIs_HandlesMovie(:)], 'Visible','off');
                obj.Handles.Values.PostProcessing.DisplayGoodROIs = 'off';
            end
            drawnow
        end


        function PostProc_DisplayBadROIs_CB(obj)
            if obj.Handles.UIElements.DisplayBadROIs.Value
                set([obj.Handles.BadROIs_HandlesMP(:)], 'Visible','on');
                set([obj.Handles.BadROIs_HandlesMovie(:)], 'Visible','on');
                obj.Handles.Values.PostProcessing.DisplayBadROIs = 'on';
            else
                set([obj.Handles.BadROIs_HandlesMP(:)], 'Visible','off');
                set([obj.Handles.BadROIs_HandlesMovie(:)], 'Visible','off');
                obj.Handles.Values.PostProcessing.DisplayBadROIs = 'off';
            end
            drawnow
        end
        function PostProc_MergeROIs_CB(obj)
            % Find index of ROIs to merge
            InxdF = find(obj.Handles.Values.SelectedROIs);
            if numel(InxdF)>1
                obj.LastOperation = [];
                % Hide merged ROIs
                NonSelected_Contour_Width = 1;
                NonSelected_Contour_Color = [0.6 0.6 0.6];
                obj.Handles.Values.MergeIndex(InxdF) = true;

                for Clicked = 1 : numel(InxdF)
                    % Deselect: revert plot to default
                    obj.Handles.Contours_HandlesMP(InxdF(Clicked)).LineWidth = NonSelected_Contour_Width;
                    obj.Handles.Contours_HandlesMovie(InxdF(Clicked)).LineWidth = NonSelected_Contour_Width;
                    obj.Handles.Contours_HandlesMP(InxdF(Clicked)).PickableParts = 'none';
                    obj.Handles.Contours_HandlesMP(InxdF(Clicked)).Visible = 'off';
                    obj.Handles.Contours_HandlesMovie(InxdF(Clicked)).PickableParts = 'none';
                    obj.Handles.Contours_HandlesMovie(InxdF(Clicked)).Visible = 'off';
                    obj.Handles.Contours_HandlesMP(InxdF(Clicked)).EdgeColor = NonSelected_Contour_Color;
                    obj.Handles.Contours_HandlesMovie(InxdF(Clicked)).EdgeColor = NonSelected_Contour_Color;
                    % Deselect: remove from selection
                    obj.Handles.Values.SelectedROIs(InxdF(Clicked)) = false;
                end

                % Merge neurons and display (no deletion yet)
                obj.Parameters.CNMFE.Neuron.manual_merge_preview(InxdF);
                obj.Handles.Values.Contours = obj.Parameters.CNMFE.Neuron.get_contours(0.8);
                Traces = obj.Parameters.CNMFE.Neuron.C_raw;
                IndxNew = size(Traces,1);
                Mins = repmat(min(Traces,[],2),1,length(Traces(1,:)));
                Maxs = repmat(max(Traces,[],2),1,length(Traces(1,:)));
                obj.Handles.Values.Traces = smoothdata((Traces - Mins)./(Maxs - Mins),2,'gaussian',5);

                % Plot contours
                obj.Handles.Contours_HandlesMP = [obj.Handles.Contours_HandlesMP fill(obj.Handles.Values.Contours{end}(1,:),obj.Handles.Values.Contours{end}(2,:),obj.Handles.Values.ColorsList(IndxNew,:),'FaceColor','none','EdgeColor',[0.6 0.6 0.6],'LineWidth',1,'ButtonDownFcn',{@(src,evt)obj.ContourClick(src,evt)},'UserData',IndxNew,'Parent',obj.Handles.UIElements.MovieDisplay(3).Box)];
                obj.Handles.Contours_HandlesMovie = [obj.Handles.Contours_HandlesMovie fill(obj.Parameters.CNMFE.Crop(1)+obj.Handles.Values.Contours{end}(1,:),obj.Parameters.CNMFE.Crop(2)+obj.Handles.Values.Contours{end}(2,:),obj.Handles.Values.ColorsList(IndxNew,:),'FaceColor','none','EdgeColor',obj.Handles.Values.ColorsList(IndxNew,:),'LineWidth',1,'ButtonDownFcn',{@(src,evt)obj.ContourClick(src,evt)},'UserData',IndxNew,'Parent',obj.Handles.UIElements.MovieDisplay(1).Box)];

                obj.Handles.Values.SelectedROIs(IndxNew) = true;

                obj.UpdateTraces;

                obj.LastOperation{1} = 'Merging';
                obj.LastOperation{2} = InxdF;
                obj.LastOperation{3} = IndxNew;
            end
        end

        function PostProc_DeleteROIs_CB(obj)
            % Find index of ROIs to delete
            InxdF = find(obj.Handles.Values.SelectedROIs);
            if ~isempty(InxdF)
                obj.LastOperation = [];
                % Hide deleted ROIs
                NonSelected_Contour_Width = 1;
                NonSelected_Contour_Color = [0.6 0.6 0.6];
                obj.Handles.Values.MergeIndex(InxdF) = true;

                for Clicked = 1 : numel(InxdF)
                    % Deselect: revert plot to default
                    obj.Handles.Contours_HandlesMP(InxdF(Clicked)).LineWidth = NonSelected_Contour_Width;
                    obj.Handles.Contours_HandlesMovie(InxdF(Clicked)).LineWidth = NonSelected_Contour_Width;
                    obj.Handles.Contours_HandlesMP(InxdF(Clicked)).PickableParts = 'none';
                    obj.Handles.Contours_HandlesMP(InxdF(Clicked)).Visible = 'off';
                    obj.Handles.Contours_HandlesMovie(InxdF(Clicked)).PickableParts = 'none';
                    obj.Handles.Contours_HandlesMovie(InxdF(Clicked)).Visible = 'off';
                    obj.Handles.Contours_HandlesMP(InxdF(Clicked)).EdgeColor = NonSelected_Contour_Color;
                    obj.Handles.Contours_HandlesMovie(InxdF(Clicked)).EdgeColor = NonSelected_Contour_Color;
                    % Deselect: remove from selection
                    obj.Handles.Values.SelectedROIs(InxdF(Clicked)) = false;
                end
                obj.UpdateTraces;

                obj.LastOperation{1} = 'Deletion';
                obj.LastOperation{2} = InxdF;
            end
        end

        function PostProc_RunIteration_CB(obj)
            % Apply deletions if needed
            if any(obj.Handles.Values.MergeIndex)
                obj.Parameters.CNMFE.Neuron.delete(find(obj.Handles.Values.MergeIndex));
            end
            obj.Handles.Values.MergeIndex = [];
            obj.LastOperation = [];
            obj.Handles.Values.PreUpdate = obj.Parameters.CNMFE.Neuron;

            obj.Parameters.CNMFE.Neuron.update_background_parallel(true);
            obj.Parameters.CNMFE.Neuron.update_spatial_parallel(true);
            obj.Parameters.CNMFE.Neuron.update_temporal_parallel(true);
            obj.Parameters.CNMFE.Neuron.update_spatial_parallel(true);
            obj.Parameters.CNMFE.Neuron.update_temporal_parallel(true);

            % Clean the traces
            if ~isempty(obj.Handles.TracesPlot)
                delete([obj.Handles.TracesPlot{:}])
                obj.Handles.TracesPlot = [];
                if isfield(obj.Handles.UIElements,'TimeLine')
                    delete(obj.Handles.UIElements.TimeLine)
                    obj.Handles.UIElements = rmfield(obj.Handles.UIElements,'TimeLine');
                end
            end

            obj.Handles.Values.Contours = obj.Parameters.CNMFE.Neuron.get_contours(0.8);
            Traces = obj.Parameters.CNMFE.Neuron.C_raw;
            Mins = repmat(min(Traces,[],2),1,length(Traces(1,:)));
            Maxs = repmat(max(Traces,[],2),1,length(Traces(1,:)));
            obj.Handles.Values.Traces = smoothdata((Traces - Mins)./(Maxs - Mins),2,'gaussian',5);
            % Plot contours
            delete(obj.Handles.Contours_HandlesMP(:))
            delete(obj.Handles.Contours_HandlesMovie(:))
            obj.Handles.Contours_HandlesMP = arrayfun(@(x) fill(obj.Handles.Values.Contours{x}(1,:),obj.Handles.Values.Contours{x}(2,:),[0.6 0.6 0.6],'FaceColor','none','EdgeColor',[0.6 0.6 0.6],'LineWidth',1,'ButtonDownFcn',{@(src,evt)obj.ContourClick(src,evt)},'UserData',x,'Parent',obj.Handles.UIElements.MovieDisplay(3).Box),1:numel(obj.Handles.Values.Contours));
            obj.Handles.Contours_HandlesMovie = arrayfun(@(x) fill(obj.Parameters.CNMFE.Crop(1)+obj.Handles.Values.Contours{x}(1,:),obj.Parameters.CNMFE.Crop(2)+obj.Handles.Values.Contours{x}(2,:),[0.6 0.6 0.6],'FaceColor','none','EdgeColor',[0.6 0.6 0.6],'LineWidth',1,'ButtonDownFcn',{@(src,evt)obj.ContourClick(src,evt)},'UserData',x,'Parent',obj.Handles.UIElements.MovieDisplay(1).Box),1:numel(obj.Handles.Values.Contours));
            % Retrieve default plot colors to use when clicking contours
            DefC = DefColors;
            obj.Handles.Values.ColorsList = repmat(DefC,[1e4 1]);
            obj.Handles.Values.SelectedROIs = false(1e4,1);
            obj.Handles.Values.MergeIndex = false(1e4,1);
            obj.LastOperation{1} = 'Iteration';

        end

        function PostProc_CancelLast_CB(obj)
            if ~isempty(obj.LastOperation)
                NonSelected_Contour_Width = 1;
                NonSelected_Contour_Color = [0.6 0.6 0.6];
                Selected_Contour_Width = 2.5;
                switch obj.LastOperation{1}
                    case 'Merging'
                        % Restore the ROIs that were merged together
                        InxdF = obj.LastOperation{2};
                        for Clicked = 1 : numel(InxdF)
                            obj.Handles.Contours_HandlesMP(InxdF(Clicked)).PickableParts = 'all';
                            obj.Handles.Contours_HandlesMP(InxdF(Clicked)).Visible = 'on';
                            obj.Handles.Contours_HandlesMP(InxdF(Clicked)).LineWidth = Selected_Contour_Width;
                            obj.Handles.Contours_HandlesMP(InxdF(Clicked)).EdgeColor = obj.Handles.Values.ColorsList(InxdF(Clicked),:);
                            obj.Handles.Contours_HandlesMovie(InxdF(Clicked)).PickableParts = 'all';
                            obj.Handles.Contours_HandlesMovie(InxdF(Clicked)).Visible = 'on';
                            obj.Handles.Contours_HandlesMovie(InxdF(Clicked)).LineWidth = Selected_Contour_Width;
                            obj.Handles.Contours_HandlesMovie(InxdF(Clicked)).EdgeColor = obj.Handles.Values.ColorsList(InxdF(Clicked),:);
                            % Add to selected list
                            obj.Handles.Values.SelectedROIs(InxdF(Clicked)) = true;
                        end
                        % Remove from the list to delete
                        obj.Handles.Values.MergeIndex(InxdF) = false;

                        % Hide the new merged ROI
                        InxdF = obj.LastOperation{3};
                        obj.Handles.Contours_HandlesMP(InxdF).LineWidth = NonSelected_Contour_Width;
                        obj.Handles.Contours_HandlesMovie(InxdF).LineWidth = NonSelected_Contour_Width;
                        obj.Handles.Contours_HandlesMP(InxdF).PickableParts = 'none';
                        obj.Handles.Contours_HandlesMP(InxdF).Visible = 'off';
                        obj.Handles.Contours_HandlesMovie(InxdF).PickableParts = 'none';
                        obj.Handles.Contours_HandlesMovie(InxdF).Visible = 'off';
                        obj.Handles.Contours_HandlesMP(InxdF).EdgeColor = NonSelected_Contour_Color;
                        obj.Handles.Contours_HandlesMovie(InxdF).EdgeColor = NonSelected_Contour_Color;
                        % Deselect: remove from selection
                        obj.Handles.Values.SelectedROIs(InxdF) = false;

                        % Remove from data: add to the list to delete at the
                        % end
                        obj.Handles.Values.MergeIndex(InxdF) = true;
                        % Remove from selection
                        obj.Handles.Values.SelectedROIs(InxdF) = false;
                        % Update traces
                        obj.UpdateTraces;

                    case 'Deletion'
                        % Restore the ROIs that were deleted
                        InxdF = obj.LastOperation{2};
                        for Clicked = 1 : numel(InxdF)
                            obj.Handles.Contours_HandlesMP(InxdF(Clicked)).PickableParts = 'all';
                            obj.Handles.Contours_HandlesMP(InxdF(Clicked)).Visible = 'on';
                            obj.Handles.Contours_HandlesMP(InxdF(Clicked)).LineWidth = Selected_Contour_Width;
                            obj.Handles.Contours_HandlesMP(InxdF(Clicked)).EdgeColor = obj.Handles.Values.ColorsList(InxdF(Clicked),:);
                            obj.Handles.Contours_HandlesMovie(InxdF(Clicked)).PickableParts = 'all';
                            obj.Handles.Contours_HandlesMovie(InxdF(Clicked)).Visible = 'on';
                            obj.Handles.Contours_HandlesMovie(InxdF(Clicked)).LineWidth = Selected_Contour_Width;
                            obj.Handles.Contours_HandlesMovie(InxdF(Clicked)).EdgeColor = obj.Handles.Values.ColorsList(InxdF(Clicked),:);
                            % Add to selected list
                            obj.Handles.Values.SelectedROIs(InxdF(Clicked)) = true;
                        end
                        % Remove from the list to delete
                        obj.Handles.Values.MergeIndex(InxdF) = false;
                        obj.UpdateTraces;

                    case 'Iteration'
                        % Restore the previous object
                        obj.Parameters.CNMFE.Neuron = obj.Handles.Values.PreUpdate;

                        % Clean the traces
                        if ~isempty(obj.Handles.TracesPlot)
                            delete([obj.Handles.TracesPlot{:}])
                            obj.Handles.TracesPlot = [];
                            if isfield(obj.Handles.UIElements,'TimeLine')
                                delete(obj.Handles.UIElements.TimeLine)
                                obj.Handles.UIElements = rmfield(obj.Handles.UIElements,'TimeLine');
                            end
                        end

                        % Re-derive the traces & contours
                        obj.Handles.Values.Contours = obj.Parameters.CNMFE.Neuron.get_contours(0.8);
                        Traces = obj.Parameters.CNMFE.Neuron.C_raw;
                        Mins = repmat(min(Traces,[],2),1,length(Traces(1,:)));
                        Maxs = repmat(max(Traces,[],2),1,length(Traces(1,:)));
                        obj.Handles.Values.Traces = smoothdata((Traces - Mins)./(Maxs - Mins),2,'gaussian',10);
                        % Plot contours
                        delete(obj.Handles.Contours_HandlesMP(:))
                        delete(obj.Handles.Contours_HandlesMovie(:))
                        obj.Handles.Contours_HandlesMP = arrayfun(@(x) fill(obj.Handles.Values.Contours{x}(1,:),obj.Handles.Values.Contours{x}(2,:),[0.6 0.6 0.6],'FaceColor','none','EdgeColor',[0.6 0.6 0.6],'LineWidth',1,'ButtonDownFcn',{@(src,evt)obj.ContourClick(src,evt)},'UserData',x,'Parent',obj.Handles.UIElements.MaxProjection.Plot),1:numel(obj.Handles.Values.Contours));
                        obj.Handles.Contours_HandlesMovie = arrayfun(@(x) fill(obj.Parameters.CNMFE.Crop(1)+obj.Handles.Values.Contours{x}(1,:),obj.Parameters.CNMFE.Crop(2)+obj.Handles.Values.Contours{x}(2,:),[0.6 0.6 0.6],'FaceColor','none','EdgeColor',[0.6 0.6 0.6],'LineWidth',1,'ButtonDownFcn',{@(src,evt)obj.ContourClick(src,evt)},'UserData',x,'Parent',obj.Handles.UIElements.MovieDisplay(1).Box),1:numel(obj.Handles.Values.Contours));
                        % Retrieve default plot colors to use when clicking contours
                        DefC = DefColors;
                        obj.Handles.Values.ColorsList = repmat(DefC,[1e4 1]);
                        obj.Handles.Values.SelectedROIs = false(1e4,1);
                        obj.Handles.Values.MergeIndex = false(1e4,1);
                end
                obj.LastOperation = [];
            end
        end

        function PostProc_Restore(obj)
            Answer = questdlg(['Are you sure you wish to restore the intial CNMFE output?' newline newline,...
                '(ALL post-processing will be lost)'],'Please choose...',...
                'YES, revert to the raw output (all post-processing will be gone).',...
                'NO, continue post-processing from the current state.',...
                'NO, continue post-processing from the current state.');
            if strcmpi(Answer,'NO, continue post-processing from the current state.')
                return
            else
                obj.Parameters.CNMFE.Neuron = obj.Parameters.CNMFE.VanillaNeuron;
                obj.PostProcessing_Layout(true);
                obj.Handles.Values.MergeIndex = [];
                obj.LastOperation = [];
            end
        end


        function PostProc_Save_CB(obj)
            % Apply deletions if needed
            if any(obj.Handles.Values.MergeIndex)
                obj.Parameters.CNMFE.Neuron.delete(find(obj.Handles.Values.MergeIndex));
            end
            obj.Handles.Values.MergeIndex = [];
            obj.LastOperation = [];

            % Save the workspace for future analysis
            obj.Parameters.CNMFE.Neuron.orderROIs('snr');

            % Update plotting
            obj.Handles.Values.Contours = obj.Parameters.CNMFE.Neuron.get_contours(0.8);
            Traces = obj.Parameters.CNMFE.Neuron.C_raw;
            Mins = repmat(min(Traces,[],2),1,length(Traces(1,:)));
            Maxs = repmat(max(Traces,[],2),1,length(Traces(1,:)));
            obj.Handles.Values.Traces = smoothdata((Traces - Mins)./(Maxs - Mins),2,'gaussian',10);
            % Plot contours
            delete(obj.Handles.Contours_HandlesMP(:))
            delete(obj.Handles.Contours_HandlesMovie(:))
            obj.Handles.Contours_HandlesMP = arrayfun(@(x) fill(obj.Handles.Values.Contours{x}(1,:),obj.Handles.Values.Contours{x}(2,:),[0.6 0.6 0.6],'FaceColor','none','EdgeColor',[0.6 0.6 0.6],'LineWidth',1,'ButtonDownFcn',{@(src,evt)obj.ContourClick(src,evt)},'UserData',x,'Parent',obj.Handles.UIElements.MovieDisplay(3).Box),1:numel(obj.Handles.Values.Contours));
            obj.Handles.Contours_HandlesMovie = arrayfun(@(x) fill(obj.Parameters.CNMFE.Crop(1)+obj.Handles.Values.Contours{x}(1,:),obj.Parameters.CNMFE.Crop(2)+obj.Handles.Values.Contours{x}(2,:),[0.6 0.6 0.6],'FaceColor','none','EdgeColor',[0.6 0.6 0.6],'LineWidth',1,'ButtonDownFcn',{@(src,evt)obj.ContourClick(src,evt)},'UserData',x,'Parent',obj.Handles.UIElements.MovieDisplay(1).Box),1:numel(obj.Handles.Values.Contours));
            % Retrieve default plot colors to use when clicking contours
            DefC = DefColors;
            obj.Handles.Values.ColorsList = repmat(DefC,[1e4 1]);
            obj.Handles.Values.SelectedROIs = false(1e4,1);
            obj.Handles.Values.MergeIndex = false(1e4,1);

            obj.UpdateTraces;
            obj.Parameters.CNMFE.Date = datetime;

            MetaFileCurrent = [...
                obj.Parameters.CNMFE.CNMFEFolder,...
                obj.Parameters.PreProcessing.Basename,...
                '_PP-' num2str(obj.Parameters.PreProcessing.ForkNumber),...
                '_PF-' num2str(obj.Parameters.PreFiltering.ForkNumber),...
                '_MC-' num2str(obj.Parameters.MotionCorrection.ForkNumber),...
                '_CNMFE-' num2str(obj.Parameters.CNMFE.ForkNumber) '_StepMeta.mat'];
            PreviousMeta = load(MetaFileCurrent);
            PreviousMeta.CNMFE = obj.Parameters.CNMFE;
            save(MetaFileCurrent,'-struct','PreviousMeta')

        end

        function SliderTracesCB(obj)
            if ~obj.Dragging
                obj.Dragging = true;
                obj.Handles.MainFigure.WindowButtonMotionFcn = @(~,~)obj.MovingTracesSlider;
                obj.Handles.MainFigure.WindowButtonUpFcn = @(~,~)obj.SliderTracesCB;
            else
                obj.Dragging = false;
                obj.Handles.MainFigure.WindowButtonMotionFcn = [];
                obj.Handles.MainFigure.WindowButtonUpFcn = [];
            end
        end

        function MovingTracesSlider(obj)
            CurrentCursor = obj.Handles.UIElements.Traces.CurrentPoint;
            if CurrentCursor(1)>= obj.Handles.UIElements.Traces.XLim(1) &&  CurrentCursor(1)<=obj.Handles.UIElements.Traces.XLim(2),
                obj.CurrentTime = CurrentCursor(1);
                CurrentPlayingStatus = obj.Playing;
                obj.Playing = false;
                obj.PlayingSingle = true;
                drawnow
                obj.PlayMovies;
                drawnow
                pause(0.01);
                obj.Playing = CurrentPlayingStatus;
                obj.PlayMovies;
            end
        end

        function PostProc_CancelSelection_CB(obj)

            NonSelected_Contour_Width = 1;
            NonSelected_Contour_Color = [0.6 0.6 0.6];
            if ~isempty(obj.Handles.TracesPlot)
                if numel(obj.Handles.TracesPlot)>1 || iscell(obj.Handles.TracesPlot)
                    delete([obj.Handles.TracesPlot{:}])
                else
                    delete(obj.Handles.TracesPlot)
                end
                delete(obj.Handles.UIElements.TimeLine)
                obj.Handles.UIElements = rmfield(obj.Handles.UIElements,'TimeLine');
            end
            InxdF = find(obj.Handles.Values.SelectedROIs);
            for Clicked = 1 : numel(InxdF)
                % Deselect: revert plot to default
                obj.Handles.Contours_HandlesMP(InxdF(Clicked)).LineWidth = NonSelected_Contour_Width;
                obj.Handles.Contours_HandlesMovie(InxdF(Clicked)).LineWidth = NonSelected_Contour_Width;
                obj.Handles.Contours_HandlesMP(InxdF(Clicked)).EdgeColor = NonSelected_Contour_Color;
                obj.Handles.Contours_HandlesMovie(InxdF(Clicked)).EdgeColor = NonSelected_Contour_Color;
                % Deselect: remove from selection
                obj.Handles.Values.SelectedROIs(InxdF(Clicked)) = false;
            end
        end

        function UpdateTraces(obj)
            hold(obj.Handles.UIElements.Traces,'on');
            PrePlot = false;
            if ~isempty(obj.Handles.TracesPlot)
                PrePlot = true;
                delete([obj.Handles.TracesPlot{:}])
                obj.Handles.TracesPlot = [];
            end
            if isfield(obj.Handles.UIElements,'TimeLine')
                delete(obj.Handles.UIElements.TimeLine)
                obj.Handles.UIElements = rmfield(obj.Handles.UIElements,'TimeLine');
            end
            if any(obj.Handles.Values.SelectedROIs)
                IndxToPlot = find(obj.Handles.Values.SelectedROIs);
                obj.Handles.TracesPlot = arrayfun(@(x) plot(obj.Parameters.PreFiltering.TimeStamps,obj.Handles.Values.Traces(IndxToPlot(x),:) + x - 0.5,'LineWidth',1.5,'Color',obj.Handles.Values.ColorsList(IndxToPlot(x),:),'Parent',obj.Handles.UIElements.Traces),1:numel(IndxToPlot),'UniformOutput',false);
                if ~PrePlot
                    obj.Handles.UIElements.Traces.XLim = obj.Parameters.PreFiltering.TimeStamps([1 end]);
                end
                obj.Handles.UIElements.Traces.YLimMode = 'auto';
                obj.Handles.UIElements.Traces.YTick = [];
                obj.Handles.UIElements.TimeLine = plot(obj.CurrentTime * [1 1],obj.Handles.UIElements.Traces.YLim,'k','Parent',obj.Handles.UIElements.Traces,'LineWidth',2,'ButtonDownFcn',{@(~,~)obj.SliderTracesCB});
            elseif ~isfield(obj.Handles.UIElements,'TimeLine') || isempty(obj.Handles.UIElements.TimeLine)
                if ~PrePlot
                    obj.Handles.UIElements.Traces.XLim = obj.Parameters.PreFiltering.TimeStamps([1 end]);
                end
                obj.Handles.UIElements.Traces.YLimMode = 'auto';
                obj.Handles.UIElements.Traces.YTick = [];
                obj.Handles.UIElements.TimeLine = plot(obj.CurrentTime * [1 1],obj.Handles.UIElements.Traces.YLim,'k','Parent',obj.Handles.UIElements.Traces,'LineWidth',2,'ButtonDownFcn',{@(~,~)obj.SliderTracesCB});
            end
        end

        function ContourClick(obj,src,~)
            Clicked = src.UserData;
            NonSelected_Contour_Width = 1;
            NonSelected_Contour_Color = [0.6 0.6 0.6];
            Selected_Contour_Width = 2.5;

            % Check if already selected or not
            if obj.Handles.Values.SelectedROIs(Clicked)
                % Deselect: revert plot to default
                obj.Handles.Contours_HandlesMP(Clicked).LineWidth = NonSelected_Contour_Width;
                obj.Handles.Contours_HandlesMovie(Clicked).LineWidth = NonSelected_Contour_Width;
                obj.Handles.Contours_HandlesMP(Clicked).EdgeColor = NonSelected_Contour_Color;
                obj.Handles.Contours_HandlesMovie(Clicked).EdgeColor = NonSelected_Contour_Color;
                % Deselect: remove from selection
                obj.Handles.Values.SelectedROIs(Clicked) = false;
            else
                % Update contour plot
                obj.Handles.Contours_HandlesMP(Clicked).LineWidth = Selected_Contour_Width;
                obj.Handles.Contours_HandlesMovie(Clicked).LineWidth = Selected_Contour_Width;
                obj.Handles.Contours_HandlesMP(Clicked).EdgeColor = obj.Handles.Values.ColorsList(Clicked,:);
                obj.Handles.Contours_HandlesMovie(Clicked).EdgeColor = obj.Handles.Values.ColorsList(Clicked,:);
                % Add to selection
                obj.Handles.Values.SelectedROIs(Clicked) = true;
            end
            % Update traces
            obj.UpdateTraces;
        end

        function PostProc_Filtering_CB(obj)
            drawnow
            obj.Handles.Values.PostProcessing.Movie1Filtered = obj.Handles.UIElements.Filtering.Value;
            obj.CurrentPlayers(1).Filter = obj.Handles.Values.PostProcessing.Movie1Filtered;
            if obj.Handles.UIElements.Filtering.Value
                obj.Handles.UIElements.Sig1PP_Edit.Enable = 'on';
                obj.Handles.UIElements.Sig2PP_Edit.Enable = 'on';
                obj.Handles.UIElements.WindowPP_Edit.Enable = 'on';
            else
                obj.Handles.UIElements.Sig1PP_Edit.Enable = 'off';
                obj.Handles.UIElements.Sig2PP_Edit.Enable = 'off';
                obj.Handles.UIElements.WindowPP_Edit.Enable = 'off';
            end
            obj.Handles.UIElements.MovieDisplay(1).Plot.CLimMode = 'manual';
            obj.UpdateDisplay;
            if obj.Handles.Values.PostProcessing.AutoCLim
                obj.Handles.UIElements.MovieDisplay(1).Plot.CLimMode = 'auto';
                if obj.Handles.Values.PostProcessing.Movie1Filtered
                    obj.Handles.Values.CaMovieFilteredCLim = obj.Handles.UIElements.MovieDisplay(1).Plot.CLim;
                    obj.Handles.Values.CaMovieFilteredInit = true;
                    obj.Handles.UIElements.ThresholdLowPostProc_Edit.String = obj.Handles.Values.CaMovieFilteredCLim(1);
                    obj.Handles.UIElements.ThresholdHighPostProc_Edit.String = obj.Handles.Values.CaMovieFilteredCLim(2);
                else
                    obj.Handles.Values.CaMovieCLim = obj.Handles.UIElements.MovieDisplay(1).Plot.CLim;
                    obj.Handles.Values.CaMovieInit = true;
                    obj.Handles.UIElements.ThresholdLowPostProc_Edit.String = obj.Handles.Values.CaMovieCLim(1);
                    obj.Handles.UIElements.ThresholdHighPostProc_Edit.String = obj.Handles.Values.CaMovieCLim(2);
                end
            else
                obj.Handles.UIElements.MovieDisplay(1).Plot.CLimMode = 'manual';
                if obj.Handles.Values.PostProcessing.Movie1Filtered
                    if obj.Handles.Values.CaMovieFilteredInit
                        obj.Handles.UIElements.MovieDisplay(1).Plot.CLim = obj.Handles.Values.CaMovieFilteredCLim;
                    else
                        obj.Handles.UIElements.MovieDisplay(1).Plot.CLimMode = 'auto';
                        obj.Handles.Values.CaMovieFilteredCLim = obj.Handles.UIElements.MovieDisplay(1).Plot.CLim;
                        obj.Handles.Values.CaMovieFilteredInit = true;
                        obj.Handles.UIElements.MovieDisplay(1).Plot.CLimMode = 'manual';
                    end
                    obj.Handles.UIElements.ThresholdLowPostProc_Edit.String = obj.Handles.Values.CaMovieFilteredCLim(1);
                    obj.Handles.UIElements.ThresholdHighPostProc_Edit.String = obj.Handles.Values.CaMovieFilteredCLim(2);
                else
                    if obj.Handles.Values.CaMovieInit
                        obj.Handles.UIElements.MovieDisplay(1).Plot.CLim = obj.Handles.Values.CaMovieCLim;
                    else
                        obj.Handles.UIElements.MovieDisplay(1).Plot.CLimMode = 'auto';
                        obj.Handles.Values.CaMovieCLim = obj.Handles.UIElements.MovieDisplay(1).Plot.CLim;
                        obj.Handles.Values.CaMovieInit = true;
                        obj.Handles.UIElements.MovieDisplay(1).Plot.CLimMode = 'manual';
                    end
                    obj.Handles.UIElements.ThresholdLowPostProc_Edit.String = obj.Handles.Values.CaMovieCLim(1);
                    obj.Handles.UIElements.ThresholdHighPostProc_Edit.String = obj.Handles.Values.CaMovieCLim(2);
                end
            end
        end

        function PostProc_AutoCLim_CB(obj)
            if obj.Handles.UIElements.AutoCLim.Value
                obj.Handles.Values.PostProcessing.AutoCLim = true;
                obj.Handles.UIElements.MovieDisplay(1).Plot.CLimMode = 'auto';
                obj.Handles.UIElements.ThresholdLowPostProc_Edit.Enable = 'off';
                obj.Handles.UIElements.ThresholdHighPostProc_Edit.Enable = 'off';
            else
                obj.Handles.Values.PostProcessing.AutoCLim = false;
                obj.Handles.UIElements.MovieDisplay(1).Plot.CLimMode = 'manual';
                obj.Handles.UIElements.ThresholdLowPostProc_Edit.Enable = 'on';
                obj.Handles.UIElements.ThresholdHighPostProc_Edit.Enable = 'on';
            end
            obj.UpdateDisplay;
        end

        function UpdatePFKernelPP(obj)
            if obj.Handles.Values.PostProcessing.Filtering.Sig1>0 && obj.Handles.Values.PostProcessing.Filtering.Sig2>0
                gaussian1 = fspecial('Gaussian', obj.Handles.Values.PostProcessing.Filtering.WindowSize, obj.Handles.Values.PostProcessing.Filtering.Sig1);
                gaussian2 = fspecial('Gaussian', obj.Handles.Values.PostProcessing.Filtering.WindowSize, obj.Handles.Values.PostProcessing.Filtering.Sig2);
                if obj.Handles.Values.PostProcessing.Filtering.Sig1>obj.Handles.Values.PostProcessing.Filtering.Sig2
                    psf = gaussian1 - gaussian2;
                else
                    psf = gaussian2 - gaussian1;
                end
            elseif obj.Handles.Values.PostProcessing.Filtering.Sig1 == 0
                psf = fspecial('Gaussian', obj.Handles.Values.PostProcessing.Filtering.WindowSize, obj.Handles.Values.PostProcessing.Filtering.Sig2);
            else
                psf = fspecial('Gaussian', obj.Handles.Values.PostProcessing.Filtering.WindowSize, obj.Handles.Values.PostProcessing.Filtering.Sig1);
            end
            ind_nonzero = (psf(:)>=max(psf(:,1)));
            psf = psf-mean(psf(ind_nonzero));
            psf(~ind_nonzero) = 0;
            obj.Handles.Values.PostProcessing.Filtering.Kernel = psf;

        end

        function Sig1PP_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.Sig1PP_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>0
                obj.Handles.Values.PostProcessing.Filtering.Sig1 = TempInput;
                obj.UpdatePFKernelPP;
            elseif TempInput == 0
                if obj.Handles.Values.PostProcessing.Filtering.Sig2 == 0
                    obj.Handles.UIElements.Sig1PP_Edit.String = num2str(obj.Handles.Values.PostProcessing.Filtering.Sig1);
                else
                    obj.Handles.Values.PostProcessing.Filtering.Sig1 = TempInput;
                    obj.UpdatePFKernelPP;
                end
            else
                obj.Handles.UIElements.Sig1PP_Edit.String = num2str(obj.Handles.Values.PostProcessing.Filtering.Sig1);
            end
            obj.UpdateDisplay;
        end

        function Sig2PP_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.Sig2PP_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>=0
                obj.Handles.Values.PostProcessing.Filtering.Sig2 = TempInput;
                obj.UpdatePFKernelPP;
            elseif TempInput == 0
                if obj.Handles.Values.PostProcessing.Filtering.Sig1 == 0
                    obj.Handles.UIElements.Sig2_EditPP.String = num2str(obj.Handles.Values.PostProcessing.Filtering.Sig2);
                else
                    obj.Handles.Values.PostProcessing.Filtering.Sig2 = TempInput;
                    obj.UpdatePFKernel;
                end
            else
                obj.Handles.UIElements.Sig2PP_Edit.String = num2str(obj.Handles.Values.PostProcessing.Filtering.Sig2);
            end
            obj.UpdateDisplay;
        end

        function WindowPP_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.WindowPP_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>=1
                obj.Handles.Values.PostProcessing.Filtering.WindowSize = TempInput;
                obj.UpdatePFKernelPP;
            else
                obj.Handles.UIElements.WindowPP_Edit.String = num2str(obj.Handles.Values.PostProcessing.Filtering.WindowSize);
            end
            obj.UpdateDisplay;
        end

        function ThresholdLowPostProc_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.ThresholdLowPostProc_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput<obj.Handles.UIElements.MovieDisplay(1).Plot.CLim(2),
                if obj.Handles.Values.PostProcessing.Movie1Filtered
                    obj.Handles.Values.CaMovieFilteredCLim(1) = TempInput;
                else
                    obj.Handles.Values.CaMovieCLim(1) = TempInput;
                end
                obj.Handles.UIElements.MovieDisplay(1).Plot.CLim(1) = TempInput;
            else
                if obj.Handles.Values.PostProcessing.Movie1Filtered
                    obj.Handles.UIElements.ThresholdLowPostProc_Edit.String = num2str(obj.Handles.Values.CaMovieFilteredCLim(1));
                else
                    obj.Handles.UIElements.ThresholdLowPostProc_Edit.String = num2str(obj.Handles.Values.CaMovieCLim(1));
                end
            end
        end

        function ThresholdHighPostProc_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.ThresholdHighPostProc_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>obj.Handles.UIElements.MovieDisplay(1).Plot.CLim(1)
                if obj.Handles.Values.PostProcessing.Movie1Filtered
                    obj.Handles.Values.CaMovieFilteredCLim(2) = TempInput;
                else
                    obj.Handles.Values.CaMovieCLim(2) = TempInput;
                end
                obj.Handles.UIElements.MovieDisplay(1).Plot.CLim(2) = TempInput;
            else
                if obj.Handles.Values.PostProcessing.Movie1Filtered
                    obj.Handles.UIElements.ThresholdHighPostProc_Edit.String = num2str(obj.Handles.Values.CaMovieFilteredCLim(2));
                else
                    obj.Handles.UIElements.ThresholdHighPostProc_Edit.String = num2str(obj.Handles.Values.CaMovieCLim(2));
                end
            end
        end

        function PostProc_Extract_CB(obj)
            % Prompt to confirm
            Answer = questdlg(['\fontsize{13}Are you sure you wish to export the current post-processed data to be used for analyses?' newline newline,...
                'This will:' newline,...
                '      1) Export the data to the session folder as CalciumData.mat file' newline,...
                '      2) Tag this particular branch of the analysis to be transfered to the longterm storage' newline,...
                '      3) Delete permanently the other branches for this particular session (i.e. the other extractions that were tried with other parameters)' newline,...
                ],'Please choose...',...
                'YES, export this version and delete the others.',...
                'NO, continue post-processing here, but choose later which one to keep.',...
                struct('Default','NO, continue post-processing here, but choose later which one to keep.','Interpreter','tex'));
            if strcmpi(Answer,'NO, continue post-processing here, but choose later which one to keep.')
                return
            else
                % As a safety measure: check that we don't have another
                % file that was tagged, and/or another exported file
                CNMFEFiles = dir([obj.Parameters.PreProcessing.MainFolder '**' filesep '*CNMFE*_StepMeta.mat']);
                CNMFEFiles = arrayfun(@(x) [CNMFEFiles(x).folder filesep CNMFEFiles(x).name],1:numel(CNMFEFiles),'UniformOutput',false);
                % We're here, so there is at least one file to find, we can
                % loop through. If the current branch was already tagged,
                % we can just prompt as for the others but with the extra
                % info
                MetaFileCurrent = [...
                    obj.Parameters.CNMFE.CNMFEFolder,...
                    obj.Parameters.PreProcessing.Basename,...
                    '_PP-' num2str(obj.Parameters.PreProcessing.ForkNumber),...
                    '_PF-' num2str(obj.Parameters.PreFiltering.ForkNumber),...
                    '_MC-' num2str(obj.Parameters.MotionCorrection.ForkNumber),...
                    '_CNMFE-' num2str(obj.Parameters.CNMFE.ForkNumber) '_StepMeta.mat'];
                CalciumDataFile = [obj.Parameters.PreProcessing.MainFolder obj.Parameters.PreProcessing.Basename '_CalciumData.mat'];

                Proceed = true;
                for CF = 1 : numel(CNMFEFiles)
                    Loaded = load(CNMFEFiles{CF});
                    if isfield(Loaded.CNMFE,'Extracted')
                        if Loaded.CNMFE.Extracted
                            % We expect to find the extracted file, but
                            % just in case... let's check
                            if exist(CalciumDataFile,'file')==2
                                if strcmpi(CNMFEFiles{CF},MetaFileCurrent)
                                    Answer = questdlg(['The current branch was already extracted at some point. Do you wish to re-extract it?' newline newline],...
                                        'Please choose...',...
                                        'YES, reexport this version.',...
                                        'NO, continue post-processing here, but choose later to reexport or not.',...
                                        'NO, continue post-processing here, but choose later to reexport or not.');
                                    if ~strcmpi(Answer,'YES, reexport this version.')
                                        Proceed = false;
                                    end
                                else
                                    Answer = questdlg(['A different branch was already extracted at some point. Do you wish replace the data with the current branch?' newline newline],...
                                        'Please choose...',...
                                        'YES, export this version instead of the other.',...
                                        'NO, continue post-processing here, and choose later to reexport or not.',...
                                        'NO, continue post-processing here, and choose later to reexport or not.');
                                    if strcmpi(Answer,'YES, export this version instead of the other.')
                                        % Change the tag of the other file
                                        Loaded.CNMFE.Extracted = false;
                                        save(CNMFEFiles{CF},'-struct','Loaded')
                                    else
                                        Proceed = false;
                                    end
                                end
                            else
                                Answer = questdlg(['A different branch (' Loaded.CNMFE.CNMFEFolder(numel(Loaded.PreProcessing.MainFolder):end) ') was already extracted at some point and is marked as such, but the actual exported file cannot be found.' newline,...
                                    'Do you still wish to export the current branch?' newline,...
                                    ],...
                                    'Please choose...',...
                                    'YES, export this version .',...
                                    'NO, continue post-processing here, and choose later to export or not.',...
                                    'NO, continue post-processing here, and choose later to export or not.');
                                if strcmpi(Answer,'YES, export this version .')
                                    % Change the tag of the other file
                                    Loaded.CNMFE.Extracted = false;
                                    save(CNMFEFiles{CF},'-struct','Loaded')
                                else
                                    Proceed = false;
                                end
                            end
                            break
                        end
                    end
                end

                if Proceed
                    % Make sure the changes were all applied, set the extract string to true and resave the file
                    obj.Parameters.CNMFE.Extracted = true;
                    obj.Parameters.CNMFE.ExtractionDate = datetime;
                    obj.PostProc_Save_CB;
                    drawnow; obj.DisableAll; drawnow

                    % Delete all the other branches, step by step
                    % Try to close all files
                    fclose('all');
                    % Preprocessing
                    PPFolders = dir([obj.Parameters.PreProcessing.MainFolder 'PP-*' ]);
                    PPFolders = PPFolders([PPFolders.isdir]);
                    PPFolders = arrayfun(@(x) [PPFolders(x).folder filesep PPFolders(x).name filesep],1:numel(PPFolders),'UniformOutput',false);
                    if isempty(PPFolders)
                        error('The folder structure has been manually edited and corrupted. Aborting.')
                    else
                        for P = 1 : numel(PPFolders)
                            if ~strcmpi(PPFolders{P},obj.Parameters.PreProcessing.PreProcessingFolder)
                                St = rmdir(PPFolders{P},'s');
                                if ~St
                                    warning(['Some folders could not be deleted, probably because some file is still open.' newline 'Aborting.'])
                                    return
                                end
                            end
                        end
                    end

                    % Prefiltering
                    PFFolders = dir([obj.Parameters.PreProcessing.PreProcessingFolder 'PF-*' ]);
                    PFFolders = PFFolders([PFFolders.isdir]);
                    PFFolders = arrayfun(@(x) [PFFolders(x).folder filesep PFFolders(x).name filesep],1:numel(PFFolders),'UniformOutput',false);
                    if isempty(PFFolders)
                        error('The folder structure has been manually edited and corrupted. Aborting.')
                    else
                        for P = 1 : numel(PFFolders)
                            if ~strcmpi(PFFolders{P},obj.Parameters.PreFiltering.PreFilteringFolder)
                                St = rmdir(PFFolders{P},'s');
                                if ~St
                                    warning(['Some folders could not be deleted, probably because some file is still open.' newline 'Aborting.'])
                                    return
                                end
                            end
                        end
                    end

                    % Motion correction
                    MCFolders = dir([obj.Parameters.PreFiltering.PreFilteringFolder 'MC-*' ]);
                    MCFolders = MCFolders([MCFolders.isdir]);
                    MCFolders = arrayfun(@(x) [MCFolders(x).folder filesep MCFolders(x).name filesep],1:numel(MCFolders),'UniformOutput',false);
                    if isempty(MCFolders)
                        error('The folder structure has been manually edited and corrupted. Aborting.')
                    else
                        for P = 1 : numel(MCFolders)
                            if ~strcmpi(MCFolders{P},obj.Parameters.MotionCorrection.MotionCorrectionFolder)
                                St = rmdir(MCFolders{P},'s');
                                if ~St
                                    warning(['Some folders could not be deleted, probably because some file is still open.' newline 'Aborting.'])
                                    return
                                end
                            end
                        end
                    end

                    % CNMFE
                    CNMFEFolders = dir([obj.Parameters.MotionCorrection.MotionCorrectionFolder 'CNMFE-*' ]);
                    CNMFEFolders = CNMFEFolders([CNMFEFolders.isdir]);
                    CNMFEFolders = arrayfun(@(x) [CNMFEFolders(x).folder filesep CNMFEFolders(x).name filesep],1:numel(CNMFEFolders),'UniformOutput',false);
                    if isempty(CNMFEFolders)
                        error('The folder structure has been manually edited and corrupted. Aborting.')
                    else
                        for P = 1 : numel(CNMFEFolders)
                            if ~strcmpi(CNMFEFolders{P},obj.Parameters.CNMFE.CNMFEFolder)
                                St = rmdir(CNMFEFolders{P},'s');
                                if ~St
                                    warning(['Some folders could not be deleted, probably because some file is still open.' newline 'Aborting.'])
                                    return
                                end
                            end
                        end
                    end

                    % Probe to see if we can modify the names
                    LastFolder = strsplit(obj.Parameters.PreProcessing.MainFolder,filesep);
                    if isempty(LastFolder{end})
                        LastFolder = LastFolder{end-1};
                    else
                        LastFolder = LastFolder{end};
                    end
                    St = system(['rename ' '"' obj.Parameters.PreProcessing.MainFolder '" '  '"' LastFolder '"'  ]);
                    Basename = obj.Parameters.PreProcessing.Basename;

                    if St~=0
                        ToSave.PreProcessing = obj.Parameters.PreProcessing;
                        ToSave.CNMFE = obj.Parameters.CNMFE;
                        % Export data to a _CalciumData.mat file (the whole structure)
                        save([obj.Parameters.PreProcessing.MainFolder Basename '_CalciumData.mat'],'-struct','ToSave')

                        warning(['Modification of some folders/files is impossible before they are currently open (or locked by Windows otherwise).' newline 'Aborting (a _CalciumData.mat file will be created for convenience, but may lead to some errors until proper export is performed).'])

                        return
                    end

                    % Set all the different step numbers to 1:
                    % Rename filenames/folders and change in metadata files
                    % Preprocessing
                    if exist(obj.Parameters.PreProcessing.PreProcessingFile_isxd,'file')==2
                        NewFN = [obj.Parameters.PreProcessing.Basename '_PP-1.isxd'];
                        St = system(['rename ' '"' obj.Parameters.PreProcessing.PreProcessingFile_isxd '" '  '"' NewFN '"' ]);
                        if St~=0 % Unlikely after the inbitial check but just in case...
                            warning(['Modification of some folders/files is impossible before they are currently open (or locked by Windows otherwise).' newline 'Aborting.'])
                            return
                        end
                    end
                    if exist(obj.Parameters.PreProcessing.PreProcessingFile,'file')==2
                        NewFN = [obj.Parameters.PreProcessing.Basename '_PP-1.h5'];
                        St = system(['rename ' '"' obj.Parameters.PreProcessing.PreProcessingFile '" '  '"' NewFN '"' ]);
                        if St~=0 % Unlikely after the inbitial check but just in case...
                            warning(['Modification of some folders/files is impossible before they are currently open (or locked by Windows otherwise).' newline 'Aborting.'])
                            return
                        end
                    end
                    % Prefiltering
                    if exist(obj.Parameters.PreFiltering.PreFilteringFile,'file')==2
                        NewFN = [obj.Parameters.PreProcessing.Basename '_PP-1_PF-1.h5'];
                        St = system(['rename ' '"' obj.Parameters.PreFiltering.PreFilteringFile '" '  '"' NewFN '"' ]);
                        if St~=0 % Unlikely after the inbitial check but just in case...
                            warning(['Modification of some folders/files is impossible before they are currently open (or locked by Windows otherwise).' newline 'Aborting.'])
                            return
                        end
                    end
                    % Motion correction
                    if exist(obj.Parameters.MotionCorrection.MotionCorrectionFile,'file')==2
                        NewFN = [obj.Parameters.PreProcessing.Basename '_PP-1_PF-1_MC-1.h5'];
                        St = system(['rename ' '"' obj.Parameters.MotionCorrection.MotionCorrectionFile '" '  '"' NewFN '"' ]);
                        if St~=0 % Unlikely after the inbitial check but just in case...
                            warning(['Modification of some folders/files is impossible before they are currently open (or locked by Windows otherwise).' newline 'Aborting.'])
                            return
                        end
                    end

                    % Folders
                    if exist(obj.Parameters.CNMFE.CNMFEFolder,'dir')==7
                        NewFN = ['CNMFE-1'];
                        St = system(['rename ' '"' obj.Parameters.CNMFE.CNMFEFolder '" '  '"' NewFN '"' ]);
                        if St~=0 % Unlikely after the inbitial check but just in case...
                            warning(['Modification of some folders/files is impossible before they are currently open (or locked by Windows otherwise).' newline 'Aborting.'])
                            return
                        end
                    end
                    if exist(obj.Parameters.MotionCorrection.MotionCorrectionFolder,'dir')==7
                        NewFN = ['MC-1'];
                        St = system(['rename ' '"' obj.Parameters.MotionCorrection.MotionCorrectionFolder '" '  '"' NewFN '"' ]);
                        if St~=0 % Unlikely after the inbitial check but just in case...
                            warning(['Modification of some folders/files is impossible before they are currently open (or locked by Windows otherwise).' newline 'Aborting.'])
                            return
                        end
                    end
                    if exist(obj.Parameters.PreFiltering.PreFilteringFolder,'dir')==7
                        NewFN = ['PF-1'];
                        St = system(['rename ' '"' obj.Parameters.PreFiltering.PreFilteringFolder '" '  '"' NewFN '"' ]);
                        if St~=0 % Unlikely after the inbitial check but just in case...
                            warning(['Modification of some folders/files is impossible before they are currently open (or locked by Windows otherwise).' newline 'Aborting.'])
                            return
                        end
                    end
                    if exist(obj.Parameters.PreProcessing.PreProcessingFolder,'dir')==7
                        NewFN = ['PP-1'];
                        St = system(['rename ' '"' obj.Parameters.PreProcessing.PreProcessingFolder '" '  '"' NewFN '"' ]);
                        if St~=0 % Unlikely after the inbitial check but just in case...
                            warning(['Modification of some folders/files is impossible before they are currently open (or locked by Windows otherwise).' newline 'Aborting.'])
                            return
                        end
                    end

                    % Metadata numbering
                    obj.Parameters.PreProcessing.PreProcessingFolder = [obj.Parameters.PreProcessing.MainFolder 'PP-1\'];
                    obj.Parameters.PreProcessing.PreProcessingFile_isxd = [obj.Parameters.PreProcessing.MainFolder 'PP-1\' Basename '_PP-1.isxd'];
                    obj.Parameters.PreProcessing.PreProcessingFile = [obj.Parameters.PreProcessing.MainFolder 'PP-1\' Basename '_PP-1.h5'];
                    obj.Parameters.PreProcessing.ForkNumber = 1;
                    obj.Parameters.PreFiltering.PreFilteringFolder = [obj.Parameters.PreProcessing.MainFolder 'PP-1\PF-1\'];
                    obj.Parameters.PreFiltering.PreFilteringFile = [obj.Parameters.PreProcessing.MainFolder 'PP-1\PF-1\' Basename '_PP-1_PF-1.h5'];
                    obj.Parameters.PreFiltering.ForkNumber = 1;
                    obj.Parameters.MotionCorrection.MotionCorrectionFolder = [obj.Parameters.PreProcessing.MainFolder 'PP-1\PF-1\MC-1\'];
                    obj.Parameters.MotionCorrection.MotionCorrectionFile = [obj.Parameters.PreProcessing.MainFolder 'PP-1\PF-1\MC-1\' Basename '_PP-1_PF-1_MC-1.h5'];
                    obj.Parameters.MotionCorrection.ForkNumber = 1;
                    obj.Parameters.CNMFE.CNMFEFolder = [obj.Parameters.PreProcessing.MainFolder 'PP-1\PF-1\MC-1\CNMFE-1\'];
                    obj.Parameters.CNMFE.DataFile = [obj.Parameters.PreProcessing.MainFolder 'PP-1\PF-1\MC-1\CNMFE-1\' Basename '.h5'];
                    obj.Parameters.CNMFE.ForkNumber = 1;
                    obj.Parameters.CNMFE.Neuron.file = obj.Parameters.CNMFE.DataFile;

                    % Delete all stepmeta files at once
                    StepMetaFiles = dir([obj.Parameters.PreProcessing.MainFolder '**' filesep '*_StepMeta.mat']);
                    arrayfun(@(x) delete([StepMetaFiles(x).folder filesep StepMetaFiles(x).name]),1:numel(StepMetaFiles))

                    ToSave.PreProcessing = obj.Parameters.PreProcessing;
                    save([obj.Parameters.PreProcessing.MainFolder 'PP-1\' Basename '_PP-1_StepMeta.mat'],'-struct','ToSave')
                    ToSave.PreFiltering = obj.Parameters.PreFiltering;
                    save([obj.Parameters.PreProcessing.MainFolder 'PP-1\PF-1\' Basename '_PP-1_PF-1_StepMeta.mat'],'-struct','ToSave')
                    ToSave.MotionCorrection = obj.Parameters.MotionCorrection;
                    save([obj.Parameters.PreProcessing.MainFolder 'PP-1\PF-1\MC-1\' Basename '_PP-1_PF-1_MC-1_StepMeta.mat'],'-struct','ToSave')
                    ToSave.CNMFE = obj.Parameters.CNMFE;
                    save([obj.Parameters.PreProcessing.MainFolder 'PP-1\PF-1\MC-1\CNMFE-1\' Basename '_PP-1_PF-1_MC-1_CNMFE-1_StepMeta.mat'],'-struct','ToSave')

                    % Export data to a _CalciumData.mat file (the whole structure)
                    save([obj.Parameters.PreProcessing.MainFolder Basename '_CalciumData.mat'],'-struct','ToSave')
                end

            end
            obj.EnableAll;

        end


        %% Single elements (reusable beween tabs)
        function SetMovieDisplay1(obj)
            % Movie display axes
            obj.Handles.UIElements.MovieDisplay(1).Plot = axes('Position',[0.195 0.45 0.36 0.45],'Color','k','UserData','1','Parent',obj.Handles.MainFigure); hold on
            obj.Handles.UIElements.MovieDisplay(1).Box = axes('Position',[0.195 0.45 0.36 0.45],'Color','none','UserData','1','Parent',obj.Handles.MainFigure); hold on
            obj.Handles.UIElements.MovieDisplay(1).Plot.Toolbar.Visible = 'off';
            obj.Handles.UIElements.MovieDisplay(1).Plot.Interactions = [];
            obj.Handles.UIElements.MovieDisplay(1).Box.Toolbar.Visible = 'on';

            obj.Handles.UIElements.MovieDisplay(1).Tb = axtoolbar(obj.Handles.UIElements.MovieDisplay(1).Box, {'zoomin', 'zoomout', 'pan'});
            Button = axtoolbarbtn(obj.Handles.UIElements.MovieDisplay(1).Tb, 'push');
            Button.Icon = 'restoreview';
            Button.ButtonPushedFcn = @obj.Restoreview;

            drawnow
            %             linkaxes([obj.Handles.UIElements.MovieDisplay(1).Plot, obj.Handles.UIElements.MovieDisplay(1).Box],'xy');
            obj.Handles.UIElements.MovieDisplay(1).Box.YAxis.Color = [0.8 0.8 0.8];
            obj.Handles.UIElements.MovieDisplay(1).Box.XAxis.Color = [0.8 0.8 0.8];
            obj.Handles.UIElements.MovieDisplay(1).Box.Box = 'on';
            obj.Handles.UIElements.MovieDisplay(1).Box.XTick = [];
            obj.Handles.UIElements.MovieDisplay(1).Box.YTick = [];
            obj.Handles.UIElements.MovieDisplay(1).Plot.TickDir = 'out';
            obj.Handles.UIElements.MovieDisplay(1).Plot.YAxis.Color = [0.8 0.8 0.8];
            obj.Handles.UIElements.MovieDisplay(1).Plot.XAxis.Color = [0.8 0.8 0.8];
            obj.Handles.UIElements.MovieDisplay(1).Plot.Colormap =  bone;
            obj.Handles.UIElements.MovieDisplay(1).Plot.YDir = 'normal';
            obj.Handles.UIElements.MovieDisplay(1).Plot.TickLength = [0 0];
            % We'll need to know the size in px to keep the movie's ratio
            obj.Handles.UIElements.MovieDisplay(1).Plot.Units = 'pixels';
            obj.Handles.UIElements.MovieDisplay(1).AbsolutePosition = obj.Handles.UIElements.MovieDisplay(1).Plot.Position;
            obj.Handles.UIElements.MovieDisplay(1).Plot.Units = 'normalized';
            % Set default limits (just to have something else than [0 1])
            obj.Handles.UIElements.MovieDisplay(1).Plot.XLim = [1 100];
            obj.Handles.UIElements.MovieDisplay(1).Plot.YLim = [1 100];

            obj.Handles.UIElements.BasenameLegend = uicontrol('Style','text',...
                'Units','Normalized','Position',[0.20 0.9 0.36 0.02],...
                'String','',...
                'FontSize',12,'FontName','Arial','FontWeight','b',...
                'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                'HorizontalAlignment','left','Parent',obj.Handles.MainFigure);

            % Player
            obj.Handles.UIElements.PlayButton = uicontrol('Style','pushbutton',...
                'Units','Normalized','Position',[0.195+0.36/2-0.075/2 0.376 0.075 0.04],...
                'String','Play',...
                'FontSize',14,'FontName','Arial','FontWeight','b',...
                'BackgroundColor',[0.15 0.15 0.15],'ForegroundColor',[0.6 0.6 0.6],...
                'HorizontalAlignment','center',...
                'Callback',{@(~,~)obj.Play_CB},'Parent',obj.Handles.MainFigure);
            obj.Handles.UIElements.SlowerButton = uicontrol('Style','pushbutton',...
                'Units','Normalized','Position',[0.195+0.36/2-0.075/2-0.03 0.376 0.03 0.04],...
                'String','< <',...
                'FontSize',14,'FontName','Arial','FontWeight','b',...
                'BackgroundColor',[0.15 0.15 0.15],'ForegroundColor',[0.6 0.6 0.6],...
                'HorizontalAlignment','center',...
                'Callback',{@(~,~)obj.Slower_CB},'Parent',obj.Handles.MainFigure);
            obj.Handles.UIElements.FasterButton = uicontrol('Style','pushbutton',...
                'Units','Normalized','Position',[0.195+0.36/2+0.075/2 0.376 0.03 0.04],...
                'String','> >',...
                'FontSize',14,'FontName','Arial','FontWeight','b',...
                'BackgroundColor',[0.15 0.15 0.15],'ForegroundColor',[0.6 0.6 0.6],...
                'HorizontalAlignment','center',...
                'Callback',{@(~,~)obj.Faster_CB},'Parent',obj.Handles.MainFigure);
            obj.Handles.UIElements.CurrentTimeEdit = uicontrol('Style','edit',...
                'Units','Normalized','Position',[0.195+0.36-0.075-0.02 0.376 0.075 0.04],...
                'String','0',...
                'FontSize',14,'FontName','Arial','FontWeight','b',...
                'BackgroundColor',[0.15 0.15 0.15],'ForegroundColor',[0.6 0.6 0.6],...
                'HorizontalAlignment','center',...
                'Callback',{@(~,~)obj.CurrentTime_CB},'Parent',obj.Handles.MainFigure);
            obj.Handles.UIElements.CurrentTimeLegend = uicontrol('Style','text',...
                'Units','Normalized','Position',[0.195+0.36-0.02 0.366 0.025 0.04],...
                'String',' s',...
                'FontSize',14,'FontName','Arial','FontWeight','b',...
                'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                'HorizontalAlignment','left','Parent',obj.Handles.MainFigure);

            obj.Handles.UIElements.PlayRateLegend = uicontrol('Style','edit',...
                'Units','Normalized','Position',[0.195+(0.36/2-0.075/2-0.03)/2-0.04/2 0.376 0.04 0.04],...
                'String','x1',...
                'FontSize',17,'FontName','Arial','FontWeight','b',...
                'BackgroundColor',[0.15 0.15 0.15],'ForegroundColor',[0.6 0.6 0.6],...
                'HorizontalAlignment','center',...
                'Enable','inactive','Parent',obj.Handles.MainFigure);

            % Slider
            obj.Handles.UIElements.MovieDisplay(1).SliderTickAxes = axes('Unit','Normalized',...
                'Position', [0.195 0.346 0.36 0.02],...
                'Color',[0.2,0.2,0.2],'Parent',obj.Handles.MainFigure); hold on
            obj.Handles.UIElements.MovieDisplay(1).SliderTickAxes.LineWidth = 1;
            obj.Handles.UIElements.MovieDisplay(1).SliderTickAxes.YAxis.Visible =  'off';
            obj.Handles.UIElements.MovieDisplay(1).SliderTickAxes.TickDir = 'out';
            obj.Handles.UIElements.MovieDisplay(1).SliderTickAxes.XColor = 'k';
            obj.Handles.UIElements.MovieDisplay(1).SliderTickAxes.YLim = [0 1];
            obj.Handles.UIElements.MovieDisplay(1).SliderTickAxes.XLim = [0 100];
            obj.Handles.UIElements.MovieDisplay(1).SliderTickAxes.FontSize = 9.5;
            obj.Handles.UIElements.MovieDisplay(1).SliderTickAxes.TickLength(1) = 0.015;

            obj.Handles.UIElements.MovieDisplay(1).SliderAxes = axes('Unit','Normalized',...
                'Position', [0.195 0.345 0.36 0.02],...
                'Color',[0.2,0.2,0.2],'LineWidth',1,'Parent',obj.Handles.MainFigure); hold on
            obj.Handles.UIElements.MovieDisplay(1).TimeLine = plot([0 100],[0 0],...
                'Parent',obj.Handles.UIElements.MovieDisplay(1).SliderAxes,...
                'LineWidth',2.5,'Color','k');
            obj.Handles.UIElements.MovieDisplay(1).Slider = plot([25 25],[-1 1],...
                'Parent',obj.Handles.UIElements.MovieDisplay(1).SliderAxes,...
                'LineWidth',4,'Color','k',...
                'ButtonDownFcn',{@(~,~)obj.SliderCB});
            obj.Handles.UIElements.MovieDisplay(1).SliderAxes.XAxis.Visible = 'off';
            obj.Handles.UIElements.MovieDisplay(1).SliderAxes.YAxis.Visible = 'off';
            if ~verLessThan('matlab','9.5')
                obj.Handles.UIElements.MovieDisplay(1).SliderAxes.Toolbar.Visible = 'off';
                disableDefaultInteractivity(obj.Handles.UIElements.MovieDisplay(1).SliderAxes);
                obj.Handles.UIElements.MovieDisplay(1).SliderTickAxes.Toolbar.Visible = 'off';
                disableDefaultInteractivity(obj.Handles.UIElements.MovieDisplay(1).SliderTickAxes);
            end
            linkaxes([obj.Handles.UIElements.MovieDisplay(1).SliderTickAxes obj.Handles.UIElements.MovieDisplay(1).SliderAxes],'x');
        end

        function SliderCB(obj)
            if ~obj.Dragging,
                if ~isempty(obj.CurrentPlayers),
                    obj.Dragging = true;
                    obj.Handles.MainFigure.WindowButtonMotionFcn = @(~,~)obj.MovingSlider;
                    obj.Handles.MainFigure.WindowButtonUpFcn = @(~,~)obj.SliderCB;
                end
            else
                obj.Dragging = false;
                obj.Handles.MainFigure.WindowButtonMotionFcn = [];
                obj.Handles.MainFigure.WindowButtonUpFcn = [];
            end
        end

        function MovingSlider(obj)
            CurrentCursor = obj.Handles.UIElements.MovieDisplay(1).SliderAxes.CurrentPoint;
            if CurrentCursor(1)>=obj.Handles.UIElements.MovieDisplay(1).SliderAxes.XLim(1) &&  CurrentCursor(1)<=obj.Handles.UIElements.MovieDisplay(1).SliderAxes.XLim(2),
                obj.CurrentTime = CurrentCursor(1);
                CurrentPlayingStatus = obj.Playing;
                obj.Playing = false;
                obj.PlayingSingle = true;
                drawnow
                obj.PlayMovies;
                drawnow
                pause(0.01);
                obj.Playing = CurrentPlayingStatus;
                obj.PlayMovies;
            end
        end

        function SetMovieDisplay2(obj)
            % Movie display axes
            obj.Handles.UIElements.MovieDisplay(2).Plot = axes('Position',[0.605 0.45 0.36 0.45],'Color','k','UserData','2','Parent',obj.Handles.MainFigure); hold on
            obj.Handles.UIElements.MovieDisplay(2).Box = axes('Position',[0.605 0.45 0.36 0.45],'Color','none','UserData','2','Parent',obj.Handles.MainFigure); hold on
            obj.Handles.UIElements.MovieDisplay(2).Tb = axtoolbar(obj.Handles.UIElements.MovieDisplay(2).Box, {'zoomin', 'zoomout', 'pan'});
            obj.Handles.UIElements.MovieDisplay(2).Plot.Toolbar.Visible = 'off';
            obj.Handles.UIElements.MovieDisplay(2).Plot.Interactions = [];
            obj.Handles.UIElements.MovieDisplay(2).Box.Toolbar.Visible = 'on';
            Button = axtoolbarbtn(obj.Handles.UIElements.MovieDisplay(2).Tb, 'push');
            Button.Icon = 'restoreview';
            Button.ButtonPushedFcn = @obj.Restoreview;

            %             linkaxes([obj.Handles.UIElements.MovieDisplay(2).Plot, obj.Handles.UIElements.MovieDisplay(2).Box],'xy');
            obj.Handles.UIElements.MovieDisplay(2).Box.YAxis.Color = [0.8 0.8 0.8];
            obj.Handles.UIElements.MovieDisplay(2).Box.XAxis.Color = [0.8 0.8 0.8];
            obj.Handles.UIElements.MovieDisplay(2).Box.Box = 'on';
            obj.Handles.UIElements.MovieDisplay(2).Box.Color = 'none';
            obj.Handles.UIElements.MovieDisplay(2).Box.XTick = [];
            obj.Handles.UIElements.MovieDisplay(2).Box.YTick = [];
            obj.Handles.UIElements.MovieDisplay(2).Plot.TickDir = 'out';
            obj.Handles.UIElements.MovieDisplay(2).Plot.YAxis.Color = [0.8 0.8 0.8];
            obj.Handles.UIElements.MovieDisplay(2).Plot.XAxis.Color = [0.8 0.8 0.8];
            obj.Handles.UIElements.MovieDisplay(2).Plot.Colormap =  bone;
            obj.Handles.UIElements.MovieDisplay(2).Plot.YDir = 'normal';
            obj.Handles.UIElements.MovieDisplay(2).Plot.TickLength = [0 0];
            % We'll need to know the size in px to keep the movie's ratio
            obj.Handles.UIElements.MovieDisplay(2).Plot.Units = 'pixels';
            obj.Handles.UIElements.MovieDisplay(2).AbsolutePosition = obj.Handles.UIElements.MovieDisplay(2).Plot.Position;
            obj.Handles.UIElements.MovieDisplay(2).Plot.Units = 'normalized';
            % Set default limits (just to have something else than [0 1])
            obj.Handles.UIElements.MovieDisplay(2).Plot.XLim = [1 100];
            obj.Handles.UIElements.MovieDisplay(2).Plot.YLim = [1 100];

            obj.Handles.UIElements.BasenameLegend = uicontrol('Style','text',...
                'Units','Normalized','Position',[0.20 0.9 0.36 0.02],...
                'String','',...
                'FontSize',12,'FontName','Arial','FontWeight','b',...
                'BackgroundColor',[0.2 0.2 0.2],'ForegroundColor',[0.6 0.6 0.6],...
                'HorizontalAlignment','left','Parent',obj.Handles.MainFigure);
        end

        function UpdateDisplay(obj)
            % To update if the movies are not playing
            if ~obj.Playing
                obj.PlayingSingle = true;
                obj.PlayMovies;
            end
        end

        function UpdateKernel_Visu(obj)
            if obj.Parameters.Visualization.Sig1>0 && obj.Parameters.Visualization.Sig2>0
                gaussian1 = fspecial('Gaussian', obj.Parameters.Visualization.WindowSize, obj.Parameters.Visualization.Sig1);
                gaussian2 = fspecial('Gaussian', obj.Parameters.Visualization.WindowSize, obj.Parameters.Visualization.Sig2);
                if obj.Parameters.Visualization.Sig1>obj.Parameters.Visualization.Sig2,
                    psf = gaussian1 - gaussian2;
                else
                    psf = gaussian2 - gaussian1;
                end
            elseif obj.Parameters.Visualization.Sig1 == 0
                psf = fspecial('Gaussian', obj.Parameters.Visualization.WindowSize, obj.Parameters.Visualization.Sig2);
            else
                psf = fspecial('Gaussian', obj.Parameters.Visualization.WindowSize, obj.Parameters.Visualization.Sig1);
            end
            ind_nonzero = (psf(:)>=max(psf(:,1)));
            psf = psf-mean(psf(ind_nonzero));
            psf(~ind_nonzero) = 0;
            obj.Parameters.Visualization.Kernel = psf;
        end

        function Sig1_Visu_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.Sig1_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>0,
                obj.Parameters.Visualization.Sig1 = TempInput;
                obj.UpdateKernel_Visu;
            elseif TempInput == 0
                if obj.Parameters.Visualization.Sig2 == 0,
                    obj.Handles.UIElements.Sig1_Edit.String = num2str(obj.Parameters.Visualization.Sig1);
                else
                    obj.Parameters.Visualization.Sig1 = TempInput;
                    obj.UpdateKernel_Visu;
                end
            else
                obj.Handles.UIElements.Sig1_Edit.String = num2str(obj.Parameters.Visualization.Sig1);
            end
            obj.UpdateDisplay;
        end

        function Sig2_Visu_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.Sig2_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>0,
                obj.Parameters.Visualization.Sig2 = TempInput;
                obj.UpdateKernel_Visu;
            elseif TempInput == 0
                if obj.Parameters.Visualization.Sig1 == 0,
                    obj.Handles.UIElements.Sig2_Edit.String = num2str(obj.Parameters.Visualization.Sig2);
                else
                    obj.Parameters.Visualization.Sig2 = TempInput;
                    obj.UpdateKernel_Visu;
                end
            else
                obj.Handles.UIElements.Sig2_Edit.String = num2str(obj.Parameters.Visualization.Sig2);
            end
            obj.UpdateDisplay;
        end

        function Window_Visu_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.Window_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>=1,
                obj.Parameters.Visualization.WindowSize = TempInput;
                obj.UpdateKernel_Visu;
            else
                obj.Handles.UIElements.Window_Edit.String = num2str(obj.Parameters.Visualization.WindowSize);
            end
            obj.UpdateDisplay;
        end

        function ThresholdLow_Visu_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.ThresholdLow_Visu_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput<obj.Handles.UIElements.MovieDisplay(1).Plot.CLim(2),
                obj.Parameters.Visualization.ThresholdLow = TempInput;
                obj.Handles.UIElements.MovieDisplay(1).Plot.CLim(1) = TempInput;
                if any(cell2mat(strfind(fieldnames(obj.Handles.UIElements.MovieDisplay(2).Plot),'Children'))) && ~isempty(obj.Handles.UIElements.MovieDisplay(2).Plot.Children), % ugly but somehow the fieldname is not "Children"... spaces are there too
                    obj.Handles.UIElements.MovieDisplay(2).Plot.CLim(1) = TempInput;
                end
            else
                obj.Handles.UIElements.ThresholdLow_Visu_Edit.String = num2str( obj.Parameters.Visualization.ThresholdLow);
            end
            obj.UpdateDisplay;
        end

        function ThresholdHigh_Visu_CB(obj)
            TempInput = str2double(obj.Handles.UIElements.ThresholdHigh_Visu_Edit.String);
            if ~isempty(TempInput) && ~isnan(TempInput) && numel(TempInput)==1 && TempInput>obj.Handles.UIElements.MovieDisplay(1).Plot.CLim(1),
                obj.Parameters.Visualization.ThresholdHigh = TempInput;
                obj.Handles.UIElements.MovieDisplay(1).Plot.CLim(2) = TempInput;
                if any(cell2mat(strfind(fieldnames(obj.Handles.UIElements.MovieDisplay(2).Plot),'Children'))) && ~isempty(obj.Handles.UIElements.MovieDisplay(2).Plot.Children), % ugly but somehow the fieldname is not "Children"... spaces are there too
                    obj.Handles.UIElements.MovieDisplay(2).Plot.CLim(2) = TempInput;
                end
            else
                obj.Handles.UIElements.ThresholdHigh_Visu_Edit.String = num2str( obj.Parameters.Visualization.ThresholdHigh);
            end
            obj.UpdateDisplay;
        end

        function PropertiesReset(obj)
            obj.CurrentPlayingStatus = false;
            obj.CurrentTime = 0;
            obj.Dragging = false;
            obj.Playing = false;
            obj.PlayingSingle = false;
        end

        function PostProc_DisplayNumbers_CB(obj)
            if obj.Handles.UIElements.DisplayNumbers.Value
                obj.Handles.Values.PostProcessing.DisplayNumbers = true;
                set(obj.Handles.Numbers_HandlesMovie,'Visible','on')
                set(obj.Handles.Numbers_HandlesMP,'Visible','on')
            else
                obj.Handles.Values.PostProcessing.DisplayNumbers = false;
                set(obj.Handles.Numbers_HandlesMovie,'Visible','off')
                set(obj.Handles.Numbers_HandlesMP,'Visible','off')
            end
        end

        function ZoomPostCallback(obj,~,src)

            src = src.Axes;
            % Check that it is an axes with UserData (e.g. one to take into
            % account here)
            if ~isempty(src.UserData)
                % Get the range for both axis
                XRange = diff(src.XLim);
                YRange = diff(src.YLim);
                XL = src.XLim;
                YL = src.YLim;
                % Keep the largest as reference: extend the other on
                % both sides to keep the ratio
                AbsolutePosition = obj.Handles.UIElements.MovieDisplay(str2double(src.UserData)).AbsolutePosition;
                if XRange >= YRange
                    % YRange will be adjusted
                    XUnitPx = XRange/AbsolutePosition(3);
                    YL_TheoRange = XUnitPx*AbsolutePosition(4);
                    Y_DiffRange = YL_TheoRange - YRange;
                    YL = [YL(1)-0.5*Y_DiffRange YL(2)+0.5*Y_DiffRange];
                else
                    % XRange will be adjusted
                    YUnitPx = YRange/AbsolutePosition(4);
                    XL_TheoRange = YUnitPx*AbsolutePosition(3);
                    X_DiffRange = XL_TheoRange - XRange;
                    XL = [XL(1)-0.5*X_DiffRange XL(2)+0.5*X_DiffRange];
                end
                for MN = 1 : numel(obj.Handles.UIElements.MovieDisplay)
                    obj.Handles.UIElements.MovieDisplay(MN).Plot.XLim = XL;
                    obj.Handles.UIElements.MovieDisplay(MN).Plot.YLim = YL;
                    if isfield(obj.Handles.UIElements.MovieDisplay(MN),'Box')
                        obj.Handles.UIElements.MovieDisplay(MN).Box.YLim = YL;
                        obj.Handles.UIElements.MovieDisplay(MN).Box.XLim = XL;
                    end
                end
            end
            drawnow
        end

        function Restoreview(obj,~,src)
            src = src.Axes;
            % Check that it is an axes with UserData (e.g. one to take into
            % account here)

            if ~isempty(src.UserData)
                for MN = 1 : numel(obj.Handles.UIElements.MovieDisplay)
                    obj.Handles.UIElements.MovieDisplay(MN).Plot.XLim = obj.Handles.UIElements.MovieDisplay(MN).SourceXLim;
                    obj.Handles.UIElements.MovieDisplay(MN).Plot.YLim = obj.Handles.UIElements.MovieDisplay(MN).SourceYLim;
                    if isfield(obj.Handles.UIElements.MovieDisplay(MN),'Box')
                        obj.Handles.UIElements.MovieDisplay(MN).Box.YLim = obj.Handles.UIElements.MovieDisplay(MN).SourceYLim;
                        obj.Handles.UIElements.MovieDisplay(MN).Box.XLim = obj.Handles.UIElements.MovieDisplay(MN).SourceXLim;
                    end
                end
            end
            drawnow
        end
    end

    methods(Static)
        % Function to call outside the GUI to transfer files to the LTS
        function ArchiveCalciumData(varargin)

            % Check optional logical inputs
            % to clean the folders from potential processing tests if
            % there is an exported file
            if nargin>0
                if ~islogical(varargin{1})
                    error('Input must be logical.')
                end
                CleanTrees = varargin{1};
            else
                CleanTrees = false;
            end
            % to export/keep only the raw data and the CNMFE
            % metadatafile
            if nargin>1
                if ~islogical(varargin{2})
                    error('Input must be logical.')
                end
                MinimalTransfer = varargin{2};
            else
                MinimalTransfer = false;
            end

            % Prompt to choose a folder to transfer
            Experimenters = DataBase.Lists.GetList('Experimenters');
            CurrentExperimenter = Experimenters(strcmpi(Experimenters(:,2),getenv('username')),1);
            Path = uigetdir(['F:\' CurrentExperimenter{1} '\Data\CalciumImaging\'],'Please choose a folder...');
            if isempty(Path) | Path == 0
                return
            end

            % Look for Extracted files and check whether corresponding data
            % was transferred
            Files = dir([Path '\**\' filesep '*_CalciumData.mat']);
            % People might have given different names to the drive, or not
            % added it at all: use the explicit name
            ArchiveDrive = '\\VNBANSRV2\LTS-Ablage\';
            if isempty(Files)
                disp('Could not find any exported session to check in the folder that was chosen.')
            else
                FilesList = arrayfun(@(x) [Files(x).folder filesep Files(x).name],1:numel(Files),'UniformOutput',false);
                for F = 1 : numel(FilesList)
                    % Load the file
                    Loaded = load(FilesList{F});
                    disp([repmat('_',[1,75]),newline,...
                        'Checking session ' num2str(F) ' out of ' num2str(numel(Files)) ':' newline,...
                        Loaded.PreProcessing.Basename, newline])

                    % Clean folders if enabled
                    if CleanTrees
                        disp(['Cleaning the arborescence...'])
                        % Delete all the other branches, step by step
                        % Try to close all files
                        fclose('all');
                        % Preprocessing
                        PPFolders = dir([Loaded.PreProcessing.MainFolder 'PP-*' ]);
                        PPFolders = PPFolders([PPFolders.isdir]);
                        PPFolders = arrayfun(@(x) [PPFolders(x).folder filesep PPFolders(x).name filesep],1:numel(PPFolders),'UniformOutput',false);
                        if isempty(PPFolders)
                            error('The folder structure has been manually edited and corrupted. Aborting.')
                        else
                            for P = 1 : numel(PPFolders)
                                if ~strcmpi(PPFolders{P},Loaded.PreProcessing.PreProcessingFolder)
                                    St = rmdir(PPFolders{P},'s');
                                    if ~St
                                        warning(['Some folders could not be deleted, probably because some file is still open.' newline 'Aborting.'])
                                        return
                                    end
                                end
                            end
                        end

                        % Prefiltering
                        PFFolders = dir([Loaded.PreProcessing.PreProcessingFolder 'PF-*' ]);
                        PFFolders = PFFolders([PFFolders.isdir]);
                        PFFolders = arrayfun(@(x) [PFFolders(x).folder filesep PFFolders(x).name filesep],1:numel(PFFolders),'UniformOutput',false);
                        if isempty(PFFolders)
                            error('The folder structure has been manually edited and corrupted. Aborting.')
                        else
                            for P = 1 : numel(PFFolders)
                                if ~strcmpi(PFFolders{P},Loaded.PreFiltering.PreFilteringFolder)
                                    St = rmdir(PFFolders{P},'s');
                                    if ~St
                                        warning(['Some folders could not be deleted, probably because some file is still open.' newline 'Aborting.'])
                                        return
                                    end
                                end
                            end
                        end

                        % Motion correction
                        MCFolders = dir([Loaded.PreFiltering.PreFilteringFolder 'MC-*' ]);
                        MCFolders = MCFolders([MCFolders.isdir]);
                        MCFolders = arrayfun(@(x) [MCFolders(x).folder filesep MCFolders(x).name filesep],1:numel(MCFolders),'UniformOutput',false);
                        if isempty(MCFolders)
                            error('The folder structure has been manually edited and corrupted. Aborting.')
                        else
                            for P = 1 : numel(MCFolders)
                                if ~strcmpi(MCFolders{P},Loaded.MotionCorrection.MotionCorrectionFolder)
                                    St = rmdir(MCFolders{P},'s');
                                    if ~St
                                        warning(['Some folders could not be deleted, probably because some file is still open.' newline 'Aborting.'])
                                        return
                                    end
                                end
                            end
                        end

                        % CNMFE
                        CNMFEFolders = dir([Loaded.MotionCorrection.MotionCorrectionFolder 'CNMFE-*' ]);
                        CNMFEFolders = CNMFEFolders([CNMFEFolders.isdir]);
                        CNMFEFolders = arrayfun(@(x) [CNMFEFolders(x).folder filesep CNMFEFolders(x).name filesep],1:numel(CNMFEFolders),'UniformOutput',false);
                        if isempty(CNMFEFolders)
                            error('The folder structure has been manually edited and corrupted. Aborting.')
                        else
                            for P = 1 : numel(CNMFEFolders)
                                if ~strcmpi(CNMFEFolders{P},Loaded.CNMFE.CNMFEFolder)
                                    St = rmdir(CNMFEFolders{P},'s');
                                    if ~St
                                        warning(['Some folders could not be deleted, probably because some file is still open.' newline 'Aborting.'])
                                        return
                                    end
                                end
                            end
                        end
                        disp(['Done!' newline])
                    end



                    % Check that we have transferred the necessary files
                    ArchivePath = CalciumImaging_Pipeline.ToArchivePath(Loaded.PreProcessing.MainFolder);
                    RawFile = Loaded.PreProcessing.RawFile;
                    ArchivedRawFile = CalciumImaging_Pipeline.ToArchivePath(RawFile);
                    ArchivedExported = [ArchivePath Loaded.PreProcessing.Basename '_CalciumData.mat'];
                    disp(['Checking that we have the complete folder structure and all the files in the archive...' newline])

                    % Folder
                    if exist(ArchivePath,'dir')~=7
                        disp(['Creating folder for the archive:' newline,...
                            ArchivePath newline])
                        mkdir(ArchivePath)
                    end

                    % Raw file and extracted data
                    if exist(ArchivedExported,'file')~=2
                        disp('Transferring final calcium data file...')
                        copyfile(FilesList{F},ArchivedExported)
                        disp('Done!')
                    else
                        TempLoad = load(ArchivedExported);
                        if isfield(TempLoad.CNMFE,'DateExported') && isfield(Loaded.CNMFE,'DateExported')
                            if TempLoad.CNMFE.DateExported ~= Loaded.CNMFE.DateExported
                                disp('Transferring final calcium data file...')
                                copyfile(FilesList{F},ArchivedExported)
                                disp('Done!')
                            end
                        else
                            disp('Transferring final calcium data file...')
                            copyfile(FilesList{F},ArchivedExported)
                            disp('Done!')
                        end
                    end
                    if exist(ArchivedRawFile,'file')~=2
                        disp('Transferring raw .isxd file...')
                        copyfile(RawFile,ArchivedRawFile)
                        disp('Done!')
                    end

                    % Last CNMFE file (same file, but to keep the
                    % arborescence; also create folders)
                    ArchivedCNMFEFolder = CalciumImaging_Pipeline.ToArchivePath(Loaded.CNMFE.CNMFEFolder);
                    ArchivedCNMFEFile = [ArchivedCNMFEFolder Loaded.PreProcessing.Basename '_PP-1_PF-1_MC-1_CNMFE-1_StepMeta.mat'];
                    if exist(ArchivedCNMFEFolder,'dir')~=7
                        mkdir(ArchivedCNMFEFolder)
                    end
                    if exist(ArchivedCNMFEFile,'file')==2
                        TempLoad = load(ArchivedCNMFEFile);
                        if isfield(TempLoad.CNMFE,'DateExported') && isfield(Loaded.CNMFE,'DateExported')
                            if TempLoad.CNMFE.DateExported ~= Loaded.CNMFE.DateExported
                                disp('Transferring final CNMFE metafile...')
                                copyfile(FilesList{F},ArchivedCNMFEFile)
                                disp('Done!')
                            end
                        else
                            disp('Transferring final CNMFE metafile...')
                            copyfile(FilesList{F},ArchivedCNMFEFile)
                            disp('Done!')
                        end
                    else
                        disp('Transferring final CNMFE metafile...')
                        copyfile(FilesList{F},ArchivedCNMFEFile)
                        disp('Done!')
                    end
                    Basename = Loaded.PreProcessing.Basename;
                    % Other folders/files if appropriate
                    if ~MinimalTransfer
                        % Preprocessing
                        PPFile = [Loaded.PreProcessing.PreProcessingFolder Loaded.PreProcessing.Basename '_PP-1_StepMeta.mat'];
                        ArchivedPPFile = CalciumImaging_Pipeline.ToArchivePath(PPFile);
                        if exist(PPFile,'file')==2
                            if exist(ArchivedPPFile,'file')==2
                                TempLoad = load(ArchivedPPFile);
                                if isfield(TempLoad.PreProcessing,'Date') && isfield(Loaded.PreProcessing,'Date')
                                    if TempLoad.PreProcessing.Date ~= Loaded.PreProcessing.Date
                                        disp('Transferring final preprocessing metafile...')
                                        copyfile(PPFile,ArchivedPPFile)
                                        disp('Done!')
                                    end
                                else
                                    disp('Transferring final preprocessing metafile...')
                                    copyfile(PPFile,ArchivedPPFile)
                                    disp('Done!')
                                end
                            else
                                disp('Transferring final preprocessing metafile...')
                                copyfile(PPFile,ArchivedPPFile)
                                disp('Done!')
                            end
                            ArchivedMetaPPFile = ArchivedPPFile;
                            PPFile = Loaded.PreProcessing.PreProcessingFile;
                            ArchivedPPFile = CalciumImaging_Pipeline.ToArchivePath(PPFile);
                            if exist(ArchivedPPFile,'file')==2
                                TempLoad = load(ArchivedMetaPPFile);
                                if isfield(TempLoad.PreProcessing,'Date') && isfield(Loaded.PreProcessing,'Date')
                                    if TempLoad.PreProcessing.Date ~= Loaded.PreProcessing.Date
                                        disp('Transferring final preprocessing .h5 file...')
                                        copyfile(PPFile,ArchivedPPFile)
                                        disp('Done!')
                                    end
                                else
                                    disp('Transferring final preprocessing .h5 file...')
                                    copyfile(PPFile,ArchivedPPFile)
                                    disp('Done!')
                                end
                            else
                                disp('Transferring final preprocessing .h5 file...')
                                copyfile(PPFile,ArchivedPPFile)
                                disp('Done!')
                            end
                        end
                        PPFile = Loaded.PreProcessing.PreProcessingFile_isxd;
                        ArchivedPPFile = CalciumImaging_Pipeline.ToArchivePath(PPFile);
                        if exist(PPFile,'file')==2
                            if exist(ArchivedPPFile,'file')==2
                                TempLoad = load(ArchivedMetaPPFile);
                                if isfield(TempLoad.PreProcessing,'Date') && isfield(Loaded.PreProcessing,'Date')
                                    if TempLoad.PreProcessing.Date ~= Loaded.PreProcessing.Date
                                        disp('Transferring final preprocessing .isxd file...')
                                        copyfile(PPFile,ArchivedPPFile)
                                        disp('Done!')
                                    end
                                else
                                    disp('Transferring final preprocessing .isxd file...')
                                    copyfile(PPFile,ArchivedPPFile)
                                    disp('Done!')
                                end
                            else
                                disp('Transferring final preprocessing .isxd file...')
                                copyfile(PPFile,ArchivedPPFile)
                                disp('Done!')
                            end
                        end

                        % Prefiltering
                        PFFile = [Loaded.PreFiltering.PreFilteringFolder Basename '_PP-1_PF-1_StepMeta.mat'];
                        if exist(PFFile,'file')==2
                            ArchivedPFFile = CalciumImaging_Pipeline.ToArchivePath(PFFile);
                            if exist(ArchivedPFFile,'file')==2
                                TempLoad = load(ArchivedPFFile);
                                if isfield(TempLoad.PreFiltering,'Date') && isfield(Loaded.PreFiltering,'Date')
                                    if TempLoad.PreFiltering.Date ~= Loaded.PreFiltering.Date
                                        disp('Transferring final PreFiltering metafile...')
                                        copyfile(PFFile,ArchivedPFFile)
                                        disp('Done!')
                                    end
                                else
                                    disp('Transferring final PreFiltering metafile...')
                                    copyfile(PFFile,ArchivedPFFile)
                                    disp('Done!')
                                end
                            else
                                disp('Transferring final PreFiltering metafile...')
                                copyfile(PFFile,ArchivedPFFile)
                                disp('Done!')
                            end
                        end
                        ArchivedMetaPFFile = ArchivedPFFile;
                        PFFile = Loaded.PreFiltering.PreFilteringFile;
                        if exist(PFFile,'file')==2
                            ArchivedPFFile = CalciumImaging_Pipeline.ToArchivePath(PFFile);
                            if exist(ArchivedPFFile,'file')==2
                                TempLoad = load(ArchivedMetaPFFile);
                                if isfield(TempLoad.PreFiltering,'Date') && isfield(Loaded.PreFiltering,'Date')
                                    if TempLoad.PreFiltering.Date ~= Loaded.PreFiltering.Date
                                        disp('Transferring final PreFiltering .h5 file...')
                                        copyfile(PFFile,ArchivedPFFile)
                                        disp('Done!')
                                    end
                                else
                                    disp('Transferring final PreFiltering .h5 file...')
                                    copyfile(PFFile,ArchivedPFFile)
                                    disp('Done!')
                                end
                            else
                                disp('Transferring final PreFiltering .h5 file...')
                                copyfile(PFFile,ArchivedPFFile)
                                disp('Done!')
                            end
                        end
                        % Motion correction
                        MCFile = [Loaded.MotionCorrection.MotionCorrectionFolder Basename '_PP-1_PF-1_MC-1_StepMeta.mat'];
                        if exist(MCFile,'file')==2
                            ArchivedMCFile = CalciumImaging_Pipeline.ToArchivePath(MCFile);
                            if exist(ArchivedMCFile,'file')==2
                                TempLoad = load(ArchivedMCFile);
                                if isfield(TempLoad.MotionCorrection,'Date') && isfield(Loaded.MotionCorrection,'Date')
                                    if TempLoad.MotionCorrection.Date ~= Loaded.MotionCorrection.Date
                                        disp('Transferring final MotionCorrection metafile...')
                                        copyfile(MCFile,ArchivedMCFile)
                                        disp('Done!')
                                    end
                                else
                                    disp('Transferring final MotionCorrection metafile...')
                                    copyfile(MCFile,ArchivedMCFile)
                                    disp('Done!')
                                end
                            else
                                disp('Transferring final MotionCorrection metafile...')
                                copyfile(MCFile,ArchivedMCFile)
                                disp('Done!')
                            end
                        end
                        ArchivedMetaMCFile = ArchivedMCFile;
                        MCFile = Loaded.MotionCorrection.MotionCorrectionFile;
                        if exist(MCFile,'file')==2
                            ArchivedMCFile = CalciumImaging_Pipeline.ToArchivePath(MCFile);
                            if exist(ArchivedMCFile,'file')==2
                                TempLoad = load(ArchivedMetaMCFile);
                                if isfield(TempLoad.MotionCorrection,'Date') && isfield(Loaded.MotionCorrection,'Date')
                                    if TempLoad.MotionCorrection.Date ~= Loaded.MotionCorrection.Date
                                        disp('Transferring final MotionCorrection .h5 file...')
                                        copyfile(MCFile,ArchivedMCFile)
                                        disp('Done!')
                                    end
                                else
                                    disp('Transferring final MotionCorrection .h5 file...')
                                    copyfile(MCFile,ArchivedMCFile)
                                    disp('Done!')
                                end
                            else
                                disp('Transferring final MotionCorrection .h5 file...')
                                copyfile(MCFile,ArchivedMCFile)
                                disp('Done!')
                            end
                        end

                        % CNMFE files
                        CNMFEFile = [Loaded.CNMFE.CNMFEFolder Basename '_PP-1_PF-1_CNMFE-1_StepMeta.mat'];
                        ArchivedMetaCNMFEFile = CalciumImaging_Pipeline.ToArchivePath(CNMFEFile);
                        CNMFEFile = Loaded.CNMFE.DataFile;
                        if exist(CNMFEFile,'file')==2
                            ArchivedCNMFEFile = CalciumImaging_Pipeline.ToArchivePath(CNMFEFile);
                            if exist(ArchivedCNMFEFile,'file')==2 && exist(ArchivedMetaCNMFEFile,'file')==2
                                TempLoad = load(ArchivedMetaCNMFEFile);
                                if isfield(TempLoad.CNMFE,'Date') && isfield(Loaded.CNMFE,'Date')
                                    if TempLoad.CNMFE.Date ~= Loaded.CNMFE.Date
                                        disp('Transferring final CNMFE .h5 file...')
                                        copyfile(CNMFEFile,ArchivedCNMFEFile)
                                        disp('Done!')
                                    end
                                else
                                    disp('Transferring final CNMFE .h5 file...')
                                    copyfile(CNMFEFile,ArchivedCNMFEFile)
                                    disp('Done!')
                                end
                            else
                                disp('Transferring final CNMFE .h5 file...')
                                copyfile(CNMFEFile,ArchivedCNMFEFile)
                                disp('Done!')
                            end
                        end
                        CNMFEFolder = [Loaded.CNMFE.CNMFEFolder Basename '_source_extraction\'];
                        if exist(CNMFEFolder,'dir')==7
                            ArchivedCNMFEFolder = CalciumImaging_Pipeline.ToArchivePath(CNMFEFolder);
                            if exist(ArchivedCNMFEFile,'dir')==7
                                TempLoad = load(ArchivedMetaCNMFEFile);
                                if isfield(TempLoad.CNMFE,'Date') && isfield(Loaded.CNMFE,'Date')
                                    if TempLoad.CNMFE.Date ~= Loaded.CNMFE.Date
                                        St = rmdir(ArchivedCNMFEFolder);
                                        if ~St
                                            warning(['Some folders could not be deleted, probably because some file is still open.' newline 'Aborting.'])
                                            return
                                        end
                                        mkdir(ArchivedCNMFEFolder)
                                        disp('Transferring final CNMFE data files...')
                                        St = system(['xcopy /e /v /d "' CNMFEFolder(1:end-1) '" "' ArchivedCNMFEFolder '"']);
                                        if St~=0
                                            warning(['Some files could not be copied.' newline 'Aborting.'])
                                            return
                                        end
                                        disp('Done!')
                                    end
                                else
                                    St = rmdir(ArchivedCNMFEFolder);
                                    if ~St
                                        warning(['Some folders could not be deleted, probably because some file is still open.' newline 'Aborting.'])
                                        return
                                    end
                                    mkdir(ArchivedCNMFEFolder)
                                    disp('Transferring final CNMFE data files...')
                                    St = system(['xcopy /e /v /d "' CNMFEFolder(1:end-1) '" "' ArchivedCNMFEFolder '"']);
                                    if St~=0
                                        warning(['Some files could not be copied.' newline 'Aborting.'])
                                        return
                                    end
                                    disp('Done!')
                                end
                            else
                                mkdir(ArchivedCNMFEFolder)
                                disp('Transferring final CNMFE data file...')
                                St = system(['xcopy /e /v /d "' CNMFEFolder(1:end-1) '" "' ArchivedCNMFEFolder '"']);
                                if St~=0
                                    warning(['Some files could not be copied.' newline 'Aborting.'])
                                    return
                                end
                                disp('Done!')
                            end
                        end
                    end
                    disp(['Done!' newline])
                    % Check what we can delete; we're going to check very
                    % carefully only the raw data file and the final file
                    disp(['Checking if some files need to be deleted...'])
                    if MinimalTransfer
                        % CNMFE source extraction folder
                        CNMFEFolder = [Loaded.CNMFE.CNMFEFolder Basename '_source_extraction\'];
                        if exist(CNMFEFolder,'dir')==7
                            St = rmdir(CNMFEFolder);
                            if ~St
                                warning(['Some folders could not be deleted, probably because some file is still open.' newline 'Aborting.'])
                                return
                            end
                        end
                        CNMFEFolder = CalciumImaging_Pipeline.ToArchivePath(CNMFEFolder);
                        if exist(CNMFEFolder,'dir')==7
                            St = rmdir(CNMFEFolder);
                            if ~St
                                warning(['Some folders could not be deleted, probably because some file is still open.' newline 'Aborting.'])
                                return
                            end
                        end
                        % CNMFE file
                        CNMFEFile = Loaded.CNMFE.DataFile;
                        CNMFEFileArchive = CalciumImaging_Pipeline.ToArchivePath(CNMFEFile);
                        if exist(CNMFEFile,'file')==2
                            delete(CNMFEFile);
                        end
                        if exist(CNMFEFileArchive,'file')==2
                            delete(CNMFEFileArchive);
                        end
                        % Prefiltering file
                        PFFile = Loaded.PreFiltering.PreFilteringFile;
                        PFFileArchive = CalciumImaging_Pipeline.ToArchivePath(PFFile);
                        if exist(PFFile,'file')==2
                            delete(PFFile);
                        end
                        if exist(PFFileArchive,'file')==2
                            delete(PFFileArchive);
                        end
                        % Preprocessing files
                        PPFile = Loaded.PreProcessing.PreProcessingFile;
                        PPFileArchive = CalciumImaging_Pipeline.ToArchivePath(PPFile);
                        if exist(PPFile,'file')==2
                            delete(PPFile);
                        end
                        if exist(PPFileArchive,'file')==2
                            delete(PPFileArchive);
                        end
                        PPFile = Loaded.PreProcessing.PreProcessingFile_isxd;
                        PPFileArchive = CalciumImaging_Pipeline.ToArchivePath(PPFile);
                        if exist(PPFile,'file')==2
                            delete(PPFile);
                        end
                        if exist(PPFileArchive,'file')==2
                            delete(PPFileArchive);
                        end
                        % Raw file
                        RawFile = Loaded.PreProcessing.RawFile;
                        DirRaw = dir(RawFile);
                        RawFileArchive = CalciumImaging_Pipeline.ToArchivePath(RawFile);
                        DirArchive = dir(RawFileArchive);
                        if exist(RawFile,'file')==2 && exist(RawFileArchive,'file')==2
                            if (DirRaw.bytes == DirArchive.bytes) && strcmpi(DirRaw.date,DirArchive.date)
                                delete(RawFile);
                            else
                                warning('Did not delete the raw file from F: drive because it does not match the archived file.')
                            end
                        end
                    else
                        % CNMFE source extraction folder
                        CNMFEFolder = [Loaded.CNMFE.CNMFEFolder Basename '_source_extraction\'];
                        if exist(CNMFEFolder,'dir')==7
                            St = rmdir(CNMFEFolder);
                            if ~St
                                warning(['Some folders could not be deleted, probably because some file is still open.'])
                            end
                        end
                        % CNMFE file
                        CNMFEFile = Loaded.CNMFE.DataFile;
                        if exist(CNMFEFile,'file')==2
                            delete(CNMFEFile);
                        end
                        % Prefiltering file
                        PFFile = Loaded.PreFiltering.PreFilteringFile;
                        if exist(PFFile,'file')==2
                            delete(PFFile);
                        end
                        % Preprocessing files
                        PPFile = Loaded.PreProcessing.PreProcessingFile;
                        if exist(PPFile,'file')==2
                            delete(PPFile);
                        end
                        PPFile = Loaded.PreProcessing.PreProcessingFile_isxd;
                        if exist(PPFile,'file')==2
                            delete(PPFile);
                        end
                        % Raw file
                        RawFile = Loaded.PreProcessing.RawFile;
                        DirRaw = dir(RawFile);
                        RawFileArchive = CalciumImaging_Pipeline.ToArchivePath(RawFile);
                        DirArchive = dir(RawFileArchive);
                        if exist(RawFile,'file')==2 && exist(RawFileArchive,'file')==2
                            if (DirRaw.bytes == DirArchive.bytes) && strcmpi(DirRaw.date,DirArchive.date)
                                delete(RawFile);
                            else
                                warning('Did not delete the raw file from F: drive because it does not match the archived file.')
                            end
                        end
                    end
                    disp(['Done!' newline])
                end
            end
        end

        function ArchivePath = ToArchivePath(InPath)
            ArchiveDrive = '\\VNBANSRV2\LTS-Ablage\';
            ArchivePath = [ArchiveDrive InPath(4:end)];
        end
    end
end