classdef Explore_Longitudinal_Alignment <handle

    properties(Access=protected)
        % Metadata
        CellReg_File
        CNMFE_Files
        Mouse
        Paradigms

        % Data
        Values
        cellRegistered
        Contours
        FootPrints
        Traces
        TracesRaw
        Times

        % CellReg parameters
        pixel_weight_threshold=0.15; % for better visualization of cells

        % Contour extraction
        Contour_Threshold = 3;

        % Plotting
        NonSelected_Contour_Width = 1;
        NonSelected_Contour_Color = [0.6 0.6 0.6];
        Selected_Contour_Width = 2.5;
        Neutral_Contour_Color = [0.929 0.6940 0.125];
        Bad_Contour_Color = [0.85 0.325 0.0980];
        Good_Contour_Color = [0 0.447 0.741];
        AllColors

        AxesPosition
        Chan_Colors = {'R','G','B'};
        Figure
        FontScaling = 1;
        Handles
        Selected_Index
        Spacing_Horizontal = 0.02;
        Spacing_Vertical = 0.05;
        SubPlots_Ratio = 0.8;

        % Interactivity
        PlotMode = 'Paradigms';
        CenterOnROI = true;
        Color_Coding = 'Alignment';
        CrossHair_Coordinates = [NaN NaN];
        DisplayContours = false;
        DisplayNumbers = false;
        DisplayTraces = false;
        EnableCrossHair = true;
        KeepTrack = false;
        LinkTraces = false;
        MissingColor = 'y';
        PauseRefreshing = false;
        Smoothing = 1;
        SingleCellMode = 'FootPrints';
        Temp
        Zoom = false;

        % Data selection
        CellCategories
        SelectedParadigms
        SelectedSingleCell
        SelectedCells
        SelectedROIs
    end

    methods
        function obj = Explore_Longitudinal_Alignment
            % Choose a CellReg output file
            Experimenters = DataBase.Lists.GetList('Experimenters');
            IndxExpe = strcmpi(Experimenters(:,2),getenv('username'));
            if any(IndxExpe)
                CurrentExperimenter = Experimenters(IndxExpe,1);
                [Path] = ['G:\' CurrentExperimenter{1} '\Data\'];
            else
                [Path] ='G:\';
            end

            [File,Path] = uigetfile([Path '*_cellRegistered.mat'],'Please select a file to explore.');
            if File==0
                warning('No file selected. Aborting.')
                return
            end
            obj.CellReg_File = fullfile(Path,File);

            % Load alignment data
            obj.cellRegistered = load(obj.CellReg_File);
            if isfield(obj.cellRegistered,'LoadedGlobal')
                % We have loaded a global alignment file
                if numel(obj.cellRegistered.LoadedGlobal)>1
                    % We have several planes and need to choose one
                    [IndexPlane] = listdlg('PromptString','Select the plane to display:',...
                        'SelectionMode','single',...
                        'ListString',num2str((1:numel(obj.cellRegistered.LoadedGlobal))'));
                    if isempty(IndexPlane)
                        return
                    end
                else
                    IndexPlane = 1;
                end
                obj.cellRegistered = obj.cellRegistered.LoadedGlobal(IndexPlane);
                if numel(obj.cellRegistered.Alignments)>1
                    % We have several alignments and need to choose one
                    AlignmentsParadigmsList = arrayfun(@(x) [strjoin(obj.cellRegistered.Alignments(x).FullParadigms(:),',')  ' (ref.: ' obj.cellRegistered.Alignments(x).Reference ')'],1:numel(obj.cellRegistered.Alignments),'UniformOutput',false);
                    [IndexAlignment] = listdlg('PromptString','Select the alignment to display:',...
                        'SelectionMode','single',...
                        'ListString',AlignmentsParadigmsList,'ListSize',[500 500]);
                    if isempty(IndexAlignment)
                        return
                    end
                else
                    IndexAlignment = 1;
                end
                obj.cellRegistered = obj.cellRegistered.Alignments(IndexAlignment);
            else
                obj.cellRegistered = obj.cellRegistered.cell_registered_struct;
            end

            % Retrieve basic infos
            if isfield(obj.cellRegistered,'FullParadigms')
                obj.Paradigms = obj.cellRegistered.FullParadigms;
            else
                obj.Paradigms = obj.cellRegistered.Paradigms;
            end
            obj.Mouse = obj.cellRegistered.Mouse;
            obj.CNMFE_Files = obj.cellRegistered.CNMFE_Files;
            Categories = cell(numel(obj.Paradigms),1);
            % Load CNMFE data (traces)
            for P = 1 : numel(obj.Paradigms)
                TempLoaded = load(obj.CNMFE_Files{P});
                % Get CNMFE folder
                TempFolder = fileparts(obj.CNMFE_Files{P});
                LoadedCNMFE = load(obj.CNMFE_Files{P});
                CatFile = [TempFolder filesep LoadedCNMFE.PreProcessing.Basename '_CellCategories.mat'];
                if exist(CatFile,'file')==2
                    CatTemp = load(CatFile);
                    Categories{P} = CatTemp.Categories;
                end
                obj.TracesRaw{P} = TempLoaded.CNMFE.Neuron.C_raw;
                obj.Times{P} = TempLoaded.PreProcessing.TimeStamps';
            end
            obj.GetTraces; % Smooth

            % Generate figure
            % Prepare the GUI
            set(0,'Units','pixels')
            Scrsz = get(0,'ScreenSize');

            obj.Figure = figure('Position',[0 45 Scrsz(3) Scrsz(4)-75],'MenuBar','none','ToolBar','figure','Renderer','painters');
            addToolbarExplorationButtons(obj.Figure)
            obj.Figure.WindowState = 'maximized'; % deals with windows taskbar so better than guessing the position
            obj.Figure.WindowKeyPressFcn = @(src,evt)obj.PressKeyCB(src,evt); % Callback to catch key presses to accomplish different actions:
            % RIGHTARROW will go to the next cell
            % LEFTARROW will mark the current cell as "good"
            % DOWNARROW will mark the current cell as "bad"

            % Font scaling for different screen resolutions
            obj.FontScaling = min([Scrsz(3)/1920,Scrsz(4)/1080]);

            LeftSpace = 0.2;
            [p,n] = numSubplots(numel(obj.Paradigms));
            numSub = numel(obj.Paradigms);
            ax_width = (1-LeftSpace - (p(2)+1)*obj.Spacing_Horizontal)/p(2);
            ax_height = (1-(p(1)+1)*obj.Spacing_Vertical)/p(1);
            Spacing = obj.Spacing_Vertical + (1-obj.SubPlots_Ratio) * ax_height;


            cells_in_all_days = sum(obj.cellRegistered.cell_to_index_map'>0)==numel(obj.Paradigms);
            obj.Selected_Index = cells_in_all_days;

            for Pr = 1 : numel(obj.Paradigms)
                % Calculate the row and column index of the current axis
                row = ceil(Pr / p(2));
                col = mod(Pr-1, p(2)) + 1;

                % Calculate the position of the current axis
                left = LeftSpace+ (col - 1) * (ax_width + obj.Spacing_Horizontal);
                bottom = Spacing + (p(1) - row) * (obj.SubPlots_Ratio*ax_height + Spacing);

                obj.Handles.SubROI(Pr).BackGroundAxis = axes('Units','normalized','Position',...
                    [left,bottom,ax_width,obj.SubPlots_Ratio * ax_height],'Color','k','XColor','none','YColor','none');
                obj.Handles.SubROI(Pr).Axis = axes('Units','normalized','Position',...
                    [left,bottom,ax_width,obj.SubPlots_Ratio * ax_height]);
                obj.AxesPosition{Pr}.Split.SubROI = obj.Handles.SubROI(Pr).Axis.Position;

                PC = permute(obj.cellRegistered.spatial_footprints_corrected{Pr},[2 3 1]);
                % Apply a cutoff per ROI
                for R = 1 : size(PC,3)
                    SubR = PC(:,:,R);
                    Max = max(SubR,[],'all');
                    IndxZ = SubR<(Max*obj.pixel_weight_threshold);
                    SubR(IndxZ) = 0;
                    PC(:,:,R) = SubR;
                end
                % Normalize intensity overall
                NZ = PC(PC~=0);
%                 MaxAll = prctile(NZ,99.5,'all');
%                 MinAll = prctile(NZ,20,'all');
                MinAll = min(NZ,[],'all');
                MaxAll = max(NZ,[],'all');
                PC(PC<MinAll) = MinAll;
%                 PC(PC>MaxAll) = MaxAll;
                PC = (PC - MinAll)/ (MaxAll - MinAll);
%                 figure; imagesc(sum(PC,3))

                obj.FootPrints(Pr).FootPrints = PC;
                obj.FootPrints(Pr).Base = repmat(sum(PC,3),[1,1,3]);

                obj.Handles.SubROI(Pr).FootPrints = [];
                obj.Handles.SubROI(Pr).Contours = [];

                % Derive the contours (CNMFE method)
                for Fp = 1 : size(obj.FootPrints(Pr).FootPrints,3)
                    obj.Contours{Pr}{Fp} = obj.Get_Contours(obj.FootPrints(Pr).FootPrints(:,:,Fp),obj.Contour_Threshold);
                end

                % Plot trace
                bottom = obj.Spacing_Vertical + (p(1) - row) * (ax_height + obj.Spacing_Vertical);
                obj.Handles.SubTrace(Pr).Axis = axes('Units','normalized','Position',...
                    [left,bottom,ax_width,(1-obj.SubPlots_Ratio) * ax_height]);
                obj.AxesPosition{Pr}.Split.SubTrace = obj.Handles.SubTrace(Pr).Axis.Position;
                obj.AxesPosition{Pr}.Full = obj.Handles.SubTrace(Pr).Axis.Position;
                obj.AxesPosition{Pr}.Full(4) = obj.AxesPosition{Pr}.Full(4) + obj.Handles.SubROI(Pr).Axis.Position(4);
                obj.Handles.SubTrace(Pr).Trace = plot([NaN NaN],[NaN NaN],'k','LineWidth',0.5,'Visible','on','Parent',obj.Handles.SubTrace(Pr).Axis);
                obj.Handles.SubTrace(Pr).Axis.TickDir = 'out';
                obj.Handles.SubTrace(Pr).Axis.YColor = 'none';
                obj.Handles.SubTrace(Pr).Axis.Box = 'off';

                % Prepare cells selection arrays
                obj.SelectedCells{Pr} = false(size(obj.FootPrints(Pr).FootPrints,3),4); % Nothing selected by default
                obj.CellCategories{Pr} = 4*ones(size(obj.FootPrints(Pr).FootPrints,3),1); % All neutral
                if ~isempty(Categories{Pr})
                    for ROI = 1 : size(Categories{Pr},1)
%                         CellIndx = obj.GetReverseMatchingCell(Pr,ROI);
                        obj.CellCategories{Pr}(ROI) = Categories{Pr}(ROI);
                    end
                end

                % Add ui elements (hidden)
                Left = obj.Handles.SubROI(Pr).Axis.Position(1);
                Bottom = obj.Handles.SubROI(Pr).Axis.Position(2) + 1.01*obj.Handles.SubROI(Pr).Axis.Position(4);
                obj.Handles.SubROI(Pr).Selection = uicontrol('Style','checkbox','String','','FontSize', obj.FontScaling*14,'FontName','Arial','FontWeight','normal',...
                    'Callback',{@(src,evt)obj.SelectParadigm(src,evt)},'Units','Normalized','Position',[Left Bottom 0.01 0.03],'UserData',Pr,'Value',1,...
                    'Enable','off','Visible','off');
                Left = obj.Handles.SubROI(Pr).Axis.Position(1)+0.01;
                uicontrol('Style','text','String',obj.Paradigms{Pr},'FontSize', obj.FontScaling*11,'FontName','Arial','FontWeight','bold',...
                    'Units','normalized','Position',[Left 1.004*Bottom obj.Handles.SubROI(Pr).Axis.Position(3)/2 0.02],'HorizontalAlignment','left');
                Left = obj.Handles.SubROI(Pr).Axis.Position(1)+0.55 * obj.Handles.SubROI(Pr).Axis.Position(3);
                obj.Handles.SubROI(Pr).R = uicontrol('Style','checkbox','String','R','FontSize', obj.FontScaling*12,'FontName','Arial','FontWeight','normal',...
                    'Callback',{@(src,evt)obj.SelectRGB(src,evt)},'Units','Normalized','Position',[Left+0.01 Bottom 0.02 0.03],'UserData',Pr,'Value',0,...
                    'Enable','off','Visible','off');
                obj.Handles.SubROI(Pr).G = uicontrol('Style','checkbox','String','G','FontSize', obj.FontScaling*12,'FontName','Arial','FontWeight','normal',...
                    'Callback',{@(src,evt)obj.SelectRGB(src,evt)},'Units','Normalized','Position',[Left+0.03 Bottom 0.02 0.03],'UserData',Pr,'Value',1,...
                    'Enable','off','Visible','off');
                obj.Handles.SubROI(Pr).B = uicontrol('Style','checkbox','String','B','FontSize', obj.FontScaling*12,'FontName','Arial','FontWeight','normal',...
                    'Callback',{@(src,evt)obj.SelectRGB(src,evt)},'Units','Normalized','Position',[Left+0.05 Bottom 0.02 0.03],'UserData',Pr,'Value',0,...
                    'Enable','off','Visible','off');
            end
            obj.SelectedSingleCell = []; % Should be [ParadigmIndx,CellIndx];

            linkaxes([obj.Handles.SubROI(:).Axis],'xy')

            obj.Temp.PreviousPlotMode = '';

            % Interactivity
            uicontrol('Style','text','String','Interactivity','FontSize', obj.FontScaling*14,'FontName','Arial','FontWeight','bold',...
                'Units','normalized','Position',[0.025 0.9 0.1 0.03],'HorizontalAlignment','left');

            uicontrol('Style','text','String','Plot mode','FontSize', obj.FontScaling*14,'FontName','Arial','FontWeight','normal',...
                'Units','normalized','Position',[0.025 0.86 0.1 0.03],'HorizontalAlignment','left');

            obj.Handles.PlotMode.Paradigms = uicontrol('Style','checkbox','String',' Paradigms','FontSize', obj.FontScaling*14,'FontName','Arial','FontWeight','normal',...
                'Enable','off','Callback',{@(src,evt)obj.SetPlotMode(src,evt)},'Units','Normalized','Position',[0.035 0.83 0.1 0.03],'UserData','Paradigms','Value',1);

            obj.Handles.PlotMode.Select_All = uicontrol('Style','PushButton','String',' Select all','FontSize', obj.FontScaling*14,'FontName','Arial','FontWeight','normal',...
                'Enable','on','Callback',{@(~,~)obj.Select_All},'Units','Normalized','Position',[0.12 0.83 0.05 0.03],'UserData','Paradigms');

            obj.Handles.PlotMode.SingleCell = uicontrol('Style','checkbox','String',' Single cell','FontSize', obj.FontScaling*14,'FontName','Arial','FontWeight','normal',...
                'Enable','off','Callback',{@(src,evt)obj.SetPlotMode(src,evt)},'Units','Normalized','Position',[0.035 0.8 0.1 0.03],'UserData','SingleCell');

            % Free cells: we can choose the color / set the current cell as
            % good or bad
            obj.Handles.PlotMode.FreeCells = uicontrol('Style','checkbox','String',' Free cells','FontSize', obj.FontScaling*14,'FontName','Arial','FontWeight','normal',...
                'Enable','off','Callback',{@(src,evt)obj.SetPlotMode(src,evt)},'Units','Normalized','Position',[0.035 0.77 0.1 0.03],'UserData','FreeCells');
            obj.Handles.PlotMode.R = uicontrol('Style','checkbox','String',' R','FontSize', obj.FontScaling*12,'FontName','Arial','FontWeight','normal',...
                'Callback',{@(src,evt)obj.SelectFreeCellsRGB(src,evt)},'Units','Normalized','Position',[0.04 0.74 0.02 0.03],'UserData','R','Value',0,...
                'Enable','off','Visible','on');
            obj.Handles.PlotMode.G = uicontrol('Style','checkbox','String',' G','FontSize', obj.FontScaling*12,'FontName','Arial','FontWeight','normal',...
                'Callback',{@(src,evt)obj.SelectFreeCellsRGB(src,evt)},'Units','Normalized','Position',[0.04 0.71 0.02 0.03],'UserData','G','Value',1,...
                'Enable','off','Visible','on');
            obj.Handles.PlotMode.B = uicontrol('Style','checkbox','String',' B','FontSize', obj.FontScaling*12,'FontName','Arial','FontWeight','normal',...
                'Callback',{@(src,evt)obj.SelectFreeCellsRGB(src,evt)},'Units','Normalized','Position',[0.04 0.68 0.02 0.03],'UserData','B','Value',0,'Value',0,...
                'Enable','off','Visible','on');

            obj.Handles.PlotMode.Reset = uicontrol('Style','pushbutton','String','Reset','FontSize', obj.FontScaling*13,'FontName','Arial','FontWeight','normal',...
                'Callback',{@(~,~)obj.Reset},'Units','Normalized','Position',[0.025 0.63 0.08 0.03],...
                'Enable','on','Visible','on');
            obj.Handles.PlotMode.Unselect_All = uicontrol('Style','pushbutton','String','Unselect all','FontSize', obj.FontScaling*13,'FontName','Arial','FontWeight','normal',...
                'Callback',{@(~,~)obj.Unselect_All},'Units','Normalized','Position',[0.025 0.60 0.08 0.03],...
                'Enable','on','Visible','on');

            obj.Handles.DisplayNumbers = uicontrol('Style','checkbox','String',' Display cell numbers','FontSize', obj.FontScaling*12,'FontName','Arial','FontWeight','normal',...
                'Callback',{@(src,evt)obj.DisplayNumbersCB(src,evt)},'Units','Normalized','Position',[0.035 0.50 0.1 0.03],'Value',0,...
                'Enable','on','Visible','on');
            obj.Handles.DisplayContours = uicontrol('Style','checkbox','String',' Display contours','FontSize', obj.FontScaling*12,'FontName','Arial','FontWeight','normal',...
                'Callback',{@(src,evt)obj.DisplayContoursCB(src,evt)},'Units','Normalized','Position',[0.035 0.47 0.1 0.03],'Value',0,...
                'Enable','on','Visible','on');

            obj.Handles.DisplayTraces = uicontrol('Style','checkbox','String',' Display Traces','FontSize', obj.FontScaling*12,'FontName','Arial','FontWeight','normal',...
                'Callback',{@(src,evt)obj.DisplayTracesCB(src,evt)},'Units','Normalized','Position',[0.035 0.44 0.1 0.03],'Value',obj.DisplayTraces,...
                'Enable','on','Visible','on');
            uicontrol('Style','text','String','Smoothing','FontSize', obj.FontScaling*12,'FontName','Arial','FontWeight','normal',...
                'Units','normalized','Position',[0.035 0.405 0.1 0.03],'HorizontalAlignment','left');
            obj.Handles.SmoothingTraces = uicontrol('Style','edit','String',num2str(obj.Smoothing),'FontSize', obj.FontScaling*12,'FontName','Arial','FontWeight','normal',...
                'Callback',{@(src,evt)obj.SetSmoothingValue(src,evt)},'Units','Normalized','Position',[0.08 0.41 0.025 0.03],...
                'Enable','on','Visible','on');

            obj.Handles.KeepTracks = uicontrol('Style','checkbox','String',' Keep track','FontSize', obj.FontScaling*12,'FontName','Arial','FontWeight','normal',...
                'Callback',{@(src,evt)obj.KeepTrackCB(src,evt)},'Units','Normalized','Position',[0.035 0.35 0.1 0.03],'Value',obj.KeepTrack,...
                'Enable','on','Visible','on');

            obj.Handles.LinkTraces = uicontrol('Style','checkbox','String',' Link traces'' x axes','FontSize', obj.FontScaling*12,'FontName','Arial','FontWeight','normal',...
                'Callback',{@(src,evt)obj.LinkTracesAxesCB(src,evt)},'Units','Normalized','Position',[0.035 0.26 0.1 0.03],'Value',obj.LinkTraces,...
                'Enable','on','Visible','on');

            obj.Handles.Export = uicontrol('Style','pushbutton','String','Export good/bad','FontSize', obj.FontScaling*12,'FontName','Arial','FontWeight','normal',...
                'Callback',{@(~,~)obj.Export},'Units','Normalized','Position',[0.035 0.2 0.1 0.03],'Value',obj.LinkTraces,...
                'Enable','on','Visible','on');



            % Prepare paradigms selection arrays
            obj.SelectedParadigms = false(numel(obj.Paradigms),4);
            obj.SelectedParadigms(:,[1 3]) = true; % All selected by default (column one) & all green (column two is true)

            obj.AllColors = [obj.Bad_Contour_Color; obj.Neutral_Contour_Color; obj.Good_Contour_Color; obj.NonSelected_Contour_Color];
%             obj.AllColors_Selection = [obj.Bad_Contour_Color; obj.Selected_Contour_Color; obj.Good_Contour_Color; ];
            obj.ApplyPlotMode;

            obj.ToggleAxes;
        end

        function Reset(obj)
            switch obj.PlotMode
                case 'Paradigms'
                    Fields = {'Selection','R','G','B'};
                    for Pr = 1 : numel(obj.Paradigms)
                        obj.SelectedCells{Pr} = false(size(obj.FootPrints(Pr).FootPrints,3),4); % Nothing selected by default
                        obj.CellCategories{Pr} = 4*ones(size(obj.FootPrints(Pr).FootPrints,3),1); % All neutral
                        set([obj.Handles.SubROI(Pr).Contours],'LineWidth',obj.NonSelected_Contour_Width)
                        set([obj.Handles.SubROI(Pr).Contours],'EdgeColor',obj.AllColors(4,:));
                        for FieldsPr = 1 : numel(Fields)
                            obj.Handles.SubROI(Pr).(Fields{FieldsPr}).Enable = 'on';
                            obj.Handles.SubROI(Pr).(Fields{FieldsPr}).Visible = 'on';
                            if strcmpi(Fields{FieldsPr},'G')
                                obj.Handles.SubROI(Pr).(Fields{FieldsPr}).Value = 1;
                            else
                                obj.Handles.SubROI(Pr).(Fields{FieldsPr}).Value = 0;
                            end
                        end
                    end
                    obj.SelectedParadigms(:,1) = 1;
                    set([obj.Handles.SubROI(:).Selection],'Value',true)
                case 'SingleCell'
                    set([obj.Handles.SubTrace.Trace],'Visible',false);
                    obj.SelectedSingleCell = [];
                case 'FreeCell'
                    for Pr = 1 : numel(obj.Paradigms)
                        obj.SelectedCells{Pr} = false(size(obj.FootPrints(Pr).FootPrints,3),4); % Nothing selected by default
                    end
                    for Chan = 1 : numel(obj.Chan_Colors)
                        if strcmpi(obj.Chan_Colors{Chan},'G')
                            obj.Handles.PlotMode.(obj.Chan_Colors{Chan}) = 1;
                        else
                            obj.Handles.PlotMode.(obj.Chan_Colors{Chan}) = 0;
                        end
                    end
            end
            obj.UpdateFootPrints;
            obj.UpdateTraces;
        end

        function Select_All(obj)
            obj.SelectedParadigms(:,1) = 1;
            set([obj.Handles.SubROI(:).Selection],'Value',true)
            obj.UpdateFootPrints;
        end

        function Unselect_All(obj)
            switch obj.PlotMode
                case 'Paradigms'
                    obj.SelectedParadigms(:,1) = 0;
                    set([obj.Handles.SubROI(:).Selection],'Value',false)
                case 'SingleCell'
                    set([obj.Handles.SubTrace.Trace],'Visible',false);
                    obj.SelectedSingleCell = [];
                case 'FreeCell'
                    for Pr = 1 : numel(obj.Paradigms)
                        obj.SelectedCells{Pr} = false(size(obj.FootPrints(Pr).FootPrints,3),4); % Nothing selected by default
                    end
                    for Chan = 1 : numel(obj.Chan_Colors)
                        if strcmpi(obj.Chan_Colors{Chan},'G')
                            obj.Handles.PlotMode.(obj.Chan_Colors{Chan}) = 1;
                        else
                            obj.Handles.PlotMode.(obj.Chan_Colors{Chan}) = 0;
                        end
                    end
            end
            obj.SelectedROIs = [];
            obj.UpdateFootPrints;
            obj.UpdateTraces;
        end

       
        function UpdateFootPrints(obj,Indx)
            if obj.PauseRefreshing
                return
            end
            if nargin>1 && ~isempty(Indx)
                Looping = Indx;
            else
                Looping = 1 : numel(obj.Paradigms);
            end
            % Prepare the color array based on paradigms/cells
            % selection and individual colors
            ColorArray = cell(numel(Looping),1);
            for SubL = 1 : numel(Looping)
                Pr = Looping(SubL);
                ColorArray{Pr} = obj.FootPrints(Pr).Base;
            end
            switch obj.PlotMode
                case 'Paradigms'
                    % We need to check which paradigms are selected for each color
                    Chan_List = {'R','G','B'};
                    for Chan = 1 : numel(Chan_List)
                        Paradigms_Chan = obj.SelectedParadigms(:,1) & obj.SelectedParadigms(:,Chan+1);
                        if ~all(Paradigms_Chan==0)
                            CellIndx = all(obj.cellRegistered.cell_to_index_map(:,Paradigms_Chan),2);
                            for SubL = 1 : numel(Looping)
                                Pr = Looping(SubL);
                                CellIndx_Pr = obj.cellRegistered.cell_to_index_map(CellIndx,Pr);
                                CellIndx_Pr = CellIndx_Pr(CellIndx_Pr~=0);
                                ColorArray{Pr}(:,:,Chan) = ColorArray{Pr}(:,:,Chan) + sum(obj.FootPrints(Pr).FootPrints(:,:,CellIndx_Pr),3);
                            end
                        end
                    end
                case 'SingleCell'
                    if strcmpi(obj.SingleCellMode,'FootPrints')
                        % We need to check which paradigms are selected for each color
                        CellIndx = obj.cellRegistered.cell_to_index_map(:,obj.SelectedSingleCell(1)) == obj.SelectedSingleCell(2);
                        for SubL = 1 : numel(Looping)
                            Pr = Looping(SubL);
                            CellIndx_Pr = obj.cellRegistered.cell_to_index_map(CellIndx,Pr);
                            CellIndx_Pr = CellIndx_Pr(CellIndx_Pr~=0);
                            ColorArray{Pr}(:,:,2) = ColorArray{Pr}(:,:,2) + sum(obj.FootPrints(Pr).FootPrints(:,:,CellIndx_Pr),3);
                        end
                    else

                        return
                    end
                case 'FreeCell'
                    % Loop through colors
                    Chan_List = {'R','G','B'};
                    for Chan = 1 : numel(Chan_List)
                        Chan_Indx = find(obj.SelectedCells(:,Chan+2));
                        UniqueParadigms = unique(obj.SelectedCells(Chan_Indx,1));
                        % Loop through paradigms
                        for Pg = 1 : numel(UniqueParadigms)
                            Paradigm_Indx = obj.SelectedCells(:,Chan+2) & obj.SelectedCells(Chan_Indx,1)==(UniqueParadigms(Pg));
                            CellIndx = obj.SelectedCells(Paradigm_Indx,2);
                            CellIndx = CellIndx(CellIndx~=0);
                            [~, CellIndx] = ismember(obj.SelectedCells(CellIndx,2),obj.cellRegistered.cell_to_index_map(:,UniqueParadigms(Pg)));
                            for SubL = 1 : numel(Looping)
                                Pr = Looping(SubL);
                                CellIndx_Pr = obj.cellRegistered.cell_to_index_map(CellIndx,Pr);
                                ColorArray{Pr}(:,:,Chan) = ColorArray{Pr}(:,:,Chan) + sum(obj.FootPrints(Pr).FootPrints(:,:,CellIndx_Pr),3);
                            end
                        end
                    end
            end
            for SubL = 1 : numel(Looping)
                Pr = Looping(SubL);
                if ~isempty(obj.Handles.SubROI(Pr).FootPrints)
                    % delete(obj.Handles.SubROI(Pr).FootPrints);

                    obj.Handles.SubROI(Pr).FootPrints.CData = ColorArray{Pr};
                else
                    obj.Handles.SubROI(Pr).FootPrints = imagesc(ColorArray{Pr},'Parent',obj.Handles.SubROI(Pr).Axis);
                    if ~ishold(obj.Handles.SubROI(Pr).Axis)
                        hold(obj.Handles.SubROI(Pr).Axis,'on')
                    end
                    if obj.DisplayContours
                        DisplayOpt = 'on';
                    else
                        DisplayOpt = 'off';
                    end              
                    
                    obj.Handles.SubROI(Pr).Contours = arrayfun(@(x) fill(obj.Contours{Pr}{x}(1,:),obj.Contours{Pr}{x}(2,:),obj.AllColors(obj.CellCategories{Pr}(x),:),'FaceColor','w','FaceAlpha',eps,'PickableParts','all','EdgeColor',obj.AllColors(obj.CellCategories{Pr}(x),:),'LineWidth',obj.NonSelected_Contour_Width,'ButtonDownFcn',{@(src,evt)obj.ContourClick(src,evt)},'UserData',[Pr,x],'Visible',DisplayOpt,'Parent',obj.Handles.SubROI(Pr).Axis),1:numel(obj.Contours{Pr}));
                    if obj.DisplayNumbers
                        DisplayOpt = 'on';
                    else
                        DisplayOpt = 'off';
                    end   
                    obj.Handles.SubROI(Pr).Numbers =  arrayfun(@(x) text(nanmean(obj.Contours{Pr}{x}(1,:)),nanmean(obj.Contours{Pr}{x}(2,:)),num2str(x),'FontName','Arial','FontSize', obj.FontScaling*10,'Color','k','FontWeight','bold','Parent',obj.Handles.SubROI(Pr).Axis,'HorizontalAlignment','center','Clipping','on','Visible',DisplayOpt),1:numel(obj.Contours{Pr}));

                    axis(obj.Handles.SubROI(Pr).Axis,'equal')
                    obj.Handles.SubROI(Pr).Axis.Color = 'k';
                    obj.Handles.SubROI(Pr).Axis.YDir = 'normal';
                    obj.Handles.SubROI(Pr).Axis.XColor = 'none';
                    obj.Handles.SubROI(Pr).Axis.XDir = 'normal';
                    obj.Handles.SubROI(Pr).Axis.YDir = 'normal';
                    obj.Handles.SubROI(Pr).Axis.YColor = 'none';
                    obj.Handles.SubROI(Pr).Axis.YLim(1) = obj.Handles.SubROI(Pr).Axis.YLim(1)-0.01*diff(obj.Handles.SubROI(Pr).Axis.YLim);
                    obj.Handles.SubROI(Pr).Axis.YLim(2) = obj.Handles.SubROI(Pr).Axis.YLim(2)+0.01*diff(obj.Handles.SubROI(Pr).Axis.YLim);
                    drawnow
                    obj.Handles.SubROI(Pr).Axis.XLimMode = 'manual';
                    obj.Handles.SubROI(Pr).Axis.YLimMode = 'manual';
                    axis(obj.Handles.SubROI(Pr).Axis,'equal')
                    if obj.EnableCrossHair
                        obj.Handles.SubROI(Pr).CrossHair = plot(obj.CrossHair_Coordinates(1),obj.CrossHair_Coordinates(2),'r.','LineWidth',1,'MarkerSize',5,'Parent',obj.Handles.SubROI(Pr).Axis,'PickableParts','none');
                        set(obj.Figure, 'WindowButtonMotionFcn', {@(~,~)obj.MoveCrossHair});
                    else
                        set(obj.Figure, 'WindowButtonMotionFcn','');
                    end
                end
            end
        end
        
        function DisplayTracesCB(obj,src,~)
            if src.Value
                obj.DisplayTraces = true;
            else
                obj.DisplayTraces = false;
            end
            obj.ToggleAxes;
            obj.UpdateTraces;
        end

        function LinkTracesAxes(obj)
            if obj.LinkTraces
                linkaxes([obj.Handles.SubTrace(:).Axis],'x')
                set([obj.Handles.SubTrace(:).Axis],'XLimMode','manual')
                GCA = gca;
                set([obj.Handles.SubTrace(:).Axis],'XLim',GCA.XLim)
                drawnow
            else
                linkaxes([obj.Handles.SubTrace(:).Axis],'')
                set([obj.Handles.SubTrace(:).Axis],'XLimMode','auto')
                drawnow
                for Pr = 1 : numel(obj.Paradigms)
                    obj.Handles.SubTrace(Pr).Axis.XLim = obj.Handles.SubTrace(Pr).Axis.UserData;
                end
            end
        end
        
        function LinkTracesAxesCB(obj,src,~)
            if src.Value
                % We need to add it
                obj.LinkTraces = 1;
            else
                % We need to remove it
                obj.LinkTraces = 0;
            end
            obj.LinkTracesAxes;
        end

        function ToggleAxes(obj)
            if obj.DisplayTraces
                set([obj.Handles.SubTrace(:).Axis],'Visible','on','HandleVisibility', 'on')
                for Pr = 1 : numel(obj.Paradigms)
                    set(obj.Handles.SubROI(Pr).Axis,'Position',obj.AxesPosition{Pr}.Split.SubROI)
                    set(obj.Handles.SubROI(Pr).BackGroundAxis,'Position',obj.AxesPosition{Pr}.Split.SubROI)
                end
            else
                set([obj.Handles.SubTrace(:).Axis],'Visible','off','HandleVisibility', 'off')
                for Pr = 1 : numel(obj.Paradigms)
                    set(obj.Handles.SubROI(Pr).Axis,'Position',obj.AxesPosition{Pr}.Full)
                    set(obj.Handles.SubROI(Pr).BackGroundAxis,'Position',obj.AxesPosition{Pr}.Full)
                end
            end
        end
        function KeepTrackCB(obj,src,~)
            if src.Value
                obj.KeepTrack = true;
            else
                obj.KeepTrack = false;
            end
        end


        function DisplayContoursCB(obj,src,~)
            if src.Value
                obj.DisplayContours = true;
                set([obj.Handles.SubROI(:).Contours],'Visible','on')
            else
                obj.DisplayContours = false;
                set([obj.Handles.SubROI(:).Contours],'Visible','off')
            end
        end

        function DisplayNumbersCB(obj,src,~)
            if src.Value
                obj.DisplayNumbers = true;
                set([obj.Handles.SubROI(:).Numbers],'Visible','on')
            else
                obj.DisplayNumbers = false;
                set([obj.Handles.SubROI(:).Numbers],'Visible','off')
            end
        end
        function MoveCrossHair(obj,~,~)
            p = get(0,'PointerLocation');
            set([obj.Handles.SubROI(:).Axis],'XLimMode','manual');
            set([obj.Handles.SubROI(:).Axis],'YLimMode','manual');
            % Compute figure offset of mouse pointer in pixels
            figPos = getpixelposition(obj.Figure);
            x = (p(1)-figPos(1));
            y = (p(2)-figPos(2));
            for Ax_Handle = 1:numel(obj.Paradigms)
                % If descendant contains the mouse pointer position, exit
                r = getpixelposition(obj.Handles.SubROI(Ax_Handle).Axis);  % Note: cache this for improved performance
                if (x>r(1)) && (x<r(1)+r(3)) && (y>r(2)) && (y<r(2)+r(4))
                    break
                end
            end
            if ~isempty(Ax_Handle)
                TempC = get(obj.Handles.SubROI(Ax_Handle).Axis,'currentpoint');
                obj.CrossHair_Coordinates = TempC(1,[1 2]);
%                 for Pr = 1 : numel(obj.Paradigms)
                    %                     delete(obj.Handles.SubROI(Pr).CrossHair)
                    %                     obj.Handles.SubROI(Pr).CrossHair = plot(obj.CrossHair_Coordinates(1),obj.CrossHair_Coordinates(2),'r+','LineWidth',0.5,'MarkerSize',8,'Parent',obj.Handles.SubROI(Pr).Axis);
                    set([obj.Handles.SubROI(:).CrossHair], 'XData', obj.CrossHair_Coordinates(1), 'YData', obj.CrossHair_Coordinates(2));
%                 end
                drawnow
            end
            if ~obj.EnableCrossHair
                set(obj.Figure, 'WindowButtonMotionFcn','');
                obj.CrossHair_Coordinates = [NaN NaN];
%                 for Pr = 1 : numel(obj.Paradigms)
                    %                     delete(obj.Handles.SubROI(Pr).CrossHair)
                    %                     obj.Handles.SubROI(Pr).CrossHair = plot(obj.CrossHair_Coordinates(1),obj.CrossHair_Coordinates(2),'r+','LineWidth',0.5,'MarkerSize',8,'Parent',obj.Handles.SubROI(Pr).Axis);
                    set([obj.Handles.SubROI(:).CrossHair], 'XData', obj.CrossHair_Coordinates(1), 'YData', obj.CrossHair_Coordinates(2));
%                 end
            end
        end

        function SetPlotMode(obj,src,~)
            Clicked = src.UserData;
            PM = fieldnames(obj.Handles.PlotMode);
            for Fi = 1 : numel(PM)
                if ~strcmpi(Clicked,PM{Fi})
                    obj.Handles.PlotMode.(PM{Fi}).Value = 0;
                else
                    obj.Handles.PlotMode.(PM{Fi}).Value = 1;
                end
            end
            obj.Temp.PreviousPlotMode = obj.PlotMode;
            obj.PlotMode = Clicked;
            obj.ApplyPlotMode;
        end

        function SelectParadigm(obj,src,~)
            Clicked = src.UserData;
            if src.Value
                % We need to add it
                obj.SelectedParadigms(Clicked,1) = 1;
            else
                % We need to remove it
                obj.SelectedParadigms(Clicked,1) = 0;
            end
            obj.UpdateFootPrints;
        end

        function SelectRGB(obj,src,~)
            Clicked = src.UserData;
            IndexChan = find(arrayfun(@(x) obj.Handles.SubROI(Clicked).(obj.Chan_Colors{x}).Value ,1:numel(obj.Chan_Colors)));
            if ~isempty(IndexChan)
                obj.SelectedParadigms(Clicked,2:end) = 0;
                obj.SelectedParadigms(Clicked,IndexChan+1) = 1;
            else
                src.Value = 1;
            end

            obj.UpdateFootPrints;
        end

        function ApplyPlotMode(obj)
            if strcmpi(obj.PlotMode,obj.Temp.PreviousPlotMode)
                return
            end

            % Clean from previous mode
            if ~isempty(obj.Temp.PreviousPlotMode)
                switch obj.Temp.PreviousPlotMode
                    case 'Paradigms'

                        Fields = {'Selection','R','G','B'};
                        for Pr = 1 : numel(obj.Paradigms)
                            % Hide tickboxes / update footprints
                            for FieldsPr = 1 : numel(Fields)
                                obj.Handles.SubROI(Pr).(Fields{FieldsPr}).Enable = 'off';
                                obj.Handles.SubROI(Pr).(Fields{FieldsPr}).Visible = 'off';
                            end
                        end

                    case 'SingleCell'

                    case 'FreeCell'

                end

            end


            % Setup new layout
            switch obj.PlotMode
                case 'Paradigms'
                    Fields = {'Selection','R','G','B'};
                    % Hide tickboxes
                    for Pr = 1 : numel(obj.Paradigms)
                        for FieldsPr = 1 : numel(Fields)
                            obj.Handles.SubROI(Pr).(Fields{FieldsPr}).Enable = 'on';
                            obj.Handles.SubROI(Pr).(Fields{FieldsPr}).Visible = 'on';
                        end
                    end
                case 'SingleCell'

                case 'FreeCell'

            end

            obj.UpdateFootPrints;
        end

        function ContourClick(obj,src,~)
            Clicked = src.UserData;
            if ~isempty(obj.SelectedROIs)
                if obj.SelectedROIs(1)==Clicked(1) && obj.SelectedROIs(2)==Clicked(2)
                    obj.DeselectCell(Clicked);
                else
                    obj.DeselectCell(obj.SelectedROIs)
                    obj.SelectCell(Clicked);
                end
            else
                obj.SelectedROIs = Clicked;
                obj.SelectCell(Clicked);
            end
            % Update traces
            obj.UpdateTraces;
        end

        function SelectCell(obj,Cell)
            obj.SelectedROIs = Cell;
            obj.Handles.SubROI(Cell(1)).Contours(Cell(2)).LineWidth = obj.Selected_Contour_Width;
            ColIndx = obj.CellCategories{Cell(1)}(Cell(2));
            if ColIndx == 4
                ColIndx = 2;
            end
            if obj.DisplayContours
                DisplayOpt = 'on';
            else
                DisplayOpt = 'off';
            end
            for Pr = 1 : numel(obj.Paradigms)
                CellID = obj.GetMatchingCell(Cell,Pr);
                if ~isempty(CellID)
                    obj.Handles.SubROI(Pr).Contours(CellID).EdgeColor = obj.AllColors(ColIndx,:);
                else
                    % We plot the ghost of the selected cell
                    obj.Handles.SubROI(Pr).GhostContours = ...
                        fill(obj.Contours{Cell(1)}{Cell(2)}(1,:),obj.Contours{Cell(1)}{Cell(2)}(2,:),obj.MissingColor,'FaceColor','w','FaceAlpha',eps,'PickableParts','all','EdgeColor',obj.MissingColor,'LineWidth',obj.NonSelected_Contour_Width,'Visible',DisplayOpt,'Parent',obj.Handles.SubROI(Pr).Axis);
                end
            end
        end

        function DeselectCell(obj,Cell)
            obj.SelectedROIs = [];
            obj.Handles.SubROI(Cell(1)).Contours(Cell(2)).LineWidth = obj.NonSelected_Contour_Width;
            if obj.CellCategories{Cell(1)}(Cell(2)) == 4 && obj.KeepTrack
                obj.CellCategories{Cell(1)}(Cell(2)) = 2;
            end
            for Pr = 1 : numel(obj.Paradigms)
                CellID = obj.GetMatchingCell(Cell,Pr);
                if ~isempty(CellID)
                    if obj.CellCategories{Pr}(CellID) == 4 && obj.KeepTrack
                        obj.CellCategories{Pr}(CellID) = 2;
                    end
                    obj.Handles.SubROI(Pr).Contours(CellID).EdgeColor = obj.AllColors(obj.CellCategories{Pr}(CellID),:);
                else
                    delete(obj.Handles.SubROI(Pr).GhostContours)
                end
            end
        end


        function NextCell(obj)
            if isempty(obj.SelectedROIs)
                return
            end
            if obj.SelectedROIs(2)==numel(obj.Contours{obj.SelectedROIs(1)})
                NewCell = [obj.SelectedROIs(1) 1];
            else
                NewCell = [obj.SelectedROIs(1) obj.SelectedROIs(2)+1];
            end
            obj.DeselectCell(obj.SelectedROIs)
            obj.SelectCell(NewCell);
            % Update traces
            obj.UpdateTraces;
        end

        function PreviousCell(obj)
            if isempty(obj.SelectedROIs)
                return
            end
            if obj.SelectedROIs(2)==1
                NewCell = [obj.SelectedROIs(1) numel(obj.Contours{obj.SelectedROIs(1)})];
            else
                NewCell = [obj.SelectedROIs(1) obj.SelectedROIs(2)-1];
            end
            obj.DeselectCell(obj.SelectedROIs)
            obj.SelectCell(NewCell);
            % Update traces
            obj.UpdateTraces;
        end

        function CategorizeCell(obj,Cat)
            switch Cat
                case 'Up'
                    if obj.CellCategories{obj.SelectedROIs(1)}(obj.SelectedROIs(2)) ~= 3
                        if obj.CellCategories{obj.SelectedROIs(1)}(obj.SelectedROIs(2)) == 4
                            obj.CellCategories{obj.SelectedROIs(1)}(obj.SelectedROIs(2)) = 3;
                        else
                            obj.CellCategories{obj.SelectedROIs(1)}(obj.SelectedROIs(2)) = obj.CellCategories{obj.SelectedROIs(1)}(obj.SelectedROIs(2))+1;
                        end
                    end
                case 'Down'
                    if obj.CellCategories{obj.SelectedROIs(1)}(obj.SelectedROIs(2)) ~= 1
                        if obj.CellCategories{obj.SelectedROIs(1)}(obj.SelectedROIs(2)) == 4
                            obj.CellCategories{obj.SelectedROIs(1)}(obj.SelectedROIs(2)) = 1;
                        else
                            obj.CellCategories{obj.SelectedROIs(1)}(obj.SelectedROIs(2)) = obj.CellCategories{obj.SelectedROIs(1)}(obj.SelectedROIs(2))-1;
                        end
                    end
            end
            obj.SetCategories(obj.SelectedROIs);
            drawnow
        end

        function SetCategories(obj,Cell)
            for Pr = 1 : numel(obj.Paradigms)
                CellID = obj.GetMatchingCell(Cell,Pr);
                obj.CellCategories{Pr}(CellID) = obj.CellCategories{Cell(1)}(Cell(2));
                if ~isempty(CellID)
                    obj.Handles.SubROI(Pr).Contours(CellID).EdgeColor = obj.AllColors(obj.CellCategories{Pr}(CellID),:);
                end
            end
        end

        function Export(obj)
            for Pr = 1 : numel(obj.Paradigms)
                % Get CNMFE folder
                TempFolder = fileparts(obj.CNMFE_Files{Pr});
                LoadedCNMFE = load(obj.CNMFE_Files{Pr});
%                 IndxPr = obj.cellRegistered.cell_to_index_map(:,Pr)>0;
%                 CellsPr = obj.cellRegistered.cell_to_index_map(IndxPr,Pr);
%                 [~,Indx] = sort(CellsPr);
%                 Categories = obj.CellCategories{Pr}(Indx,1);
                Categories = obj.CellCategories{Pr};
%                 Categories = obj.CellCategories{Pr};
                save([TempFolder filesep LoadedCNMFE.PreProcessing.Basename '_CellCategories.mat'],'Categories')
            end
        end


        function CellID = GetMatchingCell(obj,Cell,TargetParadigm)
            Indx = obj.cellRegistered.cell_to_index_map(:,Cell(1)) == Cell(2);
            CellID = obj.cellRegistered.cell_to_index_map(Indx,TargetParadigm);
            if CellID == 0
                CellID = [];
            end
        end

        function CellID = GetReverseMatchingCell(obj,Paradigm,CNMFE_Cell)
            CellID = find(obj.cellRegistered.cell_to_index_map(:,Paradigm) == CNMFE_Cell);
        end

        function PressKeyCB(obj,~,evt)
            switch evt.Key
                case 'z'
                    %                     if ~obj.Zoom
                    %                         obj.Zoom = true;
                    %                         hManager = uigetmodemanager(obj.Figure);
                    %                         [hManager.WindowListenerHandles.Enabled] = deal(false);  %% UNDOCUMENTED
                    %                         set(obj.Figure, 'WindowKeyPressFcn', [], 'KeyPressFcn',  {@(src,evt)obj.PressKeyCB(src,evt)});
                    %                         zoom on
                    %                     end
                case 'return'
                    %                     if obj.Editing
                    %                         obj.RangeEditCB
                    %                     end
                case 'space'
                    if obj.PauseRefreshing
                        obj.PauseRefreshing = false;
                    else
                        obj.PauseRefreshing = true;
                    end
                case 'backspace'

                case 'delete'
                    if ~isempty(obj.SelectedROIs)
                        obj.CategorizeCell('Neutral')
                    end
                case 'rightarrow'
                    obj.NextCell;
                case 'leftarrow'
                    obj.PreviousCell;
                case 'uparrow'
                    if ~isempty(obj.SelectedROIs)
                        obj.CategorizeCell('Up')
                    end
                case 'downarrow'
                    if ~isempty(obj.SelectedROIs)
                        obj.CategorizeCell('Down')
                    end
                case 'pagedown'
                    %                     if ~obj.Editing && ~obj.Dragging && ~isempty(obj.Reader)
                    %                         obj.PlayNextRange;
                    %                     end
                case 'pageup'
                    %                     if ~obj.Editing && ~obj.Dragging && ~isempty(obj.Reader)
                    %                         obj.PlayPreviousRange;
                    %                     end
            end
        end

        function SetSmoothingValue(obj,src,~)
            TempValue = src.String;
            if isempty(str2double(TempValue))
                src.String = obj.Smoothing;
            else
                obj.Smoothing = str2double(src.String);
                obj.GetTraces;
            end
        end

        function GetTraces(obj)
            for Pr = 1 : numel(obj.Paradigms)
                obj.Traces{Pr} = smoothdata(obj.TracesRaw{Pr},1,'gaussian',obj.Smoothing);
            end
            obj.UpdateTraces;
        end


        function UpdateTraces(obj)
            if obj.DisplayTraces
                if isempty(obj.SelectedROIs)
                    for Pr = 1 : numel(obj.Paradigms)
                        obj.Handles.SubTrace(Pr).Trace.XData = [NaN NaN];
                        obj.Handles.SubTrace(Pr).Trace.YData = [NaN NaN];
                    end
                else
                    for Pr = 1 : numel(obj.Paradigms)
                        CellID = obj.GetMatchingCell(obj.SelectedROIs,Pr);
                        if ~isempty(CellID)
                            obj.Handles.SubTrace(Pr).Trace.XData = obj.Times{Pr};
                            obj.Handles.SubTrace(Pr).Trace.YData = obj.Traces{Pr}(CellID,:);
                            obj.Handles.SubTrace(Pr).Axis.XLimMode = 'manual';
                            %                             obj.Handles.SubTrace(Pr).Axis.YLimMode = 'auto';
                            obj.Handles.SubTrace(Pr).Axis.UserData = [0 obj.Times{Pr}(end)];
                            obj.Handles.SubTrace(Pr).Axis.XLim = [0 obj.Times{Pr}(end)];
                        else
                            obj.Handles.SubTrace(Pr).Trace.XData = [NaN NaN];
                            obj.Handles.SubTrace(Pr).Trace.YData = [NaN NaN];
                            obj.Handles.SubTrace(Pr).Axis.UserData = [NaN NaN];
                        end
                    end
                    if obj.LinkTraces && ~any(isnan(obj.Handles.SubTrace(obj.SelectedROIs(1)).Axis.UserData))
                        set([obj.Handles.SubTrace(:).Axis],'XLim',obj.Handles.SubTrace(obj.SelectedROIs(1)).Axis.UserData)
                    end
                end
            end
        end
    end





    methods(Static)
        % Get contours of the all neurons (CNMFe function)
        function Coor = Get_Contours(Footprint,Threshold)
            d1 = size(Footprint,1);
            d2 = size(Footprint,2);

            A_temp = Footprint;
            % find the threshold for detecting nonzero pixels
            A_temp = A_temp(:);
            [temp,ind] = sort(A_temp(:).^2,'ascend');
            temp =  cumsum(temp);
            ff = find(temp > (1-Threshold)*temp(end),1,'first');
            thr_a = A_temp(ind(ff));

            A_temp = reshape(A_temp,[d1 d2]);

            % crop a small region for computing contours
            [tmp1, tmp2, ~] = find(A_temp);
            if isempty(tmp1)
                Coor = zeros(2,1);
                return;
            end
            rmin = max(1, min(tmp1)-3);
            rmax = min(d1, max(tmp1)+3);
            cmin = max(1, min(tmp2)-3);
            cmax = min(d2, max(tmp2)+3);
            A_temp = A_temp(rmin:rmax, cmin:cmax);

            if nnz(A_temp)>36
                l = bwlabel(medfilt2(A_temp>thr_a));
            else
                l = bwlabel(A_temp>=thr_a);
            end
            l_most = mode(l(l>0));
            if isnan(l_most)
                Coor = zeros(2, 1);
                return;
            end
            ind = (l==l_most);
            A_temp(ind) =  max(A_temp(ind), thr_a);
            A_temp(~ind) = min(A_temp(~ind), thr_a*0.99);

            pvpairs = { 'LevelList' , thr_a, 'ZData', A_temp};
            h = matlab.graphics.chart.primitive.Contour(pvpairs{:});
            temp = h.ContourMatrix;
            if isempty(temp)
                Coor = Explore_Longitudinal_Alignment.Get_Contours(Footprint,0.8*Threshold);
                return;
            else
                temp(:, 1) = temp(:, 2);
                temp = medfilt1(temp')';
                temp(:, 1) = temp(:, end);
                Coor = bsxfun(@plus, temp, [cmin-1; rmin-1]);
            end

        end
    end
end
