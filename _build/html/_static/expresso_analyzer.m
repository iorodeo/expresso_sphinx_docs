function eBrewer_r28(inputFile, outputDir, p_container, p_channel, p_close)
%% GUI-mode vs. batchmode:
% inputFile empty? guiMode : p_close? batchMode : generateMode

    GUI_MODE = 'GUI';
    BATCH_MODE = 'BATCH';
    GENERATE_MODE = 'GENERATE';
        
    currMode = [];
    
    currDir = pwd;
    if nargin < 1 || isempty(inputFile)
        currMode = GUI_MODE;
        inputFile = currDir; % '.';
    end
    if nargin < 2 || isempty(outputDir)
        outputDir = currDir; % '.';
    end
    if exist(outputDir) ~= 7 %#ok<EXIST>
        [status,message,messageid] = mkdir(outputDir);
        if ~status
            error(messageid, message);
        end
    end
    if nargin < 3
        p_container = [];
    end
    if nargin < 4
        p_channel = [];
    end
    if nargin < 5
        p_close = true;
    elseif isdeployed % delivers true, false or []
        p_close = str2num(p_close); %#ok<ST2NM> 
        p_close = true;
    end
    if isempty(p_close) || ischar(p_close) % came as string
        currMode = p_close;
    elseif isempty(currMode) % came as boolean
        if p_close
            currMode = BATCH_MODE;
        else
            currMode = GENERATE_MODE;
        end
    end
    
    common = myInitCommon();
    
    p_filter.container = p_container;
    p_filter.channel = p_channel;
    setCommon('p_filter', p_filter);
    setCommon('outputDir', outputDir);
            
    switch(upper(currMode))
        case GUI_MODE 
            fprintf('GUI-MODE!\n');
            myInitializeGUI(inputFile, outputDir);
        case BATCH_MODE 
            fprintf('BATCH-MODE!\n');
                        
            doBatch(inputFile, true);
            close(common);
        case GENERATE_MODE 
            fprintf('GENERATE-MODE!\n');
            
            doBatch(inputFile, false);
            myInitializeGUI(inputFile);
        otherwise
            fprintf('ERROR: Unknown Mode!\n');
    end

%     for fi = 1 : length(fileList)
%         doWork(fileList(fi).name);
%     end
            

% hfig = figure(   'Resize', 'off', ...
%                 'Visible', 'on', ...
%                 'CloseRequestFcn', @CloseRequestFcn, ...
%                 'Interruptible', 'off');
% 
% set(hfig, 'WindowButtonUpFcn', @buttonUpHandle);
% 
% videoPanel(1)=subplot(2,2,1, 'DrawMode','fast');
% axis off;
% videoPanel(2)=subplot(2,2,2, 'DrawMode','fast');
% axis off;
% controlPanel =subplot(2,2,3, 'DrawMode','fast');
% axis off;
% resultPanel  =subplot(2,2,4, 'DrawMode','fast');
% axis off;
% 
% set(videoPanel(1), 'pos', [0 0.5 0.5 0.5]);
% set(videoPanel(2), 'pos', [0.5 0.5 0.5 0.5]);
% set(resultPanel  , 'pos', [0.55 0.05 0.4 0.4]);
% 
% axes(videoPanel(1));
% videoImage = image(zeros(250,250,3));
% set(videoPanel(1), 'XTick',[],'YTick',[]);
% axes(videoPanel(2));
% contourImage = image(zeros(250,250,3));
% set(videoPanel(2), 'XTick',[],'YTick',[]);
% axes(resultPanel);
% set(resultPanel, 'XTick',[],'YTick',[]);
% 
% axes(controlPanel);
% buttonNext = uicontrol('Style', 'pushbutton', 'String', 'Compute', 'Callback', @funcCompute, 'ForegroundColor', 'b');
% rotationText = uicontrol('Style','text','String', 'Rotation');
% %rotationControl = uicontrol('Style','edit', 'String', '0','Callback', @funcRotate,'BackgroundColor','w');
% 
% 
% 
% g.fps = 15;
% g.curFrame = 1;
% %g.thickness = thickness;
% g.currVideo = 0;
% g.loadedVideo = 0;
% g.stat = 0;
% guidata(hfig, g);
% %SetupTimer(hfig);
% g = guidata(hfig);
% %loadNextVideo();    
    function common = myInitCommon()
        common = figure('visible','off');
        
        c.init_complete = false;
        
        c.currFileName = [];
        c.currFileData = [];
        c.currTimeInfo = [];
        c.currChannelPath = [];
        
        c.const.container = '(XP[0-9][0-9])|(device_[0-9])';
        c.const.channel = 'channel_[0-9]';
        c.const.datetimeName = 'sample_t';
        c.const.datetimePath = ['/' c.const.datetimeName];
        
        c.titlePrefix = '[eBrewer V0.5 R25]  ';
        
%        c.param.minEventVolume = 5; % [nl]
%        c.param.minEventInterMeal = 4; % [sec], otherwise merge
        c.param.minEventVolume = 4; % [nl] % new value 20140211
        c.param.minEventInterMeal = 0; % [sec], otherwise merge % new value 20140211
        c.param.minEventDuration = 4; % [sec], otherwise eliminate
        
        c.param.freqPrec = 2; % desired digits after comma
        c.param.maxLengthVariation = 0.2; % for bulk summaries, flag files lengths that are not within 20 % of other lengths
        c.param.maxFreqVariation = 0.2;
        
        c.devfakt.event = 1;
        c.devfakt.outlier = 4;
        
        c.rasterSumBinSec = 10; % CSC: find appropriate bin size, suggesting 10 sec

        % also default values for checkboxes
        c.mode.drawGUI = false;
        c.mode.saveXLS = true;
        c.mode.saveFIG = true;
        c.mode.generateFIG = false;
        c.mode.noEmpty = true;
        c.mode.noOverwrite = false;
        c.mode.exitDone = false;
        c.mode.drawLegend = false;
        c.mode.drawFIG = c.mode.saveFIG || c.mode.generateFIG;
        c.mode.collectEvents = false;
        c.mode.saveRaster = true;
        c.mode.saveRasterSum = true;
        c.mode.saveExcelSum = true;
        
        c.buttonDownXY = [];
        c.currXLimit = [];
        c.maxXLimit = [];
%        c.currYLimit = [];
        c.currEvents = [];
        c.allEvents = [];
        c.bulkEvents = [];
        
        c.currPos = [];
        
        guidata(common, c);        
    end
    function setCommon(varargin)
        c = guidata(common);
        for varI = 1 : 2 : length(varargin)
            dest = varargin{varI};
            value = varargin{varI+1};
            c.(dest) = value;
        end
        guidata(common, c);
    end
    function value = getCommon(dest)
        c = guidata(common);
        value = c.(dest);
    end
    function value = getParam(dest)
        c = getCommon('param');
        value = c.(dest);
    end
    function [dataInfo, input] = loadFileInfo(inputFile, inputDir)
        if nargin < 2
            inputDir = getCommon('inputDir');
        end
        
        input = strcat(inputDir, '/', inputFile);
        
        fprintf('loading file info ''%s''...\n', input);
        try
            dataInfo = myHdf5ToStruct(input, false);
            dataInfo.fullpath = input;
            dataInfo.filename = inputFile;
        catch
            dataInfo = [];
            fprintf('Error: Unable to read ''%s''\n', input);
            return;
        end
    end
    function [channel, timeInfo, recordingTime] = loadFileChannel(input, fieldPath, containerI)
        const = getCommon('const');
        currFileName = getCommon('currFileName');
        if strcmp(input, currFileName)
            timeInfo = getCommon('currFileTimeData');
        else
            fprintf('loading file time data''%s''...\n', input);
            try
                timeInfo = hdf5read(input, const.datetimePath);
                setCommon('currFileName', input, 'currFileTimeData', timeInfo);
            catch
                timeInfo = [];
                fprintf('Error: Unable to read time data from ''%s''\n', input);
            end
        end
        currChannelPath = getCommon('currChannelPath');
        if strcmp(input, currFileName) && strcmp(fieldPath, currChannelPath)
            channel = getCommon('currChannelData');
        else
            try
                channel = hdf5read(input, fieldPath);
                setCommon('currChannelPath', fieldPath, 'currChannelData', channel);
            catch
                channel = [];
                fprintf('Error: Unable to read channel data ''%s'' from ''%s''\n', fieldPath, input);
                return;
            end
        end
        recordingTimeStruct = hdf5read(input, '/', 'datetime');
        recordingTime = datenum(recordingTimeStruct.data, 'mm-dd-yy HH:MM:SS');
        
        lenRatio = length(timeInfo) / length(channel);
        if round(lenRatio) ~= lenRatio
            error('Error in file format: length of sample_t not a multiple of channel lengths\n');
        end
        timeInfo = timeInfo(containerI : lenRatio : end);
        
        % normalize timeInfo
        timeInfo = timeInfo - timeInfo(1);

        setCommon('currTimeInfo', timeInfo);
    end

    function dirList = getAllDirs(inputFile)
        allList = dir(inputFile);
        isdir = cat(1, allList.isdir);
        dirList = allList(isdir);
    end
    function [fileList, inputDir, dirList, b_selectFiles] = getFileList(inputFile)
        switch exist(inputFile)
            case 2
                [inputDir, fileName, fileExt] = fileparts(inputFile);
                fileList.name = strcat(fileName, fileExt);
                if nargout >= 3
                    dirList.name = '.';
                    dirlist.isdir = true;
                    if nargout >= 4
                        b_selectFiles = true;
                    end
                end
            case 7
                inputDir = inputFile;
                fileList = dir(strcat(inputFile, '\\*.hdf5'));
                isdir = cat(1, fileList.isdir);
                fileList = fileList(~isdir);
                
                if nargout >= 3
                    dirList = getAllDirs(inputFile);
                    dirList(strcmp(cat(1, {dirList.name}), '.')) = [];
                    if nargout >= 4
                        b_selectFiles = false;
                    end
                end
            otherwise
                inputDir = fileparts(inputFile);
                allList = dir(inputFile);
                isdir = cat(1, allList.isdir);
                fileList = allList(~isdir);
                if nargout >= 3
                    dirList = allList(isdir);
                    dotI = strcmp(cat(1, {dirList.name}), '.');
                    if sum(dotI) == 0
                        allDirList = getAllDirs(inputDir);
                        dotI = strcmp(cat(1, {allDirList.name}), '.');
                        dirList = [allDirList(dotI) dirList];
                    end
                    if nargout >= 4
                        b_selectFiles = true;
                    end                    
                end
        end    
        if isempty(inputDir)
            inputDir = '.';
        end        
    end
    function [allChannelList, allFileList] = loadAllFileInfos(fileList)
        allChannelList = [];
        allFileList = [];
        for fi = 1 : length(fileList)
            allFileList(fi).ID = fi;
            allFileList(fi).name = fileList(fi).name;
                        
            currInfo = loadFileInfo(fileList(fi).name);
            currList = genListFromInfo(currInfo, allFileList(fi).ID);
            allChannelList = [allChannelList currList];
            
            allFileList(fi).num = length(currList);
        end
    end
    function doBatch(inputFile, p_close)        
        setMode('batchOnly');
        setMode('generateFIG', ~p_close);

        [fileList, inputDir] = getFileList(inputFile);
        setCommon('inputDir', inputDir);
        
        listAll = loadAllFileInfos(fileList);
        setCommon('allFileInfos', listAll);        
        
        listFiltered = applyFilters(listAll);
        
        runBatch(listFiltered);
    end

    function ret = genListFromInfo(dataInfo, ID)                
        ret = [];
        retI = 1;
        containerI = 0;
        
        const = getCommon('const');
        
        if isfield(dataInfo, const.datetimeName)
            containerNames = fieldnames(dataInfo);
            for c = 1 : length(containerNames)
                if isempty(regexp(containerNames{c}, const.container, 'once' ))
                    continue;
                end
%                fprintf('processing container ''%s''\n', containerNames{c});
                containerI = containerI + 1;
                channelNames = fieldnames(dataInfo.(containerNames{c}));
                for i = 1 : length(channelNames)                    
%                    fprintf('processing data field ''%s''...\n', channelNames{i}); 
%                    doWorkChannel(data, inputFile, containerNames{c}, channelNames{i});
                    ret(retI).ID = ID;
                    ret(retI).fullpath = dataInfo.fullpath;
                    ret(retI).filename = dataInfo.filename;
                    ret(retI).fieldPath = dataInfo.(containerNames{c}).(channelNames{i});
                    ret(retI).filterValid = true;
                    ret(retI).valid = true;
                    ret(retI).containerI = containerI;
                    
                    retI = retI + 1;
                end
            end
        else
            fprintf('field ''%s'' not found in hdf5 data structure\n', const.datetimeName);
        end
    end
    function allChannelList = applyFilters(allChannelList)
        p_filter = getCommon('p_filter');
        
        for i = 1 : length(allChannelList)
            [containerName, rest] = strtok(allChannelList(i).fieldPath, '/');            
            channelName = strtok(rest, '/');
            
            if ~isempty(p_filter.container) && ~strcmp(p_filter.container, containerName)
                allChannelList(i).filterValid = false;
                continue;
            end            
            if ~isempty(p_filter.channel) && ~strcmp(p_filter.channel, channelName)
                allChannelList(i).filterValid = false;
                continue;
            end
            allChannelList(i).filterValid = true;
        end
    end
    function done = runBatch(allChannelList, xlimits)
        if nargin < 2
            xlimits = [NaN NaN]; % default limits
        end
        todoLen = length(allChannelList);
        drawWaitBar = todoLen > 1;
        if drawWaitBar
            waitbarHandle = waitbar(0, '', 'Name', 'Bulk Processing...', ...
                                           'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
            setappdata(waitbarHandle,'canceling',0);
        end
        done = {};
        for i = 1 : todoLen
            [containerName, rest] = strtok(allChannelList(i).fieldPath, '/');            
            channelName = strtok(rest, '/');
            
            if ~allChannelList(i).filterValid
                continue;
            end
            
            fprintf('processing container ''%s'', data field ''%s'' in file ''%s''...\n', containerName, channelName, allChannelList(i).filename);
            [channel, timeInfo, recordingTime] = loadFileChannel(allChannelList(i).fullpath, allChannelList(i).fieldPath, allChannelList(i).containerI);
            
            doWorkChannel(channel, timeInfo, recordingTime, allChannelList(i).filename, containerName, channelName, xlimits);
            done{i} = sprintf('%s.%s.%s', allChannelList(i).filename, containerName, channelName);
            if drawWaitBar
%                figure(getCommon('hfig'));
                waitbar(i/todoLen, waitbarHandle, sprintf('%d / %d', i, todoLen));
                if getappdata(waitbarHandle,'canceling')
                    break
                end                
            end
        end    
        if drawWaitBar
            delete(waitbarHandle);
        end
    end

    function ret = myGetEventList(reto, fh, timeInfo, it)
        if nargin == 2
            interval = fh;
        end
        if nargin < 3
            eventData = reto.data;
            reto = eventData.reto; % added original data
            fh = eventData.fh;
            timeInfo = eventData.timeInfo;
            it = eventData.it;
        end
        if nargin == 2
            interval(2) = min(interval(2), length(fh));
            reto = reto(interval(1) : min(length(reto), interval(2)+1)); % added original data
            fh = fh(interval(1) : min(length(fh), interval(2)+1));
            timeInfo = timeInfo(interval(1) : min(length(timeInfo), interval(2)+1));
            it = it(interval(1) : min(length(it), interval(2)+1));
        end        
        
        dummyEvent = struct('startIdx', [], 'endIdx', [], 'startTime', [], 'stopTime', [], ... 
                            'startVol', [], 'stopVol', [], 'startVolS', [], 'stopVolS', [], ...
                            'durationIdx', [], 'durationTime', [], 'volume', [], 'volumeS', [], ...
                            'interMeal', [], 'speed', [], 'speedS', []);
        
        dit = [it; 0] - [0; it];
%        len = length(timeInfo);
        startIdx = find(dit == 1); %startIdx(startIdx > len) = len;
        endIdx = find(dit == -1) - 1;% endIdx(endIdx > len) = len;
        
        numEvents = length(startIdx);
        ret = repmat(dummyEvent, 1, numEvents);
        
        for e = 1 : numEvents
            ret(e).startIdx = startIdx(e);
            ret(e).endIdx = endIdx(e);
            
            ret(e).startTime = timeInfo(startIdx(e));
            ret(e).stopTime  = timeInfo(min(length(timeInfo), endIdx(e) + 1));
            
            ret(e).startVol = reto(startIdx(e));
            ret(e).stopVol  = reto(min(length(reto), endIdx(e) + 1));

            ret(e).startVolS = fh(startIdx(e));
            ret(e).stopVolS  = fh(min(length(fh), endIdx(e) + 1));            
            
%             if e == 1
%                 ret(e).interMeal = NaN;
%             else
%                 ret(e).interMeal = ret(e).startTime - ret(e-1).stopTime;
%             end
            ret(e).interMeal = compInterMeal(ret, e);
            
            ret(e) = compSingleEventStats(ret(e));
        end
    end
    function ret = compInterMeal(eventList, idx)
            if idx == 1
                ret = NaN;
            else
                ret = eventList(idx).startTime - eventList(idx-1).stopTime;
            end
    end

    function event = compSingleEventStats(event)
            event.durationIdx = event.endIdx - event.startIdx + 1;
            event.durationTime = event.stopTime - event.startTime;
            
            event.volume = event.stopVol - event.startVol;
            event.speed = abs(event.volume) / event.durationTime;
            
            event.volumeS = event.stopVolS - event.startVolS; % fh(startIdx(e)) - fh(endIdx(e));                    
            event.speedS = abs(event.volumeS) / event.durationTime;
    end
    function eventList = filterEvents(eventList)            
        eventList = filterVolumeCutOff(eventList);
        eventList = mergeBigMeals(eventList);

        function eventList = filterVolumeCutOff(eventList)
%                 allVolumes = cat(1, eventList.volumeS);
%                 minEventVolume = getParam('minEventVolume');
%                 invalid = allVolumes < minEventVolume;
%            invalidIdx = find(checkStatThresh(eventList, 'volumeS', getParam('minEventVolume')));
            invalidIdx = find(checkStatThresh(eventList, 'durationIdx', getParam('minEventDuration')) | ...
                              checkStatThresh(eventList, 'volumeS',     getParam('minEventVolume')));
            % eliminate events and update interMeal values for successor events
            for e = 1 : length(invalidIdx)
                eventList(invalidIdx(e)) = [];
                invalidIdx = invalidIdx - 1;
                
                if invalidIdx(e) < length(eventList)
                    eventList(invalidIdx(e)+1).interMeal = compInterMeal(eventList, invalidIdx(e)+1);
                end                
            end
        end
        function eventList = mergeBigMeals(eventList)
            invalidIdx = find(checkStatThresh(eventList, 'interMeal', getParam('minEventInterMeal')));
            % merge events with previous event and update SingleEventStats
            for e = 1 : length(invalidIdx)
                eventList(invalidIdx(e)-1).endIdx = eventList(invalidIdx(e)).endIdx;
                eventList(invalidIdx(e)-1).stopTime = eventList(invalidIdx(e)).stopTime;
                eventList(invalidIdx(e)-1).stopVol = eventList(invalidIdx(e)).stopVol;

                eventList(invalidIdx(e)-1) = compSingleEventStats(eventList(invalidIdx(e)-1));

                eventList(invalidIdx(e)) = [];
                invalidIdx = invalidIdx - 1;
            end
        end
        function ret = checkStatThresh(eventList, statName, statThresh)
            if isempty(eventList)
                ret = [];
            else
                allStats = cat(1, eventList.(statName));
                ret = abs(allStats) < statThresh;                
            end
        end
    end
    function eventList = myGetFilteredEventList(events, interval)
        rawEventList = myGetEventList(events, interval);
        eventList = filterEvents(rawEventList);
    end


%     function ret = getLatencyFromEventSummary(currEvents)
%         if isempty(currEvents) || isempty(currEvents.eventList)
%             ret = NaN;
%         else
%             ret = currEvents.eventList(1).startTime;
%         end
%     end
    function addBulkEvent(currEvents)
        if isempty(currEvents)
            return;
        end
        bulkEvents = getCommon('bulkEvents');
        if isempty(bulkEvents)
            bulkEvents.eventList = currEvents.eventList;
%            bulkEvents.latencies = getLatencyFromEventSummary(currEvents);
            bulkEvents.stats = currEvents.stats;
            bulkEvents.info = currEvents.info;
        else
            bulkEvents.eventList = [bulkEvents.eventList currEvents.eventList];
%            bulkEvents.latencies = [bulkEvents.latencies getLatencyFromEventSummary(currEvents)];
            bulkEvents.stats = [bulkEvents.stats currEvents.stats];
            bulkEvents.info = [bulkEvents.info currEvents.info];
        end
        setCommon('bulkEvents', bulkEvents);
    end
    function adjustYLimit(xlimits, data, handle, which)
        if nargin < 3
            handle = gca;
        end
        if isempty(data)
            if nargin < 4
                return;
            else
                switch which
                    case 1
                        data = getCommon('currChannelData');
                    case 2
                        allEvents = getCommon('allEvents');
                        data = allEvents.data.fh;
                    case 3
                        return;
                end
            end
        end
        xlimits = assignBoundedValue(xlimits, [1 length(data)]);
        ylimits = [min(data(xlimits(1) : xlimits(2))) max(data(xlimits(1) : xlimits(2)))];
        if ylimits(1) == ylimits(2)
            yLimitOffset = ylimits(1) * 0.10;
        else
            yLimitOffset = (ylimits(2) - ylimits(1)) * 0.10; % +/- 10 percent
        end
        ylimits = [ylimits(1)-yLimitOffset ylimits(2)+yLimitOffset];
        set(handle, 'YLim', sort(ylimits));
    end

    function doWorkChannel(channel, timeInfo, recordingTime, inputFile, containerI, channelI, xlimits)
        
        function ret = myInterpolateIdx(data, t, idx)
            ret = data;
            if sum(idx) > 0
                xi = t(idx);
                x = t; x(idx) = [];
                Y = data; Y(idx) = [];
                yi = interp1(x, Y, xi);
                ret(idx) = yi;
            end
        end
        function ret = myInterpolateNaNs(data, t)
            ret = myInterpolateIdx(data, t, isnan(data));
        end

        function [ret, idx, rangeIdx] = myClean(data, t)
            s = size(data);
            if s(1) == 1
                data = data';
            end
            
            data(data == -1) = NaN;
            nans = isnan(data);
            if sum(nans) > 0 % ignore NaNs at beginning or end, interpolate all other NaNs
%                 startI = 1;
%                 len = length(data);
%                 while(nans(startI) && startI < len)
%                     startI = startI + 1;
%                 end
%                 endI = len;
%                 while(nans(endI) && endI > 1)
%                     endI = endI - 1;
%                 end            
                [startI, endI] = myCleanHT(data, nans);
                rangeIdx = startI : endI;
                ret = data(startI : endI);
                t = t(startI : endI);
                
                ret = myInterpolateNaNs(ret, t);
                idx = find(nans(startI : endI));            
            else % nothing else to do when no NaNs
                ret = data;
                idx = [];
                rangeIdx = 1 : length(data);
            end
        end
        function [startI, endI] = myCleanHT(data, nans)
            if nargin < 2
                nans = isnan(data);
            end
            startI = 1;
            len = length(data);
            while startI <= len && nans(startI)
                startI = startI + 1;
            end
            endI = len;
            while endI >= startI && nans(endI)
                endI = endI - 1;
            end            
        end
        
        function fret = myDenoise(data)
%            lev = 5;
            lev = 3; % new value 20140211
            fret = data;
            [startI, endI] = myCleanHT(data);
            fret(startI:endI) = wden(data(startI:endI),'sqtwolog','s','sln',lev,'sym8');
        end

        function ret = estimateStd(data)
            ret = nanstd(data);
            return;
            
%             len = length(data);
%             step = len / 10;
%             for i = 1 : 10
%                 v(i) = nanstd(data(1+(i-1)*step : i*step));
%             end
%             ret = median(v);
        end
        function [ret, v] = myStdfacts(data)
            mu = nanmean(data);
%            v = nanstd(data);
            v = estimateStd(data);
%            v = nanstd(data) / (sum(~isnan(data)));

            ret = (data - mu) ./ v;
        end
        function [it, sf, v] = myStdfactsIdx(data, devfakt)
            [sf, v] = myStdfacts(data);
            it = abs(sf) > devfakt; % it = (dfh2 > mu + v*devfakt | dfh2 < mu - v*devfakt);
        end
        function [ret, idx] = myRemoveOutliers(data, t, devfakt)
            if nargin < 3
                devfakt = 3;
            end

        % detect outliers    
%            o = myStdfactsIdx(data, devfakt) | isnan(data);
            o = myStdfactsIdx(data, devfakt) | isnan(data);
            idx = find(o);

        % interpolate values    
            ret = myInterpolateIdx(data, t, idx);
        end
        
        function [fh, it, stdRatio, patchedData] = myEvents(data, devfakt)
            if nargin < 2
                devfakt = 1; % thresh of std deviation
            end

%            v_thresh = 0.01; % find optimal threshold here
            v_thresh = 0.01; % find optimal threshold here % new value 20140211
            emptyV = 0.22; % highest observed stdev-value when no drinking event: 0.196, lowest observed whith drinking event: 0.25

            fh = data;
            if isempty(data)
                it = false(size(data));
                stdRatio = NaN(size(data));
                patchedData = data;
                return;
            end

            % detect drinking events
            dfh = diff(fh);
            cumIt = false(size(dfh));
            allV = [];
            
            [it, stdRatio, v] = myStdfactsIdx(dfh, devfakt);
            
            while(v > emptyV)
                allV = [allV v];
                cumIt = cumIt | it;
                dfh(it) = NaN;
                [it, ~, v] = myStdfactsIdx(dfh, devfakt);
            end
            
            it = cumIt;
            it(it) = stdRatio(it) < 0; % remove positive slopes
            
%            it = [false ; it];
%            stdRatio = [0 ; stdRatio];
            it = [it; false];
            stdRatio = [stdRatio; 0];

            if v == 0
                fprintf('no variance - concluding that there is no signal in this chamber\n');
                it(:) = false;
                stdRatio(:) = 0;
            elseif v < v_thresh % if too little variance -> no events at all; find optimal thresh, currently 0.01
                fprintf('variance too low - assuming there are no events at all\n');
                it(:) = false;
                stdRatio(:) = 0;
            end

            % compute drinking-corrected data
            if nargout >= 4
                itp = find(it);
                patchedData = fh;
                mu = mean(dfh);
                for iti = 1 : length(itp)
                    patchedData(itp(iti) : end) = patchedData(itp(iti) : end) - dfh(itp(iti)) + mu;
                end
            end
        end
        
        function [h, figPanel] = myInitFig(titleText, mode)
            if mode.generateFIG
                h = figure('Name', titleText,'NumberTitle','off');
            else
                h = figure('Name',titleText,'NumberTitle','off', 'visible', 'off');
            end
            figPanel(1) = subplot(3,1,1);
            figPanel(2) = subplot(3,1,2);
            figPanel(3) = subplot(3,1,3);
        end
        function xlimits = myPlotOutliers(subPanel, data, timeInfo, rangeIdx, outlierIdx, dataWithoutOutliers, xlimits, drawLegend)
            if isstruct(subPanel)
                handle = subPanel.handle;
            else
                handle = subPanel;
            end
            subplot(handle);
            noDraw = isempty(dataWithoutOutliers);
            if noDraw
                emptyPanel(); % clear Panel
            else
                axis on;
%                plot(data, 'k');
                plot(timeInfo, data, 'k');
                hold on;
%                plot(dataWithoutOutliers, 'b');
                plot(timeInfo(rangeIdx), dataWithoutOutliers, 'b');                
                fo = NaN(size(data));
                fo(outlierIdx) = data(outlierIdx);                
%                plot(fo, 'rx');                
                plot(timeInfo, fo, 'rx');                
                fo(outlierIdx) = dataWithoutOutliers(outlierIdx);
%                plot(fo, 'bx');
                plot(timeInfo, fo,'bx');                
                hold off;
                legend('original data [liquid level / time ]', 'outlier free data', 'outlier value', 'substituted value');
                if ~drawLegend
                    legend('off');                
                end
            end            
            
%                if nargin < 5 || isempty(xlimits)
%                    xlimits = get(gca, 'XLim');
%                else
                    set(gca, 'XLim', timeInfo(xlimits));
%                end
                
            if ~noDraw
                adjustYLimit(xlimits, data);
                set(gca, 'TickLength', [0.005 0.0250]); % 0.005: 2D (default: 0.010), 0.0250: 3D (default: 0.0250)
                if isstruct(subPanel)
                    labelInsideAxes(handle, subPanel.axes, subPanel.axesLabel, subPanel.offset);
                end
            end
        end
        function xlimits = myPlotEvents(subPanel, currEvents, timeInfo, xlimits, drawLegend)
            if isstruct(subPanel)
                handle = subPanel.handle;
            else
                handle = subPanel;
            end
            subplot(handle);
            noDraw = isempty(currEvents) || isempty(currEvents.data.fh);
            if noDraw
                emptyPanel(); % clear panel
%                msgbox('The channel doesn not contain valid data', 'Invalid Channel Data');
            else        
                axis on;
%                dummyX = 1 : length(currEvents.data.fh);                
                currTimeInfo = currEvents.data.timeInfo(1 : length(currEvents.data.fh));
                
                plot(currTimeInfo, currEvents.data.fh, 'b')
                    hold on;
                for e = 1 : length(currEvents.eventList)
                mk = currEvents.data.filteredIt*1.0;
                mk(:) = NaN;
                    eRange = currEvents.eventList(e).startIdx : currEvents.eventList(e).endIdx + 1;
                    mk(eRange) = currEvents.data.fh(eRange);
                    eHandle = plot(currTimeInfo, mk, 'r');
                    if e > 1 % add only first entry to legend
                        set(get(get(eHandle,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
                    end
                end
                    hold off;
%                mk(currEvents.data.filteredIt) = currEvents.data.fh(currEvents.data.filteredIt);

%                plot(currTimeInfo, currEvents.data.fh, 'b', currTimeInfo, mk, 'r')
                legend('denoised data [liquid level / time]', 'detected events');
                if ~drawLegend
                    legend('off');                
                end
            end
                
%                if nargin < 3 || isempty(xlimits)
%                    xlimits = get(gca, 'XLim');
%                else
                    set(gca, 'XLim', timeInfo(xlimits));
%                end
                
            if ~noDraw
                adjustYLimit(xlimits, currEvents.data.fh);            

                set(gca, 'TickLength', [0.005 0.0250]); % 0.005: 2D (default: 0.010), 0.0250: 3D (default: 0.0250)            
                if isstruct(subPanel)
                    labelInsideAxes(handle, subPanel.axes, subPanel.axesLabel, subPanel.offset);
                end
            end
        end
        function xlimits = myPlotDebug(subPanel, stdRatio, timeInfo, devfakt, xlimits, drawLegend)
            if isstruct(subPanel)
                handle = subPanel.handle;
            else
                handle = subPanel;
            end
            subplot(handle);
            noDraw = isempty(fh);
            if noDraw
                emptyPanel(); % clear panel
                axis on;
                text(mean(xlimits), mean(get(gca, 'YLim')), 'no data', 'fontsize', 25, 'horizontalalignment','center', 'verticalalignment','middle');
                set(gca, 'xticklabel', [], 'yticklabel', []);
            else        
                axis on;
                currTimeInfo = timeInfo(1 : length(stdRatio));
                plot(currTimeInfo, stdRatio, 'm');
                if ~isempty(stdRatio)
                    hold on;                    
                    orientationMarkHandles = plot(currTimeInfo, ones(size(stdRatio)),    'k', currTimeInfo, ones(size(stdRatio))*2,  'k', currTimeInfo, ones(size(stdRatio))*3,  'k',...
                                                  currTimeInfo, ones(size(stdRatio))*-1, 'k', currTimeInfo, ones(size(stdRatio))*-2, 'k', currTimeInfo, ones(size(stdRatio))*-3, 'k');
                    compactLegend(orientationMarkHandles);
                    
                    currThreshHandles = plot(currTimeInfo, ones(size(stdRatio))*devfakt.event,  'r', ...
                                             currTimeInfo, ones(size(stdRatio))*-devfakt.event, 'r');
                    compactLegend(currThreshHandles);
                    
                    hold off;
                end
                set(gca, 'YLim', [-4 4]);
                legend('standard deviation [\sigma / time]', 'orientation marks', 'current threshold');
                if ~drawLegend
                    legend('off');
%                else
%                    legend('on');
                end
            end
            
%                if nargin < 4 || isempty(xlimits)
%                    xlimits = get(gca, 'XLim');
%                else
                    set(gca, 'XLim', timeInfo(xlimits));
%                end
                
            if ~noDraw
                set(gca, 'TickLength', [0.005 0.0250]); % 0.005: 2D (default: 0.010), 0.0250: 3D (default: 0.0250)            
                if isstruct(subPanel)
                    labelInsideAxes(handle, subPanel.axes, subPanel.axesLabel, subPanel.offset);
                end
            end
            
            function compactLegend(handles)
                handleGroup = hggroup;
                set(handles,'Parent',handleGroup)
                set(get(get(handleGroup,'Annotation'),'LegendInformation'), 'IconDisplayStyle','on');
            end
        end
        function mySaveFIG(figureHandle, figureFileName, mode, emptyChannel)
            if emptyChannel && mode.noEmpty
                fprintf('figure not saved: empty data channel and ''noEmpty'' is on, file ''%s'' not written\n', figureFileName);
            elseif mode.noOverwrite && exist(figureFileName, 'file')
                fprintf('figure not saved: file already exists and option ''noOverwrite'' is on, file: ''%s''\n', figureFileName);
            else
                %print(figureHandle, '-dtiff', figureFileName);
                print(figureHandle, '-depsc', figureFileName);
                fprintf('figure saved in ''%s''\n', figureFileName);               
            end
        end
        function ret = myEventStats(reto, fh, timeInfo, it, rangeIdx, xlimits, inputDir, fileText, recordingTime)   
            if isempty(fh)
                ret = [];
                return;
            end
            ret.info.inputDir = inputDir;
            ret.info.fileText = fileText;
            ret.info.timeInfo = timeInfo;
            ret.info.recordingTime = recordingTime;
            
            ret.data.reto = reto;
            ret.data.fh = fh; 
            ret.data.timeInfo = timeInfo; 
            ret.data.it = it;
            currRange = max(xlimits(1), rangeIdx(1)) : min(xlimits(2), rangeIdx(end));
%             if isempty(fh)
%                 ret.rawEventList = [];
%             else
                ret.rawEventList = myGetEventList(reto(currRange), fh(currRange), timeInfo(currRange), it(currRange));
            %end
            ret.eventList = filterEvents(ret.rawEventList);
            
            ret.data.filteredIt = computeItFromEvents(ret.data.it, ret.eventList);
            
            ret.stats = myComputeEventStats(ret.eventList, timeInfo, currRange, fh, ret.data.filteredIt);
            
            function ret = computeItFromEvents(it, eventList)
                ret = false(size(it));
                for i = 1 : length(eventList)
                    ret(eventList(i).startIdx : eventList(i).endIdx) = true;
                end
            end
        end
        function ret = myComputeEventStats(eventList, timeInfo, currRange, fh, it)
%            ret.latency = getLatencyFromEventSummary(ret) - timeInfo(1);

            binSec = getCommon('rasterSumBinSec');
            rest = mod(timeInfo(end), binSec);
            numBins = floor(timeInfo(end) / 10) + (rest > 0) * 1;
            binVol = NaN(numBins, 1);

            if isempty(eventList)
                ret.latency = NaN;
                
                ret.avgInterMeal = NaN;
                ret.avgSpeed = NaN;
                
                ret.halfDrinkAt = NaN;
                
                ret.totalVol_1_3 = NaN;
                ret.totalVol_2_3 = NaN;
                ret.totalVol_3_3 = NaN;
            else
                ret.latency = eventList(1).startTime - timeInfo(1);
                
                ret.avgInterMeal = nanmean(cat(1, eventList.interMeal));
                ret.avgSpeed = nanmean(cat(1, eventList.speed));                
                
                % determine time when half the total volume is drunk
                dfh = diff(fh);
                dfh(~it) = 0;     
                ret.halfDrinkAt = timeInfo(find(cumsum(dfh) < sum(cat(1, eventList.volumeS))/2, 1)); 
   
                % determine food intake for first / middle / last third
                cl = length(currRange);
                totalVol_thirds = getTotalVol_perSegment(eventList, [floor(cl/3) floor(cl/3)*2 cl]);
                ret.totalVol_1_3 = totalVol_thirds(1);
                ret.totalVol_2_3 = totalVol_thirds(2);
                ret.totalVol_3_3 = totalVol_thirds(3);

                
                % determine volume per x-sec bin
                currIdx = 0;
                for i = 1 : numBins
                    startIdx = currIdx+1;
                    currIdx = find(timeInfo > i * binSec, 1) - 1;
                    if isempty(currIdx)
                        currIdx = length(timeInfo);
                    end
                    binVol(i) = sum(dfh(startIdx : currIdx));
                end
            end
            ret.binVol = binVol;            
            
            ret.total = length(eventList);
            ret.length = length(currRange);
            freqPrec = getParam('freqPrec'); % round freqPrec digits after comma
            ret.freq = 1 / median(diff(timeInfo(currRange)));
            ret.freq = myRound(ret.freq, freqPrec);
            %ret.freq = round(ret.freq * 10^freqPrec) / 10.0^freqPrec;                                                
            
            function ret = getTotalVol_perSegment(eventList, segments)
                ret = zeros(length(segments), 1);
                for ei = 1 : length(eventList)
                    start_si = find(eventList(ei).startIdx < segments, 1);
                    end_si   = find(eventList(ei).endIdx < segments, 1);
                    vol = abs(eventList(ei).volumeS);
                    if start_si == end_si
                        ret(start_si) = ret(start_si) + vol;
                    else
                        len_segA = segments(start_si) - eventList(ei).startIdx + 1;
                        len_segB = eventList(ei).endIdx - segments(start_si);
                        len_both = len_segA + len_segB;
                        ret(start_si) = ret(start_si) + vol * len_segA / len_both;
                        ret(end_si) = ret(end_si) + vol * len_segB / len_both;
                    end
                end
            end
        end

        function mySaveXLS(csvFileName, currEvents, recordingTime, mode)
            if isempty(currEvents) && mode.noEmpty
                fprintf('stats not written: empty data channel and ''noEmpty'' is on, file ''%s'' not written\n', csvFileName);
            elseif mode.noOverwrite && exist(csvFileName, 'file')
                fprintf('stats not written: file already exists and option ''noOverwrite'' is on, file: ''%s''\n', csvFileName);
            else
                csvHandle = fopen(csvFileName, 'w');
                if csvHandle < 0
                    msgbox(sprintf('Unable to write csv-file %s', csvFileName),'Error');
                else
                
    %                fprintf(csvHandle, 'startIdx [#], endIdx [#], durationIdx [#], startTime [s], stopTime [s], durationTime [s], volume [nl], speed [nl/s], volumeS [nl], speedS [nl/s], interMeal [s], recordedTime [HH-MM-SS  mm/dd/yyyy]\n');
                    fprintf(csvHandle, 'startIdx [#], endIdx [#], durationIdx [#], startTime [s], stopTime [s], durationTime [s], volume [nl], speed [nl/s], interMeal [s], recordedTime [HH-MM-SS  mm/dd/yyyy]\n');

                    if isempty(currEvents)
                        fprintf(csvHandle, 'no data');
                    else
                        events = currEvents.eventList;                                
                        for e = 1 : length(events)                    
    %                        fprintf(csvHandle, '%d, %d, %d, %f, %f, %f, %f, %f, %f, %f, %s, %s\n', events(e).startIdx, events(e).endIdx, events(e).durationIdx, events(e).startTime, events(e).stopTime, events(e).durationTime, abs(events(e).volume), events(e).speed, abs(events(e).volumeS), events(e).speedS, myNum2str(events(e).interMeal), datestr(recordingTime + datenum(0,0,0,0,0,double(events(e).startTime)), 'HH-MM-SS  mm/dd/yyyy'));
                            fprintf(csvHandle, '%d, %d, %d, %f, %f, %f, %f, %f, %s, %s\n', ...
                                events(e).startIdx, events(e).endIdx, events(e).durationIdx, ...
                                events(e).startTime, events(e).stopTime, events(e).durationTime, ...
                                abs(events(e).volumeS), events(e).speedS, myNum2str(events(e).interMeal), ...
                                datestr(recordingTime + datenum(0,0,0,0,0,double(events(e).startTime)), 'HH-MM-SS  mm/dd/yyyy'));
                        end

                        fprintf(csvHandle, '\n\ntotal events [#], latency [s], avgInterMeal [s], avgSpeed [s], freq [1/s], data points [#], halfDrinkAt [s], volume_1_3 [nl], volume_2_3 [nl], volume_3_3 [nl]\n');
                        if ~isempty(currEvents) && ~isempty(currEvents.stats)
                            stats = currEvents.stats;                
                            fprintf(csvHandle, '%d, %s, %f, %f, %f, %d, %f, %f, %f, %f\n', stats.total, myNum2str(stats.latency), stats.avgInterMeal, stats.avgSpeed, stats.freq, stats.length, stats.halfDrinkAt, stats.totalVol_1_3, stats.totalVol_2_3, stats.totalVol_3_3);
                        end
                    end

                    fclose(csvHandle);
                    fprintf('stats written in ''%s''\n', csvFileName);     
                end
            end
        end
        
        
        
        
        
%% body doSingleWork
%        const = getCommon('const');
        devfakt = getCommon('devfakt');
        mode = getCommon('mode');

%        channel = data.(containerI).(channelI);
%        timeInfo = data.(const.datetimeName);
        
        fprintf('data length: %d entries\n', length(channel));        
        
%% remove outliers
        [cleanChannel, cleanIdx, rangeIdx] = myClean(channel, timeInfo);
        [reto, outlierIdx] = myRemoveOutliers(cleanChannel, timeInfo, devfakt.outlier);
        
        inputDir = getCommon('inputDir');
        fileText = sprintf('''%s''.%s.%s', inputFile, containerI, channelI);
        if mode.drawFIG || mode.drawGUI
            titleText = sprintf('Data from file: ''%s/%s''.%s.%s', inputDir, inputFile, containerI, channelI);
        end
        xlimitDefault = [1 length(cleanChannel)];
        xlimits(isnan(xlimits)) = xlimitDefault(isnan(xlimits));        
        if mode.drawFIG
            [h, figPanel] = myInitFig([getCommon('titlePrefix') titleText], mode);
            myPlotOutliers(figPanel(1), channel, timeInfo, rangeIdx, [cleanIdx outlierIdx], reto, xlimits, mode.drawLegend);
        end

%% denoise
        fret = myDenoise(reto);

%% detect drinking events
        [fh, it, stdRatio] = myEvents(fret, devfakt.event);
        
        currEvents = myEventStats(reto, fh, timeInfo, it, rangeIdx, xlimits, inputDir, fileText, recordingTime);
        if mode.collectEvents
            addBulkEvent(currEvents);
        end        

        if mode.drawFIG
            myPlotEvents(figPanel(2), currEvents, timeInfo, xlimits, mode.drawLegend);
            myPlotDebug(figPanel(3), stdRatio, timeInfo, devfakt, xlimits, mode.drawLegend);
        end
                
%% save figure to file
        if mode.saveFIG
            outputDir = getCommon('outputDir');
            figureFileName = genOutFileName('.eps', xlimits, all(xlimits == xlimitDefault));
            mySaveFIG(h, figureFileName, mode, isempty(fh));
        end
        if mode.drawFIG && (~mode.generateFIG || (mode.noEmpty && isempty(fh)))
            close(h);
        end        
        
%% draw GUI-figs
        if mode.drawGUI
            hfig = getCommon('hfig');
            set(hfig, 'Name', [getCommon('titlePrefix') titleText]);
            
%            panel = getCommon('panel');            
            xlimitsPanelData = getCommon('xlimitPanelsData');            
            xlimits = getCommon('currXLimit');
            if isempty(xlimits)
                xlimits = [1 length(channel)];
                setCommon('allEvents', currEvents);
                setCommon('maxXLimit', xlimits);
            end
            myPlotOutliers(xlimitsPanelData(1), channel, timeInfo, rangeIdx, [cleanIdx outlierIdx], reto, xlimits, mode.drawLegend);
            setCommon('currXLimit', xlimits);
            myPlotEvents(xlimitsPanelData(2), currEvents, timeInfo, xlimits, mode.drawLegend);
            myPlotDebug(xlimitsPanelData(3), stdRatio, timeInfo, devfakt, xlimits, mode.drawLegend);
            
            setCommon('currEvents', currEvents);
        end

%% create .csv file                
        if mode.saveXLS
%            outputDir = getCommon('outputDir');
%            xlsFileName = genOutFileName('.csv', outputDir, inputFile, containerI, channelI, xlimits, all(xlimits == xlimitDefault));
            xlsFileName = genOutFileName('.csv', xlimits, all(xlimits == xlimitDefault));
            mySaveXLS(xlsFileName, currEvents, recordingTime, mode);
        end        
        
        function ret = genOutFileName(extension, xlimits, includeRange)
            if includeRange
                ret = sprintf(['%s/%s.%s.%s' extension],        outputDir, inputFile, containerI, channelI);
            else
                ret = sprintf(['%s/%s.%s.%s_%dto%d' extension], outputDir, inputFile, containerI, channelI, xlimits(1), xlimits(2));
            end            
        end
    end   
function ret = myNum2str(value)
    if isnan(value)
        ret = '';
    else
        ret = sprintf('%f', value);
    end
end








% function buttonUpHandle(~, ~)
%     markerWidth = 5;
%     clickXY = get(gcbf,'CurrentPoint');
%     figpos  = get(hfig,'Position');
%     pp = 0;
%     if clickXY(1) > figpos(3)/2 &&  clickXY(2) > figpos(4)/2
%         pp = clickXY(1) - figpos(3)/2;
%     end
%     if clickXY(1) < figpos(3)/2 &&  clickXY(2) > figpos(4)/2
%         pp = clickXY(1);
%     end
%     if pp ~= 0
%         g = guidata(hfig);
%         g.pos = [g.pos pp];
% 
%         axes(videoPanel(1));
%         x(1) = line([pp pp],[0 figpos(4)/2],'Color','r');
%         f = find(g.contr(:, pp, g.curFrame));
%         x(2) = line([pp-markerWidth pp+1+markerWidth],[min(f(:)) min(f(:))],'Color','y');
%         x(3) = line([pp-markerWidth pp+1+markerWidth],[max(f(:)) max(f(:))],'Color','y');
% 
%         axes(videoPanel(2));
%         x(4) = line([pp pp],[0 figpos(4)/2],'Color','r');
%         [pw, freq] = calcStatisticForYdiff(g.ydiff,pp);
%         g.posHandler = [g.posHandler; x];
%         figure(hfig);
%         axes(resultPanel);
%         g.stat = plot(freq*15,pw);
%         grid on;
%         guidata(hfig, g);
%     end
% end


function KeypressFcn(~,~)
end
function ResizeFcn(~,~)
    setCommon('init_complete', false);
    
    hfig = getCommon('hfig');    
    gui = get(hfig, 'userdata');
    guiBackup = gui;
    
	[pos, gui] = calcGUIpostions(initGUIConst(hfig, false));
    
    gui.panel = guiBackup.panel;
    gui.handles = guiBackup.handles;
    gui.logo = guiBackup.logo;
    gui.pos = pos;    
    
    updateSubPanels(gui);
    gui.logo.pos = updateLogoPos(gui);
    
    allHandleNames = fieldnames(gui.handles);
    for hi = 1 : length(allHandleNames)
        set(gui.handles.(allHandleNames{hi}), 'Position', gui.pos.(allHandleNames{hi}));
    end
    
	set(hfig, 'userdata', gui);
    guidata(hfig, gui);
    
    setCommon('init_complete', true);
end
function CloseRequestFcn(~,~)
    g = guidata(common);
    
    if ~isempty(g) && isfield(g, 'htimer')
        stop(g.htimer); % Shut off timer if running
    end
    
    delete(gcbf);
    close(common);
end

function [panelID, columnID] = getHPanelByYPosition(y, gui)
%    panelID = max(0, ceil((y-figureHeight*(gui.videoHeightFactor+gui.infoHeightFactor))/(figureHeight*gui.panelHeightFactor)));
    if length(y) == 2
        x = y(1);
        y = y(2);
        columnID = find(x > unique(gui.panelOffsetX, 'rows'), 1, 'last');        
        if isempty(columnID)
            columnID = 1;
        end
    else
        columnID = 1;
    end
    panelID = find(y < gui.panelOffsetY(:, columnID), 1);
    if isempty(panelID)
        panelID = 1;
    end
end
function buttonID = getButtonID(f)
    button = get(f, 'SelectionType');
    if strcmp(button, 'normal')
        buttonID = 1;
    elseif strcmp(button, 'alt')
        buttonID = 2;
    elseif strcmp(button, 'extend')
        buttonID = 3;
    elseif strcmp(button, 'open')
        buttonID = 4;
    else
        buttonID = -1;
    end
end
function buttonDownHandle(hObject, ~)
    init_complete = getCommon('init_complete');
    if init_complete
        clickXY = get(hObject,'CurrentPoint');
        clickB = getButtonID(hObject);
        switch clickB
            case 1 % left mouse button - zoom in
%                panelID = getHPanelByYPosition(clickXY(2), gui);
                setCommon('buttonDownXY', clickXY);
            case 2 % right mouse button - zoom out
                markedStart = 1;
                markedEnd = length(getCommon('currChannelData'));
%                timeInfo = getCommon('currTimeInfo');
%                markedEnd = length(timeInfo);
                updateLimits([markedStart markedEnd], getCommon('xlimitPanelsData'));
%                updateLimits(timeInfo([markedStart markedEnd]), getCommon('xlimitPanelsData'));
            case 3 % middle moudse button - toggle on/off
                gui = guidata(hObject);
                panelID = getHPanelByYPosition(clickXY(2), gui);
                panelSelect(panelID);
        end
    end
end
function [vpi, hpi] = panelID2VpiHpi(panelID, gui)
    vpi = floor(panelID / gui.numVertPanel) + 1;
    hpi = mod(panelID, gui.numVertPanel) + 1;
end 
function [retMarkedStart, retMarkedEnd] = clickCoords2Idx(hStart, hEnd, panelID, gui)
    if nargin < 4
        gui = guidata(getCommon('hfig'));
        if nargin < 3
            panelID = 0;
        end
    end
    if length(panelID) == 2
        vpi = panelID(1);
        hpi = panelID(2);
    else
        [vpi, hpi] = panelID2VpiHpi(panelID, gui);
    end    

%    xlim = get(gui.panel(vpi, hpi), 'XLim');
    xlim = getCommon('currXLimit');
    if isempty(xlim)
        xlim = get(gui.panel(vpi, hpi), 'XLim');
    end
    currIntervalLength = xlim(end) - xlim(1) + 1;
    currPanelLength = gui.panelDimX(vpi, hpi);
%    currXlim = getCommon('currXLimit');
    markedStartOld = xlim(1); % currXlim(1); % markedEndOld = markedEnd;
    retMarkedStart = max(1,markedStartOld-1 + round(min(1,(hStart-gui.panelOffsetX(vpi, hpi)) / currPanelLength) * currIntervalLength));
    retMarkedEnd =   max(2,markedStartOld-1 + round(min(1,(hEnd-  gui.panelOffsetX(vpi, hpi)) / currPanelLength) * currIntervalLength));
    if retMarkedStart == retMarkedEnd
        retMarkedEnd = retMarkedStart + 1;
    end

%                     xlim = get(dataPanel(1), 'XLim');
%                     currIntervalLength = xlim(end) - xlim(1) + 1;
%                     currPanelLength = figureWidth * gui.panelWidthFactor;
%                     markedStartOld = markedStart; % markedEndOld = markedEnd;
%                     retMarkedStart = max(1,markedStartOld-1 + round(min(1,(hStart-figureWidth*gui.panelOffset) / currPanelLength) * currIntervalLength));
%                     retMarkedEnd = max(2,markedStartOld-1 + round(min(1,hEnd / currPanelLength) * currIntervalLength));
%                     if retMarkedStart == retMarkedEnd
%                             retMarkedEnd = retMarkedStart + 1;
%                     end                            
end
function leftClickHandle(pos)
    setCommon('currPos', pos);
%    set(getUIHandle('zoomButtonEditPos'), 'String', num2str(pos));
    timeInfo = getCommon('currTimeInfo');
    set(getUIHandle('zoomButtonEditPos'), 'String', num2str(timeInfo(pos)));
    setGreenLine(pos);
    updateSelectedEventSummary(getCommon('allEvents'), pos);
end        
function buttonUpHandle(hObject, ~)
    init_complete = getCommon('init_complete');    
    if init_complete
%        msgbox({'Mouselocation:',num2str(get(gcbf,'CurrentPoint')),get(gcbf, 'SelectionType')});
%        fprintf('mouse: %s %s\n',num2str(get(gcbf,'CurrentPoint')),get(gcbf, 'SelectionType'));
        clickXY = get(gcbf,'CurrentPoint');
        clickB = getButtonID(gcbf);
        if clickB == 1
            gui = guidata(hObject);
            [panelID, columnID] = getHPanelByYPosition(clickXY, gui);
            if panelID > 0 && ~isempty(find(getCommon('xlimitPanels') == gui.panel(panelID, columnID), 1))
                buttonDownXY = getCommon('buttonDownXY');
                if ~isempty(buttonDownXY)
                    hStart = buttonDownXY(1);
                    hEnd = clickXY(1);
                    if hStart > hEnd
                        h = hStart;
                        hStart = hEnd;
                        hEnd = h;
                    end
                    if hStart ~= hEnd
                        % set drag & dropped zoom interval
                        [markedStart, markedEnd] = clickCoords2Idx(hStart, hEnd, [panelID columnID], gui);
                        updateLimits([markedStart markedEnd], getCommon('xlimitPanelsData'));
                    else
                        % jump to clicked frame
%                         setCurrFrame(clickCoords2Idx(hStart, NaN));
%                         guidata(hfig, g);
%                         videoPic.frames(1) = loadFrame(g.frames.curr, videoPic.frames(1));
%                         updateImages(hfig);
                        pos = clickCoords2Idx(hStart, NaN);
                        leftClickHandle(pos);
                    end
                    setCommon('buttonDownXY', []);
%                    guidata(gcbf, g);
                end
            end
        end
    end
end

function setMousePointer(handle, newPointer, freeze)
    if strcmp(newPointer, 'back')
        oldpointer = getCommon('oldpointer');
        set(handle, 'pointer', oldpointer);
        drawnow;
    else
        oldpointer = get(handle, 'pointer');
        set(handle, 'pointer', newPointer)
        drawnow;
        setCommon('oldpointer', oldpointer);
    end
    if nargin > 2 && ~isempty(freeze)
        warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
        if freeze % disable / enable whole figure
            jFigPeer = get(handle,'JavaFrame');
            jFigPeer.getFigurePanelContainer.setEnabled(false)
        else
            jFigPeer = get(handle,'JavaFrame');
            jFigPeer.getFigurePanelContainer.setEnabled(true);
        end
    end
end
function outSelectAllCallback(hObject, ~)
    [channelListBoxData, clbh] = getChannelListBoxData(hObject);
    set(clbh, 'value', find(channelListBoxData.ischannel));
    channelListBoxCallback(clbh);
end
function outButtonSaveCallback(hObject, ~)
    channelListBoxData = getChannelListBoxData(hObject);
    xlimits = getCommon('currXLimit');    
    setMode('batchProcess');% setBatchProcessMode();
    
    setMousePointer(get(hObject, 'parent'), 'watch');
    runBatch(channelListBoxData.currChannels(channelListBoxData.currChannel), xlimits);
    setMousePointer(get(hObject, 'parent'), 'back');
end
function outButtonSaveAllCallback(hObject, ~)
	channelListBoxData = getChannelListBoxData(hObject);
    
    setMode('batchProcess');% setBatchProcessMode();
    setCommon('bulkEvents', []);
    
	setMousePointer(get(hObject, 'parent'), 'watch', true);
    done = runBatch(channelListBoxData.currChannels);
    setMousePointer(get(hObject, 'parent'), 'back', false);
    
    if get(getUIHandle('optExitDone'), 'Value');
        close(get(hObject, 'parent')); % close GUI
    end
    
    [warnLengthIdx, warnFreqIdx] = updateBulkEventSummary();
    flagWarnings(warnLengthIdx, done, 'Data sample length of \n   %s \nnot within 20 %% of other sample lengths');
    flagWarnings(warnFreqIdx, done, 'Data sampling frequency of \n   %s \nnot within 20 %% of other sample frequencies');
    
    writeRasterPlots();
    writeExcelSummary();
    
    function flagWarnings(idx, done, text)
        for wi = 1 : length(idx)
            msgbox(sprintf(text, done{wi}), 'Warning');
        end
    end
end

function dirButton(editHandle, editCallback)
%        editHandle = getUIhandle(hObject, 'outEditDir');
    start_path = get(editHandle, 'String');
    folder_name = uigetdir(start_path);
    set(editHandle, 'String', folder_name);
    editCallback(editHandle);
%    outEditDirCallback(editHandle);
end
function dirEditButtonCallback(hObject, ~)
    dirButton(getUIHandle(hObject, 'dirEditBox'), @dirEditBoxCallback);    
end
function outButtonDirCallback(hObject, ~)
    dirButton(getUIHandle(hObject, 'outEditDir'), @outEditDirCallback);
%     editHandle = getUIhandle(hObject, 'outEditDir');
%     start_path = get(editHandle, 'String');
%     folder_name = uigetdir(start_path);
%     set(editHandle, 'String', folder_name);
%     outEditDirCallback(editHandle);
end
function outEditDirCallback(hObject, ~)
    outputDir = get(hObject, 'String'); 
    if exist(outputDir) ~= 7 %#ok<EXIST>        
        choice = questdlg('Directory does not exist. Create it?', ...
                            'Output Directory', ...
                            'Yes','No','No');
        switch choice
            case 'Yes'
                [status,message,messageid] = mkdir(outputDir);
                if ~status
                    error(messageid, message);
                    outputDir = [];
                end                
            case 'No'
                fprintf('Use chose not to create directory\n');
%                outputDir = [];
            otherwise
                fprintf('Warning: This should never happen: Create Output Directory returns value other than ''Yes'' or ''No''; going with ''No''\n');
%                outputDir = [];
        end
    end
    if isempty(outputDir)
        outputDir = pwd;
    end
    set(hObject, 'String', outputDir);
	setCommon('outputDir', outputDir);
end

function bulkUISet(handleNames, fieldName, value)
    parent = getCommon('hfig');
    gui = get(parent, 'userdata');
    for hi = 1 : length(handleNames)
        set(gui.handles.(handleNames{hi}), fieldName, value);
    end
end
function ret = getUIHandle(hObject, handleName, p_getParent)
    if nargin == 1
        parent = getCommon('hfig');
        handleName = hObject;
    else
        if nargin < 3
            p_getParent = true;
        end
        if p_getParent
            parent = get(hObject, 'Parent');
        else
            parent = hObject;
        end
    end
 	gui = get(parent, 'userdata');
	ret = gui.handles.(handleName);
end
function [ret, clbh] = getChannelListBoxData(hObject)
    clbh = getUIHandle(hObject, 'channelListBox');
    ret = get(clbh, 'userdata');
end
function setChannelListBoxData(hObject, channelListBoxData)
    set(getUIHandle(hObject, 'channelListBox'), 'userdata', channelListBoxData)
end

function setMode(item, value)
    mode = getCommon('mode');
    if nargin == 2
        if isfield(mode, item)
            mode.(item) = value;
%            mode.drawFIG = mode.saveFIG || mode.generateFIG;
        else
            mode = value;
        end
        if ~mode.saveFIG && ~mode.saveXLS && ~mode.generateFIG && ~mode.saveRaster && ~mode.saveRasterSum && ~mode.saveExcelSum
            bulkUISet({'outButtonSave', 'outButtonSaveAll'}, 'enable', 'off');
        else
            bulkUISet({'outButtonSave', 'outButtonSaveAll'}, 'enable', 'on');
        end                    
    else
        switch(item)
            case 'guiProcess'
                	mode.drawGUI = true;
                    mode.saveXLS = false;
                    mode.saveFIG = false;
                    mode.generateFIG = false;
                    mode.saveRaster = false;
                    mode.saveRasterSum = false;
                    mode.saveExcelSum = false;
%                    mode.drawFIG = false;
                    mode.noEmpty = true;
                    mode.noOverwrite = false;
                    mode.exitDone = false;                    
                    mode.collectEvents = false;
            case 'batchProcess'
                    mode.drawGUI = false;
                    mode.saveXLS = get(getUIHandle('optSaveExcel'), 'Value');
                    mode.saveFIG = get(getUIHandle('optSaveFig'), 'Value');
                    mode.generateFIG = get(getUIHandle('optGenFig'), 'Value');
                    mode.saveRaster = get(getUIHandle('optSaveRaster'), 'Value');
                    mode.saveRasterSum = get(getUIHandle('optSaveRasterSum'), 'Value');
                    mode.saveExcelSum = get(getUIHandle('optSaveExcelSum'), 'Value');
%                    mode.drawFIG = mode.saveFIG || mode.generateFIG;
                    mode.noEmpty = get(getUIHandle('optNoEmpty'), 'Value');
                    mode.noOverwrite = get(getUIHandle('optNoOverwrite'), 'Value');
                    mode.exitDone = get(getUIHandle('optExitDone'), 'Value');
                    mode.collectEvents = true;
            case 'batchOnly'
                    mode.drawGUI = false;
                    mode.saveXLS = true;
                    mode.saveFIG = true;                    
                    mode.generateFIG = false; % ~p_close;
                    mode.saveRaster = true;
                    mode.saveRasterSum = true;
                    mode.saveExcelSum = true;
%                    mode.drawFIG = true;
                    mode.noEmpty = true;
                    mode.noOverwrite = false;
                    mode.exitDone = true;
                    mode.collectEvents = true;
        end
    end
	mode.drawFIG = mode.saveFIG || mode.generateFIG;        
    setCommon('mode', mode);
end
function optNoEmptyCallback(hObject, ~)
    setMode('noEmpty', get(hObject, 'value'));
end
function optNoOverwriteCallback(hObject, ~)
    setMode('noOverwrite', get(hObject, 'value'));
end
function optSaveRasterCallback(hObject, ~)
    setMode('saveRaster', get(hObject, 'value'));
end
function optSaveRasterSumCallback(hObject, ~)
    setMode('saveRasterSum', get(hObject, 'value'));
end
function optSaveExcelSumCallback(hObject, ~)
    setMode('saveExcelSum', get(hObject, 'value'));
end
function optExitDoneCallback(hObject, ~)
    value = get(hObject, 'value');
    setMode('exitDone', value);
    if value
        ogfh = getUIHandle(hObject, 'optGenFig');
        set(ogfh, 'value', false, 'enable', 'off');
        optGenFigCallback(ogfh);
    else
        set(getUIHandle(hObject, 'optGenFig'), 'enable', 'on');
    end
end
function optDrawLegendCallback(hObject, ~)
    value = get(hObject, 'value');
    setMode('drawLegend', value);
    
	legendHandles = getCommon('xlimitPanelsData');
    for pi = 1 : length(legendHandles)
        legend(legendHandles(pi).handle, 'toggle');
    end
    
%    if value
%        valueStr = 'on';
%    else
%        valueStr = 'off';
%    end
%    for pi = 1 : length(legendHandles)
%        legend(legendHandles(pi).handle, 'visible', valueStr);
%        set(get(get(legendHandles(pi), 'annotation'), 'legendInformation'), 'IconDisplayStyle', valueStr)
%    end
end

function optSaveFigCallback(hObject, ~)
%     mode = getCommon('mode');
%     mode.saveFIG = get(hObject, 'value');
%     setCommon('mode', mode);
    setMode('saveFIG', get(hObject, 'value'));
end
function optSaveExcelCallback(hObject, ~)
    setMode('saveXLS', get(hObject, 'value'));
end
function optGenFigCallback(hObject, ~)
    setMode('generateFIG', get(hObject, 'value'));
end

function updatePanels(channelListBoxHandle)
    channelListBoxCallback(channelListBoxHandle);
end
function ret = assignBoundedValue(value, limits)
    ret = NaN(size(value));
    for vi = 1 : length(value)
        ret(vi) = max(limits(1), min(value(vi), limits(2)));
    end
end
function updateLimits(newLimit, xlimitPanelsData)
    if nargin < 2
        xlimitPanelsData = [];
        if nargin < 1
            newLimit = [NaN NaN];
        end
    end
    
    xlim = getCommon('currXLimit');
    timeInfo = getCommon('currTimeInfo');
    xlim(~isnan(newLimit)) = newLimit(~isnan(newLimit));
%     if ~isnan(newLimit(1))
%         xlim(1) = newLimit(1);
%     end
%     if ~isnan(newLimit(2))
%         xlim(2) = newLimit(2);
%     end
    if isempty(xlim)
         bulkUISet({'zoomButtonEditFrom', 'zoomButtonEditTo'}, 'String', '');
    else
        setCommon('currXLimit', xlim);

        % set zoomButtonEditFrom and zoomButtonEditTo
%         set(getUIHandle('zoomButtonEditFrom'), 'String', num2str(xlim(1)))
%         set(getUIHandle('zoomButtonEditTo'), 'String', num2str(xlim(2)))
         set(getUIHandle('zoomButtonEditFrom'), 'String', num2str(timeInfo(xlim(1))))
         set(getUIHandle('zoomButtonEditTo'), 'String', num2str(timeInfo(xlim(2))))
        for hi = 1 : length(xlimitPanelsData)
            set(xlimitPanelsData(hi).handle, 'XLim', timeInfo(xlim));
            adjustYLimit(xlim, [], xlimitPanelsData(hi).handle, hi);
            labelInsideAxes(xlimitPanelsData(hi).handle, xlimitPanelsData(hi).axes, xlimitPanelsData(hi).axesLabel, xlimitPanelsData(hi).offset);
        end        
        currPos = getCommon('currPos');
        if ~isempty(currPos) && currPos >= xlim(1) && currPos <= xlim(2)
            setGreenLine(currPos);
        end
        updateZoomButtons(xlim);
        
        updateCurrEventSummary(getCommon('currEvents'), xlim);
    end
end
function updateCurrEventSummary(events, interval)
    if isempty(events)
        bulkUISet({'eventCurrLatencyV', 'eventCurrDurationV', 'eventCurrVolumeV', 'eNavEdit'}, 'String', '');
        bulkUISet({'eventCurrLatencyV', 'eventCurrDurationV', 'eventCurrVolumeV', 'eventCurrLatencyT', 'eventCurrDurationT', 'eventCurrVolumeT', 'eventCurrText'}, 'visible', 'off');
    else
        
        eventList = myGetFilteredEventList(events, interval);
        
        bulkUISet({'eventCurrLatencyV', 'eventCurrDurationV', 'eventCurrVolumeV', 'eventCurrLatencyT', 'eventCurrDurationT', 'eventCurrVolumeT', 'eventCurrText'}, 'visible', 'on');
        if isempty(eventList)
            bulkUISet({'eventCurrLatencyV', 'eventCurrDurationV', 'eventCurrVolumeV'}, 'String', '')
        else
            set(getUIHandle('eventCurrLatencyV'), 'String', formatUIValue(eventList(1).startTime, 's'));
            set(getUIHandle('eventCurrDurationV'), 'String', formatUIValue(sum(cat(1, eventList.durationTime)), 's'));
            set(getUIHandle('eventCurrVolumeV'), 'String', formatUIValue(abs(sum(cat(1, eventList.volumeS))), 'nl'));
        end
        set(getUIHandle('eventCurrText'), 'String', getCurrEventHeader(length(eventList)));
    end
end
function updateSelectedEventSummary(events, pos)
    function emptySelectEventSummary()
        bulkUISet({'eventSelectedLatencyV', 'eventSelectedDurationV', 'eventSelectedVolumeV'}, 'String', '');
        bulkUISet({'eventSelectedLatencyV', 'eventSelectedDurationV', 'eventSelectedVolumeV', ...
                   'eventSelectedLatencyT', 'eventSelectedDurationT', 'eventSelectedVolumeT', ...
                   'eventSelectedText'}, 'visible', 'off');        
    end
    
    if isempty(events)
        emptySelectEventSummary();
    else
        if nargin < 2
            pos = getCommon('currPos');
        end
        eIdx = find(pos >= cat(1, events.eventList.startIdx), 1, 'last');
        bulkUISet({'eNavButtonLeft', 'eNavButtonRight'}, 'enable', 'on');        
        if isempty(eIdx) || pos > events.eventList(eIdx).endIdx
            emptySelectEventSummary();
            if isempty(eIdx)
                eIdx = 0;
            end
            set(getUIHandle('eNavEdit'), 'String', sprintf('%d+', eIdx));
            if eIdx == 0
                bulkUISet({'eNavButtonLeft'}, 'enable', 'off');
            end
        else
            set(getUIHandle('eventSelectedLatencyV'), 'String', formatUIValue(events.eventList(eIdx).startTime, 's'));
            set(getUIHandle('eventSelectedDurationV'), 'String', formatUIValue(events.eventList(eIdx).durationTime, 's'));
            set(getUIHandle('eventSelectedVolumeV'), 'String', formatUIValue(abs(events.eventList(eIdx).volumeS), 'nl'));
            bulkUISet({'eventSelectedLatencyV', 'eventSelectedDurationV', 'eventSelectedVolumeV', ...
                       'eventSelectedLatencyT', 'eventSelectedDurationT', 'eventSelectedVolumeT', ...
                       'eventSelectedText'}, 'visible', 'on');
            set(getUIHandle('eNavEdit'), 'String', sprintf('%d', eIdx));
            if eIdx == 1
                bulkUISet({'eNavButtonLeft'}, 'enable', 'off');
            end
        end        
        if eIdx == length(events.eventList)
        	bulkUISet({'eNavButtonRight'}, 'enable', 'off');                
        end
    end
end
function [totalVolumeS, totalDurationS] = calcGroupStats(eventList)
        totalDurationS = sum(cat(1, eventList.durationTime));
        totalVolumeS = abs(sum(cat(1, eventList.volumeS)));    
end
function [warnLengthIdx, warnFreqIdx] = updateBulkEventSummary(bulkEvents)
    function [groupStat, warnLengthIdx, warnFreqIdx]= getEventGroupStats(bulkEvents)
        groupStat.totalChannels = length(bulkEvents.stats);
        
        eventList = bulkEvents.eventList;
        groupStat.totalEvents = length(eventList);        
        allLatencies = cat(1, bulkEvents.stats.latency);
        groupStat.latency(1) = nanmean(allLatencies);
        groupStat.latency(2) = sum(~isnan(allLatencies));
%        groupStat.durationS = sum(cat(1, eventList.durationTime));
%        groupStat.volumeS = abs(sum(cat(1, eventList.volumeS)));
        [groupStat.volumeS, groupStat.durationS] = calcGroupStats(eventList);

        allLengths = cat(1, bulkEvents.stats.length);
        groupStat.feedingPerc = groupStat.durationS / sum(allLengths) * 100;
        allInterMeals = cat(1, bulkEvents.stats.avgInterMeal);
        groupStat.bulkAvgInterMeal(1) = nanmean(allInterMeals);
        groupStat.bulkAvgInterMeal(2) = sum(~isnan(allInterMeals));
        groupStat.bulkAvgSpeed = nanmean(cat(1, bulkEvents.stats.avgSpeed));
        
        % check range of data sample lengths, flag outliers of more than 20 %
        warnLengthIdx = find(LOOcheck(allLengths, getParam('maxLengthVariation')));
        warnFreqIdx = find(LOOcheck(cat(1, bulkEvents.stats.freq), getParam('maxFreqVariation')));
        	
        function invalid = LOOcheck(allValues, maxThresh)
            invalid = false(size(allValues));
            for i = 1 : length(allValues)
                currValue = allValues(i);
                currValueSet = allValues; currValueSet(i) = [];
                avgValue = nanmean(currValueSet);
                
                invalid(i) = abs(currValue-avgValue) / avgValue >= maxThresh;
            end
        end
    end
    
    if nargin < 2
        bulkEvents = getCommon('bulkEvents');
    end
    eventBulkUINames = {'eventBulkLatencyV', 'eventBulkDurationV', 'eventBulkVolumeV', ...
                   'eventBulkLatencyT', 'eventBulkDurationT', 'eventBulkVolumeT', ...
                   'eventBulkText', 'eventBulkText2' ...
                   'eventBulkInterEvtT', 'eventBulkEvtPercT', 'eventBulkSpeedT', ...
                   'eventBulkInterEvtV', 'eventBulkEvtPercV', 'eventBulkSpeedV', ...                   
                   };
    if isempty(bulkEvents)
        bulkUISet({'eventBulkLatencyV', 'eventBulkDurationV', 'eventBulkVolumeV'}, 'String', '');        
        bulkUISet(eventBulkUINames, 'visible', 'off');
        
        warnLengthIdx = []; warnFreqIdx = [];
    else
        bulkUISet(eventBulkUINames, 'visible', 'on');
        
        [groupStat, warnLengthIdx, warnFreqIdx] = getEventGroupStats(bulkEvents);
        
        set(getUIHandle('eventBulkText'), 'String', getBulkEventHeader(groupStat.totalEvents, groupStat.totalChannels));
        set(getUIHandle('eventBulkLatencyV'), 'String', sprintf('%s, %d ch', formatUIValue(groupStat.latency(1), 's'), groupStat.latency(2)));
        set(getUIHandle('eventBulkDurationV'), 'String', formatUIValue(groupStat.durationS, 's'));
        set(getUIHandle('eventBulkVolumeV'), 'String', formatUIValue(groupStat.volumeS, 'nl'));
        
        set(getUIHandle('eventBulkInterEvtV'), 'String', sprintf('%s, %d ch', formatUIValue(groupStat.bulkAvgInterMeal(1), 's'), groupStat.bulkAvgInterMeal(2)));
        set(getUIHandle('eventBulkEvtPercV'), 'String', formatUIValue(groupStat.feedingPerc, '%'));
        set(getUIHandle('eventBulkSpeedV'), 'String', formatUIValue(groupStat.bulkAvgSpeed, 'nl/s'));
    end
end
function writeRasterPlots(bulkEvents)
	if nargin < 2
        bulkEvents = getCommon('bulkEvents');
    end
    mode = getCommon('mode');
    
    rasterFileName = sprintf('%s/%f.eps', outputDir, now);
    
%                 bar(markedStart : markedEnd, dataSelector(panelID).values2(markedStart : markedEnd)*-1,'barWidth',1,'LineWidth',0.01,'LineStyle','none','facecolor','r','edgecolor','r','ShowBaseLine','off');
%                 hold on;
%                 bar(markedStart : markedEnd, dataSelector(panelID).values1(markedStart : markedEnd),'barWidth',1,'LineWidth',0.01,'LineStyle','none','facecolor','b','edgecolor','b','ShowBaseLine','off');
%                 hold off; 
%                 set(dataPanel(panelID),'ylim',[-1 1]);

%% init pic with bgColor
%    bgColor = [240 240 240];
    rowWidth = 10;
    fillWidth = 6;
    
    numRows = length(bulkEvents.stats);    
    minLength = min(cat(1, bulkEvents.stats.length));
%    maxHeight = (rowWidth+fillWidth) * numRows;
    
    M = NaN(numRows, minLength);
    
%     pic = zeros(maxHeight, minLength, 3);
%     for ci = 1 : length(bgColor)
%         pic(:,:,ci) = bgColor(ci);
%     end        
%% draw event bars
    if mode.saveRaster
%        rasterPlotHandle = figure('Name', 'Raster Plot', 'visible', 'off');
        rasterPlotHandle = figure('Name', [getCommon('titlePrefix') 'Raster Plot']);
    end

    eventStart = 1;
    tickLabels = {};
    minMaxTime = Inf;
    for rowI = 1 : numRows
        eventEnd = bulkEvents.stats(rowI).total;
        currEvents = bulkEvents.eventList(eventStart : eventStart + eventEnd - 1);
        
        row = NaN(1, bulkEvents.stats(rowI).length);
%        rowV = (rowWidth + fillWidth) * (rowI - 1);
        rowV = numRows - rowI + 1;
        for eventI = 1 : length(currEvents)
%            pic = drawEventBar(pic, currEvents(eventI).startIdx, currEvents(eventI).endIdx, rowI, rowWidth, fillWidth, [0 0 255]);
            row(currEvents(eventI).startIdx : currEvents(eventI).endIdx) = rowV;
        end
        M(rowI, :) = ~isnan(row(1:minLength));
%        row = row(1 : minLength);

        eventStart = eventStart + eventEnd;
        
        if mode.saveRaster
            if rowI > 1
                hold on;
            end

            plot(bulkEvents.info(rowI).timeInfo, row, 'b', 'linewidth', rowWidth);
            %bar(bulkEvents.stats(rowI).timeInfo, row, 'b', 'barwidth', 1);

            currLabel = bulkEvents.info(rowI).fileText; %sprintf('test%d', rowI);
            tickLabels = [{currLabel} tickLabels];

            minMaxTime = min(minMaxTime, bulkEvents.info(rowI).timeInfo(end));
        end        
    end
    if mode.saveRaster
        set(gca, 'YLim', [-0.5 numRows+1])
        set(gca, 'YTick', 1 : numRows);
        set(gca, 'YTickLabel', [tickLabels]);
        set(gca, 'XLim', [0 minMaxTime]);
        set(gca, 'Position', [0.33 0 0.67 1]);
        set(gca, 'TickLength', [0.005 0.0250]);
        labelInsideAxes(gca, 'x', {'[sec]'});    

        hfig = getCommon('hfig');    
        gui = get(hfig, 'userdata');
        rasterHeight = (numRows+1.5)*(fillWidth+rowWidth);
        set(rasterPlotHandle, 'Position', [0 gui.screenHeight-rasterHeight-100 gui.screenWidth rasterHeight]);

        %print(rasterPlotHandle, '-dtiff', rasterFileName);
        print(rasterPlotHandle, '-depsc', rasterFileName);
        fprintf('raster plots saved in ''%s''\n', rasterFileName);
    end
    
%    figure, imshow(pic);
    if mode.saveRasterSum
        rasterSummaryHandle = figure('Name', [getCommon('titlePrefix') 'Raster Summary'], 'visible', 'off');
%        % old raster summary plot
%         bar(sum(M, 1));
%         ylim([0 numRows]);
%         xlim([0 minLength]);
%         set(gca, 'ytick', 1:numRows)

        % new raster summary plot
        minBins = Inf;
        for i = 1 : length(bulkEvents.stats)
            minBins = min(minBins, length(bulkEvents.stats(i).binVol));
        end
        cumBins = zeros(minBins, 1);
        for i = 1 : length(bulkEvents.stats)
            cumBins = nansum(cat(2, cumBins, bulkEvents.stats(i).binVol(1 : minBins)), 2);
        end
        binSec = getCommon('rasterSumBinSec');
        plot((1 : minBins) * binSec + binSec/2, abs(cumBins)./length(bulkEvents.stats), '-');
                     rasterSummaryName = [rasterFileName '.summary.eps'];
                     print(rasterSummaryHandle, '-depsc', rasterSummaryName);
                     fprintf('raster summary saved in ''%s''\n', rasterSummaryName);
         close(rasterSummaryHandle);

        
    end

    
    function row = addEventBar(row, startIdx, endIdx, rowI, rowWidth, fillWidth)
        rowOffset = (rowWidth + fillWidth) * (rowI - 1);
        row(startIdx : endIdx)
    end

    function pic = drawEventBar(pic, startIdx, endIdx, rowI, rowWidth, fillWidth, color)
        rowOffset = (rowWidth + fillWidth) * (rowI - 1) + 1;
        for c = 1 : length(color)
            pic(rowOffset : rowOffset + rowWidth - 1, startIdx : endIdx, c) = color(c);
        end
    end
end
function writeExcelSummary(bulkEvents)
    
    mode = getCommon('mode');

    if mode.saveExcelSum
        if nargin < 2
            bulkEvents = getCommon('bulkEvents');
        end
        excelSumFileName = sprintf('%s/%f.xls', outputDir, now);
        header = {'filename','XP','channel','number of events [#]', 'total volume [nl]','total duration [s]', 'latency [s]', 'avg interMeal [s]', 'halfDrinkAt [s]', 'volume_1_3 [nl]', 'volume_2_3 [nl]', 'volume_3_3 [nl]'};

        C = cell(length(bulkEvents.stats), length(header));
        
        D = {'filename', 'XP', 'channel', 'startIdx [#]', 'endIdx [#]', 'durationIdx [#]', 'startTime [s]', 'stopTime [s]', 'durationTime [s]', 'volume [nl]', 'speed [nl/s]', 'interMeal [s]', 'recordedTime [HH-MM-SS  mm/dd/yyyy]'};
        DnewT = cell(1, length(D));
        
        currEvtStartIdx = 1;
        for fi = 1 : length(bulkEvents.stats)
            total = bulkEvents.stats(fi).total;
            currEvtStopIdx = currEvtStartIdx + total - 1;
            currEvtList = bulkEvents.eventList(currEvtStartIdx : currEvtStopIdx);
            currEvtStartIdx = currEvtStopIdx + 1;
            
            [C{fi, 1}, C{fi, 2}, C{fi, 3}] = parseFileText(bulkEvents.info(fi).fileText);                        
            
            [totalVolumeS, totalDurationS] = calcGroupStats(currEvtList);
            
            C{fi, 4} = total;
            C{fi, 5} = totalVolumeS;
            C{fi, 6} = totalDurationS;
            C{fi, 7} = bulkEvents.stats(fi).latency;
            C{fi, 8} = bulkEvents.stats(fi).avgInterMeal;            
            C{fi, 9} = bulkEvents.stats(fi).halfDrinkAt;            

            C{fi, 10} = bulkEvents.stats(fi).totalVol_1_3;
            C{fi, 11} = bulkEvents.stats(fi).totalVol_2_3;
            C{fi, 12} = bulkEvents.stats(fi).totalVol_3_3;
            
            Dnew = DnewT;
            Dnew(1:3) = C(fi, 1:3);
            D = cat(1, D, Dnew, events2cell(currEvtList, bulkEvents.info(fi).recordingTime));
        end
        xlswrite(excelSumFileName, cat(1, header, C), 'summary');
        xlswrite(excelSumFileName, D, 'events');                
    end       
    
    function [filename, XP, channel] = parseFileText(t)
        Cs = strsplit(t, '''');
        filename = Cs{2};
        Css = strsplit(Cs{3}, '.');
        XP = Css{2};
        channel = Css{3};
    end
    function D = events2cell(events, recordingTime)
        D = cell(length(events), length(DnewT));
        for e = 1 : length(events)            
            D{e, 4}  = events(e).startIdx;
            D{e, 5}  = events(e).endIdx;
            D{e, 6}  = events(e).durationIdx;
            D{e, 7}  = events(e).startTime;
            D{e, 8}  = events(e).stopTime;
            D{e, 9}  = events(e).durationTime;
            D{e, 10} = abs(events(e).volumeS);
            D{e, 11} = events(e).speedS;
            D{e, 12} = myNum2str(events(e).interMeal);
            D{e, 13} = datestr(recordingTime + datenum(0,0,0,0,0,double(events(e).startTime)), 'HH-MM-SS  mm/dd/yyyy');
        end
    end
end

    function ret = myRound(value, prec)
        ret = round(value * 10^prec) / 10.0^prec;
    end
    function ret = formatUIValue(value, unit, prec)
        if nargin < 3
            prec = 3;
        end
        if nargin < 2 || isempty(unit)
            formatSpec = sprintf('%%.%df', prec);
            ret = sprintf(formatSpec, value);
        else
            formatSpec = sprintf('%%.%df [%%s]', prec);            
            ret = sprintf(formatSpec, value, unit);
        end
    end

function zoomButtonEdit(hObject, newLimitI)    
    value = get(hObject, 'String');
%    [iValue, vStatus] = str2num(value);    
	[tValue, vStatus] = str2num(value);        
    iValue = findClosestIdx(tValue);
    
    if ~vStatus
        iValue = 1;
        set(hObject, 'String', iValue);
    end
    newLimits = [NaN NaN];
    newLimits(newLimitI) = iValue;
    updateLimits(newLimits, getCommon('xlimitPanelsData'));
end
function zoomButtonEditFromCallback(hObject, ~)    
    zoomButtonEdit(hObject, 1);
end
function zoomButtonEditToCallback(hObject, ~)    
    zoomButtonEdit(hObject, 2);
%     value = get(hObject, 'String');
%     [iValue, vStatus] = str2num(value);
%     
%     if ~vStatus
%         iValue = 1;
%         set(hObject, 'String', iValue);
%     end
%     updateLimits([iValue NaN], getCommon('xlimitPanelsData'));
end
                function zoomButtonMinus2Callback(hObject,~)
                    gui = guidata(get(hObject, 'Parent'));
                    doZoom(gui.zoomFactor2);
                end
                function zoomButtonMinus3Callback(hObject,~)
                    gui = guidata(get(hObject, 'Parent'));
                    doZoom(gui.zoomFactor3);
                end
                function zoomButtonPlus2Callback(hObject,~)
                    gui = guidata(get(hObject, 'Parent'));
                    doZoom(1 / gui.zoomFactor2);
                end
                function zoomButtonPlus3Callback(hObject,~)
                    gui = guidata(get(hObject, 'Parent'));
                    doZoom(1 / gui.zoomFactor3);
                end
    function checkAndUpdateLimits(xlimits)
                    maxLen = length(getCommon('currChannelData'));
                    if xlimits(2) - xlimits(1) + 1 > maxLen
                        xlimits = [1 maxLen];
                    end
                    if xlimits(1) < 1
                        xlimits = xlimits + (1 - xlimits(1));
                    end
                    if xlimits(2) > maxLen
                        xlimits = xlimits + (maxLen - xlimits(2));
                    end        
                    updateLimits(xlimits, getCommon('xlimitPanelsData'));
    end
                function doZoom(pFactor)
                    xlimits = getCommon('currXLimit');

                    len = xlimits(2) - xlimits(1);
                    destLen = len * pFactor;
                    xlimits(1) = round(xlimits(1) - (destLen - len)/2);
                    xlimits(2) = round(xlimits(2) + (destLen - len)/2);
                    
                    checkAndUpdateLimits(xlimits);
                end

                function zoomButtonLeftCallback(hObject,~)
                    gui = guidata(get(hObject, 'Parent'));
                    doShift(-gui.zoomShiftFactor1);
                end
                function zoomButtonRightCallback(hObject,~)
                    gui = guidata(get(hObject, 'Parent'));
                    doShift(gui.zoomShiftFactor1);
                end
                function doShift(pFactor)
                    xlimits = getCommon('currXLimit');
                    offset = (xlimits(2) - xlimits(1)) * pFactor;
                    
                    xlimits = xlimits + round(offset);                    
                    checkAndUpdateLimits(xlimits);
                end



function navButtonEditCallback(hObject, ~)    
    value = get(hObject, 'String');
    [iValue, vStatus] = str2num(value);
    
    if vStatus
        setCurrChannel(hObject, false, iValue)
    else
        setCurrChannel(hObject, false, 1);
    end
% 	channelListBoxData = getChannelListBoxData(hObject);
%     if vStatus
%         channelListBoxData.currChannel = assignBoundedValue(iValue, [1 length(channelListBoxData.currChannels)]);
%     else
%         channelListBoxData.currChannel = 1;
%     end
% 	set(hObject, 'String', num2str(channelListBoxData.currChannel));
%     setChannelListBoxData(hObject, channelListBoxData);
%     
%     drawPreview(channelListBoxData.currChannels(channelListBoxData.currChannel));
end
function updateNavButtons(curr, limit)
	bulkUISet({'navButtonLeft1', 'navButtonLeft10', 'navButtonRight1', 'navButtonRight10'}, 'enable', 'on');
    if curr <= limit(1)
        bulkUISet({'navButtonLeft1', 'navButtonLeft10'}, 'enable', 'off');
    end
    if curr >= limit(2)
        bulkUISet({'navButtonRight1', 'navButtonRight10'}, 'enable', 'off');
    end
end
function updateZoomButtons(xlim)
    limit = getCommon('maxXLimit');
	bulkUISet({'zoomButtonMinus3', 'zoomButtonMinus2', 'zoomButtonLeft', 'zoomButtonRight'}, 'enable', 'on');
    if xlim(1) <= limit(1)
        bulkUISet({'zoomButtonLeft'}, 'enable', 'off');
    end
    if xlim(2) >= limit(2)
        bulkUISet({'zoomButtonRight'}, 'enable', 'off');
    end
    if all(xlim == limit)
    	bulkUISet({'zoomButtonMinus3', 'zoomButtonMinus2'}, 'enable', 'off');
    end
end
function setCurrChannel(hObject, pRelative, value)
	channelListBoxData = getChannelListBoxData(hObject);
    limits = [1 length(channelListBoxData.currChannels)];
	channelListBoxData.currChannel = assignBoundedValue(channelListBoxData.currChannel * pRelative + value, limits);
    updateNavButtons(channelListBoxData.currChannel, limits);
	set(getUIHandle(hObject, 'navButtonEdit'), 'String', num2str(channelListBoxData.currChannel));
    setChannelListBoxData(hObject, channelListBoxData);
    
    drawPreview(channelListBoxData.currChannels(channelListBoxData.currChannel));
end
function navButtonLeft1Callback(hObject, ~)        
    setCurrChannel(hObject, true, -1);
end
function navButtonLeft10Callback(hObject, ~)        
    setCurrChannel(hObject, true, -10);
end
function navButtonRight1Callback(hObject, ~)        
    setCurrChannel(hObject, true, +1);
end
function navButtonRight10Callback(hObject, ~)        
    setCurrChannel(hObject, true, +10);
end
function setSelectedEvent(hObject, pRelative, value)
    events = getCommon('allEvents');
    if ~isempty(events)
        currStrValue = get(getUIHandle(hObject, 'eNavEdit'), 'String');
        currValue = sscanf(currStrValue, '%d');
        if isempty(currValue)
            currValue = 0;
        end
        currValue = assignBoundedValue(currValue * pRelative + value, [1 length(events.eventList)]);
        pos = events.eventList(currValue).startIdx;
        leftClickHandle(pos);
%        set(navEditHandle, 'String', num2str(currValue));
    end
end
function eNavButtonLeftCallback(hObject, ~)
    currStrValue = get(getUIHandle(hObject, 'eNavEdit'), 'String');
    if isempty(strfind(currStrValue, '+'))
        setSelectedEvent(hObject, true, -1);
    else
        setSelectedEvent(hObject, true, 0);        
    end
end
function eNavButtonRightCallback(hObject, ~)
    setSelectedEvent(hObject, true, +1);    
end
function eNavEditCallback(hObject, ~)
    value = get(hObject, 'String');
    [iValue, vStatus] = str2num(value);
    
    if vStatus
        setSelectedEvent(hObject, false, iValue);
    else
        currPos = getCommon('currPos');
        if isempty(currPos)
            currPos = 0;
        end
        leftClickHandle(currPos);
    end
end
function zoomButtonEditPosCallback(hObject, ~)
    value = get(hObject, 'String');
%    [iValue, vStatus] = str2num(value);
    [tValue, vStatus] = str2num(value);        
    iValue = findClosestIdx(tValue);
    
    if vStatus
        leftClickHandle(iValue)
    else
        leftClickHandle(getCommon('currPos'))
    end    
end
function iValue = findClosestIdx(tValue)
    timeInfo = getCommon('currTimeInfo');
    dist = abs(timeInfo - tValue);
    iValue = find(dist == min(dist));
end
function debugValThreshEditCallback(hObject, ~)
    devfakt = getCommon('devfakt');
    value = get(hObject, 'String');
    [iValue, vStatus] = str2num(value);
    
    if vStatus
        devfakt.event = iValue;
        setCommon('devfakt', devfakt);
        
%         parent = get(hObject, 'parent');
%         gui = get(parent, 'userdata');
%         clbh = gui.handles.channelListBox;    
%        clbh = getUIHandle(hObject, 'channelListBox');
        updatePanels(getUIHandle(hObject, 'channelListBox'));        
    else
        set(hObject, 'String', num2str(devfakt.event));
    end    
end
function updateChannelListBox(dirListBoxHandle)
%     parent = get(dirListBoxHandle, 'parent');
%     gui = get(parent, 'userdata');
%     clbh = gui.handles.channelListBox;
    clbh = getUIHandle(dirListBoxHandle, 'channelListBox');
    [channelListBoxEntries, channelListBoxData] = prepareChannelListBox(dirListBoxHandle);
    set(clbh, 'String', channelListBoxEntries, 'Max', length(channelListBoxEntries), 'userdata', channelListBoxData, 'value', find(channelListBoxData.ischannel, 1));
    
    updatePanels(clbh);
end

function updateChannelListBoxFilter(filterHandle)
% 	parent = get(filterHandle, 'Parent');
%  	gui = get(parent, 'userdata');
%   clbh = gui.handles.channelListBox;
%    clbh = getUIHandle(filterHandle, 'channelListBox');
	[channelListBoxData, clbh] = getChannelListBoxData(filterHandle);
    [channelListBoxEntries, channelListBoxData] = prepareChannelListBox2(channelListBoxData.allChannelList, channelListBoxData.allFileList);
	set(clbh, 'String', channelListBoxEntries, 'Max', length(channelListBoxEntries), 'userdata', channelListBoxData, 'value', find(channelListBoxData.ischannel, 1));
end

function updateFilter(filterName, filterHandle)
    const = getCommon('const');
	p_filter = getCommon('p_filter');
    
    newValue = get(filterHandle, 'String');
    if isempty(regexp(newValue, const.(filterName), 'once' ))        
        set(filterHandle, 'String', '');
        newValue = [];
    end
    p_filter.(filterName) = newValue;
    setCommon('p_filter', p_filter);
    
    updateChannelListBoxFilter(filterHandle);
end
function channelEditBox1Callback(hObject, ~)
    updateFilter('container', hObject);
end
function channelEditBox2Callback(hObject, ~)
    updateFilter('channel', hObject);
end

    
function dirEditBoxCallback(hObject, ~)
    inputDirValue = get(hObject, 'String'); 
%   parent = get(hObject, 'Parent');
%  	gui = get(parent, 'userdata');
% 	dlbh = gui.handles.dirListBox;    
    dlbh = getUIHandle(hObject, 'dirListBox');
         
    [dirListBoxEntries, dirListBoxData] = prepareDirListBox(inputDirValue);
    set(dlbh, 'String', dirListBoxEntries, 'Max', length(dirListBoxEntries), 'userdata', dirListBoxData, 'value', 1);
	updateChannelListBox(dlbh); % clear channel listbox
end
function dirListBoxCallback(hObject, ~)    
    index_selected = get(hObject,'Value');    
    parent = get(hObject, 'Parent');
    selectionType = get(parent,'selectionType');    
%    list = get(hObject,'String');
    
    dirListBoxData = get(hObject, 'userdata');
    
    isdir = index_selected <= dirListBoxData.numDirEntries;
    isfile = ~isdir;
    if sum(isdir) >= 1 
        if sum(isfile) >= 1 % directories and files selected -> remove directories
            index_selected(isdir) = [];
            set(hObject, 'Value', index_selected);            
        else % only directories selected
            index_selected(1 : end-1) = [];
            set(hObject, 'Value', index_selected);            
            
            item_selected = dirListBoxData.dirNames{index_selected};
            fprintf('selected %d: %s %s\n', index_selected, item_selected, selectionType);
            
            if strcmp(selectionType, 'open')
                % change directory
                newInputDir = getNewDirectory(item_selected);
                fprintf('new inputDir: ''%s''\n', newInputDir);
%                gui = get(parent, 'userdata');
%                debh = gui.handles.dirEditBox;
                debh = getUIHandle(hObject, 'dirEditBox');
                set(debh, 'String', newInputDir);
                dirEditBoxCallback(debh);
                return;
            end            
            updateChannelListBox(hObject); % clear channel listbox
        end
    end
    if sum(isfile) >= 1 % only files selected (directories have been removed)
        files_selected = index_selected - dirListBoxData.numDirEntries;
        fileList = dirListBoxData.fileNames(files_selected);
        for li = 1 : length(files_selected)
            fprintf('selected %d: %s\n', files_selected(li), fileList{li});
        end
%         % update channel listbox
%         gui = get(parent, 'userdata');
%         clbh = gui.handles.channelListBox;
%         [channelListBoxEntries, channelListBoxData] = prepareChannelListBox(hObject);
%         set(clbh, 'String', channelListBoxEntries, 'Max', length(channelListBoxEntries), 'userdata', channelListBoxData, 'value', find(channelListBoxData.ischannel, 1));        
        updateChannelListBox(hObject);
    end
        
    function ret = getNewDirectory(dirName, inputDir)
        if nargin < 2
            inputDir = getCommon('inputDir');
        end

        switch dirName
            case '.'
                ret = inputDir;
            case '..'
                pathSepI = sort([strfind(inputDir, '\') strfind(inputDir, '/')]);
                if ~isempty(pathSepI)
                    ret = inputDir(1 : pathSepI(end)-1);
                    if (ispc && length(ret)==2) || isunix % root
                        ret = [ret filesep];
                    end
                else
                    ret = pwd;
                end
            otherwise
                if inputDir(end) ~= filesep
                    ret = [inputDir filesep dirName];
                else
                    ret = [inputDir dirName];
                end
        end
    end
end
function channelListBoxCallback(hObject, ~)    
    index_selected = get(hObject,'Value');
    channelListBoxData = get(hObject, 'userdata');
    
    newSelection = [];
    for li = 1 : length(index_selected)
        newSelection = [newSelection; channelListBoxData.imap(index_selected(li)).channelI];
    end
    newSelection = unique(newSelection);    
    newSelectionV = newSelection;
    for li = 1 : length(newSelection)
        newSelectionV(li) = channelListBoxData.smap(newSelection(li));
    end
    set(hObject, 'Value', newSelectionV);
    
	channelListBoxData.currChannels = channelListBoxData.channelsFiltered(newSelection);    
    numChannels = length(channelListBoxData.currChannels);
    if numChannels == 0
        channelListBoxData.currChannel = 0;
        set(getUIHandle(hObject, 'navButtonText'), 'String', '');
        set(getUIHandle(hObject, 'navButtonEdit'), 'String', '', 'Enable', 'off');
        set(getUIHandle(hObject, 'outSelectAll'), 'Enable', 'off');
    else
        channelListBoxData.currChannel = 1;
        set(getUIHandle(hObject, 'navButtonText'), 'String', num2str(numChannels));
        set(getUIHandle(hObject, 'navButtonEdit'), 'String', num2str(channelListBoxData.currChannel), 'Enable', 'on');
        set(getUIHandle(hObject, 'outSelectAll'), 'Enable', 'on');        
    end
	set(hObject, 'userdata', channelListBoxData);
    
%     if numChannels > 1    
%         bulkUISet({'navButtonLeft1', 'navButtonLeft10', 'navButtonRight1', 'navButtonRight10'}, ...
%                    'enable', 'on');        
% %         set(getUIHandle(hObject, 'navButtonLeft1'), 'Enable', 'on');
% %         set(getUIHandle(hObject, 'navButtonLeft10'), 'Enable', 'on');
% %         set(getUIHandle(hObject, 'navButtonRight1'), 'Enable', 'on');
% %         set(getUIHandle(hObject, 'navButtonRight10'), 'Enable', 'on');        
%     else
%         bulkUISet({'navButtonLeft1', 'navButtonLeft10', 'navButtonRight1', 'navButtonRight10'}, ...
%                    'enable', 'off');                
% %         set(getUIHandle(hObject, 'navButtonLeft1'), 'Enable', 'off');
% %         set(getUIHandle(hObject, 'navButtonLeft10'), 'Enable', 'off');
% %         set(getUIHandle(hObject, 'navButtonRight1'), 'Enable', 'off');
% %         set(getUIHandle(hObject, 'navButtonRight10'), 'Enable', 'off');
%     end
    updateNavButtons(channelListBoxData.currChannel, [1 numChannels]);
    setCommon('bulkEvents', []);
    updateBulkEventSummary([]);
    
    for li = 1 : length(newSelection)
%        if channelListBoxData.ischannel(index_selected(li))
%        	channelListBoxData.allChannelList(channelListBoxData.imap(5).channelI)
%        end
             item_selected = channelListBoxData.allChannelList(newSelection(li));
             fprintf('selected %d: ''%s'' in ''%s''\n', li, item_selected.fieldPath, item_selected.fullpath);
%         end
    end        
    
    if numel(newSelection) == 1
%        setCommon('currYLimit', []);
        drawPreview(channelListBoxData.channelsFiltered(newSelection(li)))
    elseif numel(newSelection) == 0
        clearPreview();
    end
end
function setGreenLine(f)
    xlim = getCommon('currXLimit');
    if ~isempty(xlim)
        timeInfo = getCommon('currTimeInfo');
        if f >= xlim(1) && f <= xlim(2)
            x = timeInfo(f);
%            x = f;
        else
            x = 0;
        end
        greenLinePanels = getCommon('xlimitPanels');
        greenLineHandles = getCommon('greenLineHandles');
        for gi = 1 : length(greenLineHandles)
            set(greenLineHandles(gi),'XData',[x x],'YData',get(greenLinePanels(gi), 'YLim'));
        end
    end
end
function initGreenLine(value)
	greenLinePanels = getCommon('xlimitPanels');
    greenLineHandles = NaN(length(greenLinePanels), 1);
    for pi = 1 : length(greenLinePanels)
        subplot(greenLinePanels(pi));
        greenLineHandles(pi) = line('XData',[0 0],'YData',[0 0], 'color','g');
        set(get(get(greenLineHandles(pi),'Annotation'),'LegendInformation'), 'IconDisplayStyle','off'); % turn off legend for green line
    end        
    setCommon('greenLineHandles', greenLineHandles);    
    if nargin >= 1 && ~isempty(value)
        setGreenLine(value);
    else
        currPos = getCommon('currPos');
        if ~isempty(currPos)
            setGreenLine(currPos);
        end
    end
end
function drawPreview(channelListItem)
    setMode('guiProcess');
    
    setCommon('currXLimit', []);
    setCommon('currEvents', []);
    setCommon('allEvents', []);
    setCommon('currPos', []);
    
	setMousePointer(getCommon('hfig'), 'watch');
	runBatch(channelListItem);
    setMousePointer(getCommon('hfig'), 'back');    
    
	setMode('batchProcess');% setBatchProcessMode();
    
    bulkUISet({'zoomButtonMinus3', 'zoomButtonMinus2', 'zoomButtonEditFrom', 'zoomButtonEditTo', 'zoomButtonPlus2', 'zoomButtonPlus3', ...
               'zoomButtonEditPos', 'zoomButtonLeft', 'zoomButtonRight', ...
               'eNavButtonLeft', 'eNavButtonRight', 'eNavEdit'}, ...               
               'enable', 'on');
    bulkUISet({'zoomButtonEditPos', 'eNavEdit'}, 'String', '');
           
    initGreenLine();
    updateLimits();        
    leftClickHandle(1);
%	setCommon('currPos', 1);    
%    updateSelectedEventSummary(getCommon('allEvents'));
end
function clearPreview
    hfig = getCommon('hfig');
    set(hfig, 'Name', getCommon('titlePrefix'));

    setCommon('currXLimit', []);
    setCommon('currEvents', []);    
    setCommon('allEvents', []);    
    
	setCommon('greenLineHandles', []);        
    panel = getCommon('panel');
	emptyPanel(panel(1,1));
    emptyPanel(panel(2,1));
    emptyPanel(panel(3,1));
        
    setCommon('currXLimit', []);
    updateLimits();
    updateCurrEventSummary([]);
    updateSelectedEventSummary([]);
    bulkUISet({'zoomButtonMinus3', 'zoomButtonMinus2', 'zoomButtonEditFrom', 'zoomButtonEditTo', 'zoomButtonPlus2', 'zoomButtonPlus3', ...
               'zoomButtonEditPos', 'zoomButtonLeft', 'zoomButtonRight', ...
               'eNavButtonLeft', 'eNavButtonRight', 'eNavEdit'}, ...
               'enable', 'off');
    bulkUISet({'zoomButtonEditPos', 'eNavEdit'}, 'String', '');
end
function emptyPanel(subPanel)
    if nargin >= 1 && ~isempty(subPanel)
        subplot(subPanel);
    end
    cla;
    if nargin < 2
        plot(NaN, NaN); % clear panel
    else
        plot(NaN, xlimits(1):xlimits(2));
    end
    axis off;
end

function ret = htmlItem(text, htmlColor)
    function ret = patchEntry(text, htmlColor)
        ret = sprintf('<HTML><FONT color="%s">%s</FONT></HTML>', htmlColor, text);
    end
    if ischar(text)
        ret = patchEntry(text, htmlColor);
    elseif iscell(text)
        ret = text;
        for ti = 1 : length(text)
            ret{ti} = patchEntry(text{ti}, htmlColor);
        end
    else
        ret = text;
    end        
end
function [dirListBoxEntries, dirListBoxData] = prepareDirListBox(inputDir)
    if nargin < 1
        inputDir = getCommon('inputDir');
    end
    
    [fileList, inputDir, dirList, b_selectFiles] = getFileList(inputDir);
    setCommon('inputDir', inputDir);
    
    dirNames = cat(1, {dirList.name});
    fileNames = cat(1, {fileList.name});
        
    dirListBoxData.numDirEntries = length(dirNames);
    dirListBoxData.dirNames = dirNames;
    dirListBoxData.fileNames = fileNames;
    
	dirListBoxEntries = htmlItem(dirNames, 'blue');
    dirListBoxEntries = [dirListBoxEntries fileNames];
end
function [channelListBoxEntries, channelListBoxData] = prepareChannelListBox(dirListBoxHandle)     
    index_selected = get(dirListBoxHandle,'Value'); %value ue dirListBox Handle
    dirListBoxData = get(dirListBoxHandle, 'userdata');
    
    isdir = index_selected <= dirListBoxData.numDirEntries;
    index_selected(isdir) = [];    
    files_selected = index_selected - dirListBoxData.numDirEntries;
	fileList = dirListBoxData.fileNames(files_selected);
    
    for li = 1 : length(files_selected)
        fprintf('selected %d: %s\n', files_selected(li), fileList{li});
    end
    [allChannelList, allFileList] = loadAllFileInfos(cell2struct(fileList, 'name'));
    
%	setCommon('allFileInfos', allChannelList);
%   setCommon('allFileList', allFileList);
    [channelListBoxEntries, channelListBoxData] = prepareChannelListBox2(allChannelList, allFileList);
end
function [channelListBoxEntries, channelListBoxData] = prepareChannelListBox2(allChannelList, allFileList)
    channelListBoxEntries = []; %{'<HTML>test1</HTML>','<HTML><FONT color="blue">test2</FONT></HTML>','test3'};
    
    channelListBoxData.allChannelList = allChannelList;
    channelListBoxData.allFileList = allFileList;
    channelListBoxData.imap = [];
    channelListBoxData.smap = [];
    channelListBoxData.ischannel = [];
    channelListBoxData.channelsFiltered = [];
%    channelListBoxData.channels = [];
    channelListBoxData.currChannels = [];

    if ~isempty(allChannelList)
        channelsFiltered = applyFilters(allChannelList);
        isvalid = cat(1, channelsFiltered.filterValid);
        channelsFiltered(~isvalid) = [];

        allIDs = cat(1, channelsFiltered.ID);

        lastID = 0;
        ei = 1;
        channelListBoxData.channelsFiltered = channelsFiltered;
%        channelListBoxData.channels = allChannelList; 
        channelListBoxEntries = cell(length(channelsFiltered) + length(unique(allIDs)), 1);
        for i = 1 : length(channelsFiltered)
            currChannel = channelsFiltered(i);
            if currChannel.ID ~= lastID % add filename entry
                channelListBoxEntries{ei} = htmlItem(channelsFiltered(i).filename, 'black');
                channelListBoxData.imap(ei).channelI = find(allIDs == currChannel.ID);
                channelListBoxData.ischannel(ei) = false;                
                ei = ei + 1;
                lastID = currChannel.ID;
            end
            channelListBoxEntries{ei} = htmlItem(['&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;' channelsFiltered(i).fieldPath], 'green');
            channelListBoxData.imap(ei).channelI = i;
            channelListBoxData.ischannel(ei) = true;
            channelListBoxData.smap(i) = ei;
            ei = ei + 1;
        end
    end
end
function ret = getCurrEventHeader(totalEvents)
    ret = sprintf('Meal Summary: Curr View (Total: %d)', totalEvents);
end
function ret = getBulkEventHeader(totalEvents, totalChannels)
    ret = sprintf('Meal Summary: Bulk (Total: %d in %d channels)', totalEvents, totalChannels);
end
function ret = pos2pixel(pos, dim)
    ret = NaN(size(pos));
    ret([1 3]) = pos([1 3])*dim(1);
    ret([2 4]) = pos([2 4])*dim(2);
end
function ret = pixel2pos(pixel, dim)
    ret = pos2pixel(pixel, 1./dim);
end
function logoPos = updateLogoPos(gui)
    logoPos = getLogoPos(gui);
    set(gui.logo.handle, 'Position', logoPos);
end
function logoPos = getLogoPos(gui)
    whichX = gui.logo.whichX; whichY = gui.logo.whichY;
    logoPos = [sum(gui.panelWidthFactor(1:whichX)) sum(gui.panelHeightFactor(1:whichY)) sum(gui.panelWidthFactor(whichX+1)) sum(gui.panelHeightFactor(whichY+1:end))];

    figureDim = [gui.figureWidth gui.figureHeight];
    panelDim = logoPos(3:4) .* figureDim;
    imageDim = gui.logo.imageDim; 
    %dimOff = (panelDim - imageDim)/2;
%    dimOff = [(panelDim(1) - imageDim(1))/2 (panelDim(2) - imageDim(2))*2/3];
    dimOff = [(panelDim(1) - imageDim(1))/2 (panelDim(2) - imageDim(2))/2];

    logoPixel = pos2pixel(logoPos, figureDim);
    logoPos = pixel2pos([logoPixel(1)+dimOff(1) logoPixel(2)+dimOff(2) imageDim(1) imageDim(2)], figureDim);
end
function gui = loadLogo(gui, whichX, whichY)
%I=imread('Expresso_logo_3.tif'); 
I=imread('expressoGreyE.tif'); 

gui.logo.whichX = whichX;
gui.logo.whichY = whichY;
gui.logo.imageDim = [size(I, 2) size(I, 1)];

gui.logo.handle = axes('units','normalized');

gui.logo.pos = updateLogoPos(gui);
uistack(gui.logo.handle,'bottom');

hi = imagesc(I);

set(gui.logo.handle,'handlevisibility','off', 'visible','off');
end
function [pos, gui] = calcGUIpostions(gui)
	pos.hfig = gui.figPos;

    vpi = 1; hpi = 2;

    gui.panelXDelimiter = gui.delimiter * 2;
    pos.dirListBox      = [gui.panelOffsetX(vpi,hpi)+gui.panelXDelimiter              gui.figureHeight-gui.panelOffsetY(vpi,hpi)                                             (gui.panelDimX(vpi,hpi)-gui.panelXDelimiter)/2       gui.panelDimY(vpi,hpi)-(gui.textDimY+gui.editDimY+gui.delimiter*3)];
    pos.dirEditBox      = [pos.dirListBox(1)                                          pos.dirListBox(2)+gui.panelDimY(vpi,hpi)-(gui.textDimY+gui.delimiter*2)-(gui.textDimY) pos.dirListBox(3)-gui.iconDimX                       gui.editDimY];
    pos.dirEditButton   = [pos.dirListBox(1)+pos.dirEditBox(3)-1                      pos.dirEditBox(2)                                                                      gui.iconDimX                                         gui.editDimY];
    pos.dirTextBox      = [pos.dirEditBox(1)                                          pos.dirEditBox(2)+gui.editDimY+gui.delimiter                                           pos.dirEditBox(3)                                    gui.textDimY];
    pos.channelListBox  = [pos.dirListBox(1)+pos.dirListBox(3)                        pos.dirListBox(2)                                                                      pos.dirListBox(3)                                    pos.dirListBox(4)];
    pos.channelEditBox1 = [pos.channelListBox(1)                                      pos.dirEditBox(2)                                                                      (pos.channelListBox(3)-gui.iconDimX-gui.delimiter)/2 pos.dirEditBox(4)];
    pos.channelTextBox1 = [pos.channelEditBox1(1)                                     pos.dirTextBox(2)                                                                      pos.channelEditBox1(3)                               pos.dirTextBox(4)];
    pos.channelEditBox2 = [pos.channelListBox(1)+pos.channelEditBox1(3)+gui.delimiter pos.channelEditBox1(2)                                                                 pos.channelEditBox1(3)                               pos.channelEditBox1(4)];
    pos.channelTextBox2 = [pos.channelEditBox2(1)                                     pos.dirTextBox(2)                                                                      pos.channelEditBox1(3)                               pos.channelTextBox1(4)];        

    
	vpi = 2; hpi = 2;
    
    gui.flexEditBoxX = pos.dirTextBox(3)-gui.iconDimX*4-gui.delimiter*5;
    gui.flexEditBoxXHalf = (gui.flexEditBoxX-gui.delimiter)/2; %gui.buttonDimX;
    gui.checkDimX = pos.dirTextBox(3);
    gui.checkDimXHalf = gui.checkDimX/2; % (dirTextBoxPosition(3)-gui.delimiter)/2;

    pos.navText            = [gui.panelOffsetX(vpi,hpi)+gui.panelXDelimiter               gui.figureHeight-gui.panelOffsetY(vpi,hpi)+gui.panelDimY(vpi,hpi)-gui.textDimY pos.dirTextBox(3)    gui.textDimY];
    pos.navButtonLeft10    = [pos.navText(1)+gui.delimiter                                pos.navText(2)-pos.navText(4)-gui.delimiter                                    gui.iconDimX         gui.iconDimY];
    pos.navButtonLeft1     = [pos.navButtonLeft10(1)+pos.navButtonLeft10(3)+gui.delimiter pos.navButtonLeft10(2)                                                         gui.iconDimX         gui.iconDimY];
    pos.navButtonEdit      = [pos.navButtonLeft1(1)+pos.navButtonLeft1(3)+gui.delimiter   pos.navButtonLeft10(2)+(gui.iconDimY-gui.editDimY)/2                           gui.flexEditBoxXHalf gui.editDimY];
    pos.navButtonTextC     = [pos.navButtonEdit(1)+pos.navButtonEdit(3)                   pos.navButtonLeft10(2)+(gui.iconDimY-gui.textDimY)/2                           gui.delimiter        gui.textDimY];
    pos.navButtonText      = [pos.navButtonTextC(1)+pos.navButtonTextC(3)                 pos.navButtonLeft10(2)+(gui.iconDimY-gui.textDimY)/2                           gui.flexEditBoxXHalf gui.textDimY];
    pos.navButtonRight1    = [pos.navButtonText(1)+pos.navButtonText(3)+gui.delimiter     pos.navButtonLeft10(2)                                                         gui.iconDimX         gui.iconDimY];
    pos.navButtonRight10   = [pos.navButtonRight1(1)+pos.navButtonRight1(3)+gui.delimiter pos.navButtonLeft10(2)                                                         gui.iconDimX         gui.iconDimY];

    pos.zoomText           = [pos.navText(1)                                                pos.navButtonLeft10(2)-pos.navButtonLeft10(4)-gui.delimiter pos.navText(3)       pos.navText(4)];
    pos.zoomButtonMinus3   = [pos.zoomText(1)+gui.delimiter                                 pos.zoomText(2)-pos.zoomText(4)-gui.delimiter               gui.iconDimX         gui.iconDimY];
    pos.zoomButtonMinus2   = [pos.zoomButtonMinus3(1)+pos.zoomButtonMinus3(3)+gui.delimiter pos.zoomButtonMinus3(2)                                     gui.iconDimX         gui.iconDimY];
    pos.zoomButtonEditFrom = [pos.zoomButtonMinus2(1)+pos.zoomButtonMinus2(3)+gui.delimiter pos.zoomButtonMinus3(2)+(gui.iconDimY-gui.editDimY)/2       gui.flexEditBoxXHalf gui.editDimY];
    pos.zoomButtonTextC    = [pos.zoomButtonEditFrom(1)+pos.zoomButtonEditFrom(3)           pos.zoomButtonMinus3(2)+(gui.iconDimY-gui.textDimY)/2       gui.delimiter        gui.textDimY];
    pos.zoomButtonEditTo   = [pos.zoomButtonTextC(1)+pos.zoomButtonTextC(3)                 pos.zoomButtonMinus3(2)+(gui.iconDimY-gui.editDimY)/2       gui.flexEditBoxXHalf gui.editDimY];
    pos.zoomButtonPlus2    = [pos.zoomButtonEditTo(1)+pos.zoomButtonEditTo(3)+gui.delimiter pos.zoomButtonMinus3(2)                                     gui.iconDimX         gui.iconDimY];
    pos.zoomButtonPlus3    = [pos.zoomButtonPlus2(1)+pos.zoomButtonPlus2(3)+gui.delimiter   pos.zoomButtonMinus3(2)                                     gui.iconDimX         gui.iconDimY];
    
	pos.zoomButtonLeft     = [pos.zoomButtonMinus3(1)   pos.zoomButtonMinus3(2)-pos.zoomButtonMinus3(4)-gui.delimiter gui.iconDimX     gui.iconDimY];
    pos.zoomButtonRight    = [pos.zoomButtonPlus3(1)    pos.zoomButtonLeft(2)                                         gui.iconDimX     gui.iconDimY];
    pos.zoomButtonEditPos  = [pos.zoomButtonEditFrom(1) pos.zoomButtonLeft(2)+(gui.iconDimY-gui.editDimY)/2           gui.flexEditBoxX gui.editDimY];
    
    pos.eNavText           = [pos.navText(1)            pos.zoomButtonLeft(2)-pos.zoomButtonLeft(4)-gui.delimiter pos.navText(3)           pos.navText(4)];
    pos.eNavButtonLeft     = [pos.navButtonLeft1(1)     pos.eNavText(2)-pos.eNavText(4)-gui.delimiter             pos.navButtonLeft1(3)    pos.navButtonLeft1(4)];
    pos.eNavEdit           = [pos.zoomButtonEditPos(1)  pos.eNavButtonLeft(2)                                     pos.zoomButtonEditPos(3) pos.zoomButtonEditPos(4)];
    pos.eNavButtonRight    = [pos.navButtonRight1(1)    pos.eNavButtonLeft(2)                                     pos.eNavButtonLeft(3)    pos.eNavButtonLeft(4)];
        
    pos.outText            = [pos.channelListBox(1)                 pos.navText(2)                                                            pos.navText(3) pos.navText(4)];
    pos.outEditDir         = [pos.outText(1)                        pos.outText(2)-pos.outText(4)-gui.delimiter+(gui.iconDimY-gui.editDimY)/2 pos.outText(3) gui.editDimY];
	pos.outButtonDir       = [pos.outEditDir(1)+pos.outEditDir(3)-1 pos.outEditDir(2)                                                         gui.iconDimX   gui.iconDimY];

	pos.optText           = [pos.outText(1)                      pos.outEditDir(2)-gui.iconDimY-gui.delimiter pos.outText(3)    pos.outText(4)];
    pos.optSaveFig        = [pos.optText(1)                      pos.optText(2)-gui.buttonDimY-gui.delimiter  gui.checkDimXHalf gui.buttonDimY];
    pos.optSaveExcel      = [pos.optSaveFig(1)+pos.optSaveFig(3) pos.optSaveFig(2)                            gui.checkDimXHalf pos.optSaveFig(4)];
    pos.optNoOverwrite    = [pos.optSaveFig(1)                   pos.optSaveFig(2)-pos.optSaveFig(4)          gui.checkDimX     pos.optSaveFig(4)];    
    pos.optSaveRaster     = pos.optNoOverwrite;
    pos.optSaveExcelSum   = [pos.optSaveExcel(1)                 pos.optNoOverwrite(2)                        gui.checkDimXHalf pos.optSaveFig(4)];    
    pos.optExitDone       = [pos.optSaveFig(1)                   pos.optSaveRaster(2)-pos.optSaveRaster(4) gui.checkDimXHalf pos.optSaveFig(4)];
    pos.optSaveRasterSum  = pos.optExitDone;
%    pos.optSaveRaster     = [pos.optSaveFig(1)                   pos.optSaveExcelSum(2)-pos.optSaveExcelSum(4) gui.checkDimXHalf pos.optSaveFig(4)];
    pos.optGenFig         = [pos.optSaveExcel(1)                 pos.optExitDone(2)                           gui.checkDimXHalf pos.optSaveFig(4)];
    pos.optNoEmpty        = pos.optGenFig;
%   pos.optNoEmpty        = [pos.optSaveExcel(1)                 pos.optExitDone(2)                           gui.checkDimXHalf pos.optSaveFig(4)];    
	pos.optText2          = [pos.optSaveFig(1)                   pos.optExitDone(2)-pos.optExitDone(4)        gui.checkDimXHalf pos.optSaveFig(4)];
    pos.optDrawLegend     = [pos.optSaveExcel(1)                 pos.optExitDone(2)-pos.optExitDone(4)        gui.checkDimXHalf pos.optSaveFig(4)];    
    
    if gui.figPos(4) > 700
        gui.outButtonDimY = gui.buttonDimY*2;
    else
        gui.outButtonDimY = gui.buttonDimY;
    end
    pos.outButtonSave     = [pos.zoomButtonEditFrom(1)                     pos.optDrawLegend(2)-pos.optDrawLegend(4)-gui.outButtonDimY*1.5-gui.delimiter*3 pos.zoomButtonEditTo(1)+pos.zoomButtonEditTo(3)-pos.zoomButtonEditFrom(1) gui.outButtonDimY];
    pos.outButtonSaveAll  = [pos.outButtonSave(1)+gui.panelDimX(vpi,hpi)/2 pos.outButtonSave(2)                                                            pos.outButtonSave(3)                                                      pos.outButtonSave(4)];
    pos.outSelectAll      = [pos.outButtonSaveAll(1)                       pos.outButtonSaveAll(2)+pos.outButtonSaveAll(4)+gui.delimiter                   pos.outButtonSaveAll(3)                                                   pos.outButtonSaveAll(4)];
    
        
    vpi = 3; hpi = 2;
    
    gui.flexEditBoxX = pos.dirTextBox(3);
    gui.flexEditBoxXHalf = (gui.flexEditBoxX-gui.delimiter/2)/2; %gui.buttonDimX;
    
    pos.eventCurrText      = [pos.dirListBox(1)                                         gui.figureHeight-gui.panelOffsetY(vpi,hpi)+gui.panelDimY(vpi,hpi)-gui.textDimY gui.flexEditBoxX          gui.textDimY];
    pos.eventCurrLatencyT  = [pos.eventCurrText(1)                                      pos.eventCurrText(2)-pos.eventCurrText(4)-gui.delimiter/2                      gui.flexEditBoxXHalf      pos.eventCurrText(4)];
    pos.eventCurrLatencyV  = [pos.eventCurrText(1)+gui.flexEditBoxXHalf+gui.delimiter/2 pos.eventCurrLatencyT(2)                                                       pos.eventCurrLatencyT(3)  pos.eventCurrLatencyT(4)];
    pos.eventCurrDurationT = [pos.eventCurrLatencyT(1)                                  pos.eventCurrLatencyT(2)-pos.eventCurrLatencyT(4)                              pos.eventCurrLatencyT(3)  pos.eventCurrLatencyT(4)];
    pos.eventCurrDurationV = [pos.eventCurrLatencyV(1)                                  pos.eventCurrDurationT(2)                                                      pos.eventCurrDurationT(3) pos.eventCurrDurationT(4)];
    pos.eventCurrVolumeT   = [pos.eventCurrDurationT(1)                                 pos.eventCurrDurationT(2)-pos.eventCurrDurationT(4)                            pos.eventCurrDurationT(3) pos.eventCurrDurationT(4)];
    pos.eventCurrVolumeV   = [pos.eventCurrDurationV(1)                                 pos.eventCurrVolumeT(2)                                                        pos.eventCurrVolumeT(3)   pos.eventCurrVolumeT(4)];        
        
	pos.eventSelectedText      = [pos.channelListBox(1)                                         pos.eventCurrText(2)                                              pos.eventCurrText(3)          pos.eventCurrText(4)];
    pos.eventSelectedLatencyT  = [pos.eventSelectedText(1)                                      pos.eventSelectedText(2)-pos.eventSelectedText(4)-gui.delimiter/2 gui.flexEditBoxXHalf          pos.eventSelectedText(4)];
    pos.eventSelectedLatencyV  = [pos.eventSelectedText(1)+gui.flexEditBoxXHalf+gui.delimiter/2 pos.eventSelectedLatencyT(2)                                      pos.eventSelectedLatencyT(3)  pos.eventSelectedLatencyT(4)];
    pos.eventSelectedDurationT = [pos.eventSelectedLatencyT(1)                                  pos.eventSelectedLatencyT(2)-pos.eventSelectedLatencyT(4)         pos.eventSelectedLatencyT(3)  pos.eventSelectedLatencyT(4)];
    pos.eventSelectedDurationV = [pos.eventSelectedLatencyV(1)                                  pos.eventSelectedDurationT(2)                                     pos.eventSelectedDurationT(3) pos.eventSelectedDurationT(4)];
    pos.eventSelectedVolumeT   = [pos.eventSelectedDurationT(1)                                 pos.eventSelectedDurationT(2)-pos.eventSelectedDurationT(4)       pos.eventSelectedDurationT(3) pos.eventSelectedDurationT(4)];
    pos.eventSelectedVolumeV   = [pos.eventSelectedDurationV(1)                                 pos.eventSelectedVolumeT(2)                                       pos.eventSelectedVolumeT(3)   pos.eventSelectedVolumeT(4)];    

	pos.eventBulkText      = [pos.eventCurrText(1)                                      pos.eventCurrVolumeT(2)-pos.eventCurrVolumeT(4)-gui.delimiter         pos.eventCurrText(3)      pos.eventCurrText(4)];
    pos.eventBulkLatencyT  = [pos.eventBulkText(1)                                      pos.eventBulkText(2)-pos.eventBulkText(4)-gui.delimiter/2             gui.flexEditBoxXHalf      pos.eventBulkText(4)];
    pos.eventBulkLatencyV  = [pos.eventBulkText(1)+gui.flexEditBoxXHalf+gui.delimiter/2 pos.eventBulkLatencyT(2)                                              pos.eventBulkLatencyT(3)  pos.eventBulkLatencyT(4)];
    pos.eventBulkDurationT = [pos.eventBulkLatencyT(1)                                  pos.eventBulkLatencyT(2)-pos.eventBulkLatencyT(4)                     pos.eventBulkLatencyT(3)  pos.eventBulkLatencyT(4)];
    pos.eventBulkDurationV = [pos.eventBulkLatencyV(1)                                  pos.eventBulkDurationT(2)                                             pos.eventBulkDurationT(3) pos.eventBulkDurationT(4)];
    pos.eventBulkVolumeT   = [pos.eventBulkDurationT(1)                                 pos.eventBulkDurationT(2)-pos.eventBulkDurationT(4)                   pos.eventBulkDurationT(3) pos.eventBulkDurationT(4)];
    pos.eventBulkVolumeV   = [pos.eventBulkDurationV(1)                                 pos.eventBulkVolumeT(2)                                               pos.eventBulkVolumeT(3)   pos.eventBulkVolumeT(4)];    
    
	pos.eventBulkText2     = [pos.channelListBox(1)                                                     pos.eventBulkText(2)      pos.eventBulkText(3)      pos.eventBulkText(4)];
    pos.eventBulkInterEvtT = [pos.eventBulkText2(1)                                                     pos.eventBulkLatencyT(2)  pos.eventBulkLatencyT(3)  pos.eventBulkLatencyT(4)];
    pos.eventBulkInterEvtV = [pos.eventBulkText2(1)+(pos.eventBulkLatencyV(1)-pos.eventBulkLatencyT(1)) pos.eventBulkLatencyV(2)  pos.eventBulkLatencyV(3)  pos.eventBulkLatencyV(4)];
    pos.eventBulkEvtPercT  = [pos.eventBulkInterEvtT(1)                                                 pos.eventBulkDurationT(2) pos.eventBulkInterEvtT(3) pos.eventBulkInterEvtT(4)];
    pos.eventBulkEvtPercV  = [pos.eventBulkInterEvtV(1)                                                 pos.eventBulkDurationV(2) pos.eventBulkInterEvtV(3) pos.eventBulkInterEvtV(4)];
	pos.eventBulkSpeedT    = [pos.eventBulkEvtPercT(1)                                                  pos.eventBulkVolumeT(2)   pos.eventBulkEvtPercT(3)  pos.eventBulkEvtPercT(4)];
	pos.eventBulkSpeedV    = [pos.eventBulkEvtPercV(1)                                                  pos.eventBulkVolumeV(2)   pos.eventBulkEvtPercV(3)  pos.eventBulkEvtPercV(4)];

        
	pos.debugMethText = [pos.channelListBox(1)  pos.eventCurrText(2)                              pos.eventCurrText(3) pos.eventCurrText(4)];
    pos.debugMethOpt1 = [pos.debugMethText(1)-1 pos.debugMethText(2)-gui.buttonDimY-gui.delimiter gui.flexEditBoxX     gui.buttonDimY];
    pos.debugMethOpt2 = [pos.debugMethOpt1(1)   pos.debugMethOpt1(2)-pos.debugMethOpt1(4)         pos.debugMethOpt1(3) pos.debugMethOpt1(4)];
    pos.debugMethOpt3 = [pos.debugMethOpt1(1)   pos.debugMethOpt2(2)-pos.debugMethOpt2(4)         pos.debugMethOpt1(3) pos.debugMethOpt1(4)];
    
    pos.debugValText       = [pos.debugMethText(1)                                              pos.debugMethOpt3(2)-pos.debugMethOpt3(4)-gui.delimiter     pos.debugMethText(3)      pos.debugMethText(4)];
    pos.debugValThreshTextOffsY = pos.debugValText(2)-pos.debugValText(4)-gui.delimiter;
    pos.debugValThreshText = [pos.debugValText(1)                                               pos.debugValThreshTextOffsY+(gui.buttonDimY-gui.textDimY)/2 gui.flexEditBoxXHalf      gui.textDimY];
    pos.debugValThreshEdit = [pos.debugValThreshText(1)+pos.debugValThreshText(3)+gui.delimiter pos.debugValThreshTextOffsY+(gui.buttonDimY-gui.editDimY)/2 pos.debugValThreshText(3) gui.editDimY];
    
%	pos.optText2      = [pos.debugValText(1) pos.debugValThreshText(2)-pos.debugValThreshText(4)-gui.delimiter pos.debugValText(3) pos.debugValText(4)];
%    pos.optDrawLegend = [pos.optText2(1)     pos.optText2(2)-gui.buttonDimY-gui.delimiter                  gui.checkDimX       gui.buttonDimY];    
    
    
    [pos.debugMeth, pos.radioPosRAW] = abs2perc([gui.figureWidth, gui.figureHeight], [pos.debugMethOpt3(1) pos.debugMethOpt3(2) pos.debugMethOpt1(1)+pos.debugMethOpt1(3) pos.debugMethOpt1(2)+pos.debugMethOpt1(4)]);
    pos.debugMethOpt1 = adjustRadioPos(pos.radioPosRAW, pos.debugMethOpt1);
    pos.debugMethOpt2 = adjustRadioPos(pos.radioPosRAW, pos.debugMethOpt2);
    pos.debugMethOpt3 = adjustRadioPos(pos.radioPosRAW, pos.debugMethOpt3);
       
    
    function [ret, retRAW] = abs2perc(pdim, l, o, r, u)
        if nargin == 2
            o = l(2);
            r = l(3);
            u = l(4);
            l = l(1);
        end
        ret(1) = l;
        ret(2) = o;
        ret(3) = abs(r-l);
        ret(4) = abs(o-u);
        retRAW = ret;
        ret([1 3]) = ret([1 3]) ./ pdim(1);
        ret([2 4]) = ret([2 4]) ./ pdim(2);
    end
    function ret = adjustRadioPos(radioRAW, opt)
        ret(1) = radioRAW(1) - opt(1);
        ret(2) = -(radioRAW(2) - opt(2));
        ret([3 4]) = opt([3 4]);
    end    
end
function setFigureMinSize(hfig, minSize)
    jFrame = get(handle(hfig), 'JavaFrame');
    try
        jProx = jFrame.fFigureClient.getWindow();
    catch
        jProx = jFrame.fHG1Client.getWindow(); 
    end
    jProx.setMinimumSize(java.awt.Dimension(minSize(1), minSize(2)));    
end
function gui = initGUIConst(hfig, firstInvocation)
    
    screenSize = get(0, 'ScreenSize');
    gui.screenWidth = screenSize(3);
    gui.screenHeight = screenSize(4);    
    
    if nargin >= 2 && firstInvocation
        set(hfig, 'Position', [5 0 gui.screenWidth gui.screenHeight]);    
    end
	gui.figPos = get(hfig, 'Position');
    
	gui.hfig = hfig;

    gui.figureWidth = gui.figPos(3);
    gui.figureHeight = gui.figPos(4);

    gui.defaultPanelHeightFactor = [0.33 0.33 0.34];
    gui.minPanelHeight = [NaN NaN NaN];
    [gui.panelHeightFactor, totalHeight] = getConformFactor(gui.defaultPanelHeightFactor, gui.minPanelHeight, gui.figureHeight);
    
    gui.defaultPanelWidthFactor = [0.60 0.40];
    gui.minPanelWidth = [522 352];
    [gui.panelWidthFactor, totalWidth] = getConformFactor(gui.defaultPanelWidthFactor, gui.minPanelWidth, gui.figureWidth);
    
    gui.minPanelDim = 20;
    setFigureMinSize(hfig, [sum(gui.minPanelWidth(~isnan(gui.minPanelWidth)))+sum(isnan(gui.minPanelWidth)*gui.minPanelDim) max(680, sum(gui.minPanelHeight(~isnan(gui.minPanelHeight)))+sum(isnan(gui.minPanelHeight)*gui.minPanelDim))]);
        
    if totalHeight > gui.figureHeight || totalWidth > gui.figureWidth
        newFigPos = gui.figPos;
        newFigPos(3) = totalWidth;
        newFigPos(4) = totalHeight;
        if newFigPos(1) + newFigPos(3) > gui.screenWidth
            newFigPos(1) = gui.screenWidth - newFigPos(3);
        end
        if newFigPos(2) + newFigPos(4) > gui.screenHeight
            newFigPos(2) = gui.screenHeight - newFigPos(4);
        end
        set(hfig, 'Position', newFigPos);
        gui.figPos = newFigPos;
    end

    gui.numVertPanel = length(gui.panelHeightFactor);
    gui.numHorizPanel = length(gui.panelWidthFactor);

    gui.panelHeightFactor = [0 gui.panelHeightFactor];
    gui.panelWidthFactor = [0 gui.panelWidthFactor];


    gui.panelOffset = 0;%0.05
    %gui.panelWidthFactor = 1-gui.panelOffset;
    gui.sepHeight = 0.0;

    gui.delimiter = 4;
    gui.iconDimX = 20;
    gui.iconDimY = 20;
    gui.buttonDimX = gui.iconDimX * 2 + gui.delimiter;
    gui.buttonDimY = gui.iconDimY;
    gui.textDimY = 16;
    gui.editDimY = 18;
    gui.zoomFactor3 = 3;
    gui.zoomFactor2 = 2;
    gui.zoomShiftFactor1 = 0.5; % scroll half the panel size per button click
        
    for vpi = 1 : gui.numVertPanel
        for hpi = 1 : gui.numHorizPanel
            gui.panelOffsetX(vpi, hpi) = sum(gui.panelWidthFactor(1:hpi)) * gui.figureWidth;
            gui.panelOffsetY(vpi, hpi) = (sum(gui.panelHeightFactor(end-vpi+1:end))-gui.sepHeight) * gui.figureHeight;
            gui.panelDimX(vpi, hpi) = gui.panelWidthFactor(hpi+1) * gui.figureWidth;
            gui.panelDimY(vpi, hpi) = gui.panelHeightFactor(end-vpi+1) * gui.figureHeight;
        end
    end
    
	function [ret, newTotal] = getConformFactor(defaultFac, minVal, totalVal)
        function ret = checkConstraint(vals, minVal)
            ret = false(size(vals));
            ret(isnan(minVal)) = true;
            ret(minVal <= vals) = true;
        end
        newTotal = totalVal;
        noConstraint = isnan(minVal);        
        if all(noConstraint)
            ret = defaultFac;            
        elseif sum(minVal) > totalVal
            vals = minVal;
            vals(noConstraint) = 0;
            newTotal = sum(vals);
            ret = vals ./ newTotal;
        else
            takeMin = false(size(defaultFac));
            ret = defaultFac; NaN(size(defaultFac));
            
            vals = totalVal .* defaultFac;
            
            %loop            
            conform = checkConstraint(vals, minVal);
            while(any(~conform))
            
                takeMin = takeMin | ~conform;
                vals(takeMin) = minVal(takeMin);
    %            restVal = totalVal - vals(takeMin);
                ret(takeMin) = vals(takeMin) ./ totalVal;
                restFac = 1 - sum(ret(takeMin));
    %            restVals = vals(~takeMin);

                ret(~takeMin) = (ret(~takeMin) ./ sum(ret(~takeMin))) .* restFac;
                vals(~takeMin) = ret(~takeMin) .* totalVal;
                conform = checkConstraint(vals, minVal);
            end
            
            
%             for vi = 1 : length(vals)
%                 if flexible(vi)
%                     continue;
%                 end
%                 if vals(vi) < minVal(vi)
%                     vals(vi) = minVal(vi);
%                     totalVal = totalVal - minVal(vi);
%                 end
%             end
        end                
    end
end
function updateSubPanels(gui)
    for vpi = 1 : gui.numVertPanel
        for hpi = 1 : gui.numHorizPanel
            set(gui.panel(gui.numVertPanel-vpi+1,hpi), 'pos', [sum(gui.panelWidthFactor(1:hpi)) sum(gui.panelHeightFactor(1:vpi))-gui.sepHeight gui.panelWidthFactor(hpi+1) gui.panelHeightFactor(vpi+1)]);
        end
    end    
end
function hfig = myInitializeGUI(inputDir, outputDir)
% initialize display
    setCommon('init_complete', false);
    
    [pos, gui] = calcGUIpostions(initGUIConst(figure(), true));
    hfig = gui.hfig;
    
    set(hfig, 'KeyPressFcn',  @KeypressFcn, ...
            'Visible','on', ...
            'ResizeFcn',@ResizeFcn, ...
            'CloseRequestFcn',@CloseRequestFcn, ...
            'WindowButtonDownFcn', @buttonDownHandle, ...
            'WindowButtonUpFcn', @buttonUpHandle, ...
            'Interruptible','Off', ...
            'Name', getCommon('titlePrefix'));
    setCommon('hfig', hfig);    
  
        
    gui.panel = zeros(gui.numVertPanel, gui.numHorizPanel);       

	gui = loadLogo(gui, 1, 2);
%	loadLogo(gui.panelWidthFactor, 1);
    for vpi = 1 : gui.numVertPanel
        for hpi = 1 : gui.numHorizPanel
            gui.panel(vpi, hpi) = subplot(gui.numVertPanel, gui.numHorizPanel, hpi + (vpi-1) * gui.numHorizPanel, 'DrawMode','fast');
            axis off;        
        end
    end
    setCommon('panel', gui.panel);
    updateSubPanels(gui);
    xlimitPanelsData(1).handle = gui.panel(1,1); xlimitPanelsData(1).axesLabel = {'[sec]','[nl]'}; xlimitPanelsData(1).offset = [0 0]; xlimitPanelsData(1).axes = 'xy';
    xlimitPanelsData(2).handle = gui.panel(2,1); xlimitPanelsData(2).axesLabel = {'[sec]','[nl]'}; xlimitPanelsData(2).offset = [0 0]; xlimitPanelsData(2).axes = 'xy';
    xlimitPanelsData(3).handle = gui.panel(3,1); xlimitPanelsData(3).axesLabel = {'','[\sigma]'}; xlimitPanelsData(3).offset = [0 0.03]; xlimitPanelsData(3).axes = 'y';
    setCommon('xlimitPanels', cat(2, xlimitPanelsData.handle));
    setCommon('xlimitPanelsData', xlimitPanelsData);
%    setCommon('xlimitPanels', [gui.panel(1,1) gui.panel(2,1) gui.panel(3,1)]);
%    setCommon('labelInsideOffsets', [0 0; 0 0; 0 0.3]);
%    setCommon('labelInsideAxes', {'xy', 'xy', 'y'});
    setCommon('greenLineHandles', []);    

        
    % subplot 1,2 (file controls)
    vpi = 1; hpi = 2; subplot(gui.panel(vpi,hpi));

    handles.dirTextBox = uicontrol('Style','text', 'String', 'Input Directory', 'HorizontalAlignment', 'center', 'Position', pos.dirTextBox);
    handles.dirEditBox = uicontrol('Style','edit', 'String', inputDir, 'Callback', @dirEditBoxCallback,'BackgroundColor','w', 'HorizontalAlignment', 'left', 'Position', pos.dirEditBox);
    handles.dirEditButton = uicontrol('Style','pushbutton', 'String', '...', 'Callback', @dirEditButtonCallback, 'Position', pos.dirEditButton);
    [dirListBoxEntries, dirListBoxData] = prepareDirListBox(inputDir);
    handles.dirListBox = uicontrol('Style', 'listbox', 'String', dirListBoxEntries, 'Callback', @dirListBoxCallback, 'Min', 1, 'Max', length(dirListBoxEntries), 'BackgroundColor','w', 'Position', pos.dirListBox, 'userdata', dirListBoxData);

    p_filter = getCommon('p_filter');
    handles.channelTextBox1 = uicontrol('Style','text', 'String', 'Filter:  XPxx', 'HorizontalAlignment', 'left', 'Position', pos.channelTextBox1);
    handles.channelEditBox1 = uicontrol('Style','edit', 'String', p_filter.container, 'Callback', @channelEditBox1Callback,'BackgroundColor','w', 'Position', pos.channelEditBox1);
    handles.channelTextBox2 = uicontrol('Style','text', 'String', 'Filter:  channel_x', 'HorizontalAlignment', 'left', 'Position', pos.channelTextBox2);
    handles.channelEditBox2 = uicontrol('Style','edit', 'String', p_filter.channel, 'Callback', @channelEditBox2Callback,'BackgroundColor','w', 'Position', pos.channelEditBox2, 'HorizontalAlignment', 'right');
    [channelListBoxEntries, channelListBoxData] = prepareChannelListBox(handles.dirListBox);
    handles.channelListBox  = uicontrol('Style', 'listbox', 'String', channelListBoxEntries, 'Callback', @channelListBoxCallback, 'Min', 0, 'Max', length(channelListBoxEntries), 'BackgroundColor','w', 'Position', pos.channelListBox, 'userdata', channelListBoxData);
    

    % subplot 2,2 (navigation control / output)
    
    vpi = 2; hpi = 2; subplot(gui.panel(vpi,hpi));        
    
    handles.navText          = uicontrol('Style', 'text',      'String', 'File Navigation', 'HorizontalAlignment', 'center', 'Position', pos.navText);
    handles.navButtonLeft10  = uicontrol('Style','pushbutton', 'String', '-10', 'Callback', @navButtonLeft10Callback, 'Position', pos.navButtonLeft10, 'enable', 'off');
    handles.navButtonLeft1   = uicontrol('Style','pushbutton', 'String', '-1',  'Callback', @navButtonLeft1Callback,  'Position', pos.navButtonLeft1, 'enable', 'off');
    handles.navButtonEdit    = uicontrol('Style','edit', 'String', '',    'Callback', @navButtonEditCallback,  'HorizontalAlignment', 'left', 'BackgroundColor','w', 'Position', pos.navButtonEdit, 'enable', 'off');
    handles.navButtonTextC   = uicontrol('Style','text', 'String', '/',                                         'HorizontalAlignment', 'center', 'Position', pos.navButtonTextC, 'enable', 'off');
    handles.navButtonText    = uicontrol('Style','text', 'String', '',                                         'HorizontalAlignment', 'right', 'Position', pos.navButtonText);
    handles.navButtonRight1  = uicontrol('Style','pushbutton', 'String', '+1',  'Callback', @navButtonRight1Callback,  'Position', pos.navButtonRight1, 'enable', 'off');
    handles.navButtonRight10 = uicontrol('Style','pushbutton', 'String', '+10', 'Callback', @navButtonRight10Callback,  'Position', pos.navButtonRight10, 'enable', 'off');

    handles.zoomText           = uicontrol('Style', 'text',      'String', 'Zoom', 'HorizontalAlignment', 'center', 'Position', pos.zoomText);
    handles.zoomButtonMinus3   = uicontrol('Style','pushbutton', 'String', num2str(-gui.zoomFactor3, '%+d'), 'Callback', @zoomButtonMinus3Callback, 'Position', pos.zoomButtonMinus3, 'enable', 'off');
    handles.zoomButtonMinus2   = uicontrol('Style','pushbutton', 'String', num2str(-gui.zoomFactor2, '%+d'),  'Callback', @zoomButtonMinus2Callback,  'Position', pos.zoomButtonMinus2, 'enable', 'off');
    handles.zoomButtonEditFrom = uicontrol('Style','edit', 'String', '',    'Callback', @zoomButtonEditFromCallback,  'HorizontalAlignment', 'left', 'BackgroundColor','w', 'Position', pos.zoomButtonEditFrom, 'enable', 'off');
    handles.zoomButtonTextC    = uicontrol('Style','text', 'String', '-',                                          'HorizontalAlignment', 'center', 'Position', pos.zoomButtonTextC);
    handles.zoomButtonEditTo   = uicontrol('Style','edit', 'String', '',  'Callback', @zoomButtonEditToCallback,  'HorizontalAlignment', 'right', 'BackgroundColor','w', 'Position', pos.zoomButtonEditTo, 'enable', 'off');
    handles.zoomButtonPlus2    = uicontrol('Style','pushbutton', 'String', num2str(gui.zoomFactor2, '%+d'),  'Callback', @zoomButtonPlus2Callback,  'Position', pos.zoomButtonPlus2, 'enable', 'off');
    handles.zoomButtonPlus3    = uicontrol('Style','pushbutton', 'String', num2str(gui.zoomFactor3, '%+d'), 'Callback', @zoomButtonPlus3Callback,  'Position', pos.zoomButtonPlus3, 'enable', 'off');
    
    handles.zoomButtonEditPos  = uicontrol('Style','edit', 'String', '',  'Callback', @zoomButtonEditPosCallback,  'HorizontalAlignment', 'right', 'BackgroundColor','w', 'Position', pos.zoomButtonEditPos);    
    handles.zoomButtonLeft  = uicontrol('Style','pushbutton', 'String', '<-', 'Callback', @zoomButtonLeftCallback,  'Position', pos.zoomButtonLeft, 'enable', 'off');
    handles.zoomButtonRight = uicontrol('Style','pushbutton', 'String', '->', 'Callback', @zoomButtonRightCallback,  'Position', pos.zoomButtonRight, 'enable', 'off');
    
	handles.eNavText        = uicontrol('Style', 'text',      'String', 'Meal Navigation', 'HorizontalAlignment', 'center', 'Position', pos.eNavText);
    handles.eNavButtonLeft  = uicontrol('Style','pushbutton', 'String', '<',  'Callback', @eNavButtonLeftCallback,  'Position', pos.eNavButtonLeft, 'enable', 'off');
    handles.eNavEdit  = uicontrol('Style','edit', 'String', '',    'Callback', @eNavEditCallback,  'HorizontalAlignment', 'right', 'BackgroundColor','w', 'Position', pos.eNavEdit, 'enable', 'off');
    handles.eNavButtonRight = uicontrol('Style','pushbutton', 'String', '>',  'Callback', @eNavButtonRightCallback,  'Position', pos.eNavButtonRight, 'enable', 'off');
    
    handles.outText      = uicontrol('Style', 'text',      'String', 'Output Directory', 'HorizontalAlignment', 'center', 'Position', pos.outText);
    handles.outEditDir = uicontrol('Style','edit', 'String', outputDir, 'Callback', @outEditDirCallback, 'BackgroundColor','w', 'HorizontalAlignment', 'left', 'Position', pos.outEditDir);
    handles.outButtonDir = uicontrol('Style','pushButton', 'String', '..', 'Callback', @outButtonDirCallback, 'Position', pos.outButtonDir);
    handles.outButtonSave  = uicontrol('Style','pushbutton', 'String', 'Save Curr', 'Callback', @outButtonSaveCallback,  'Position', pos.outButtonSave);
    handles.outButtonSaveAll = uicontrol('Style','pushbutton', 'String', 'Save All', 'Callback', @outButtonSaveAllCallback,  'Position', pos.outButtonSaveAll);
    handles.outSelectAll = uicontrol('Style','pushbutton', 'String', 'Select All', 'Callback', @outSelectAllCallback, 'Position', pos.outSelectAll);

    mode = getCommon('mode');
    handles.optText          = uicontrol('Style', 'text',      'String', 'Batch Processing Options / Plot Options', 'HorizontalAlignment', 'center', 'Position', pos.optText);
    handles.optSaveFig    = uicontrol('Style','checkBox', 'String', 'save figure', 'Callback', @optSaveFigCallback, 'Position', pos.optSaveFig, 'value', mode.saveFIG);
    handles.optSaveExcel  = uicontrol('Style','checkBox', 'String', 'save excel', 'Callback', @optSaveExcelCallback, 'Position', pos.optSaveExcel, 'value', mode.saveXLS);
    handles.optNoEmpty  = uicontrol('Style','checkBox', 'String', 'skip empty input', 'Callback', @optNoEmptyCallback, 'Position', pos.optNoEmpty, 'value', mode.noEmpty);
    handles.optText2      = uicontrol('Style','text',     'String', '', 'HorizontalAlignment', 'center', 'Position', pos.optText2);
	handles.optDrawLegend = uicontrol('Style','checkBox', 'String', htmlItem('insert legend', 'black'), 'Callback', @optDrawLegendCallback, 'Position', pos.optDrawLegend, 'value', mode.drawLegend);
    
    handles.optSaveRaster  = uicontrol('Style','checkBox', 'String', 'save raster', 'Callback', @optSaveRasterCallback, 'Position', pos.optSaveRaster, 'value', mode.saveRaster);
    handles.optSaveRasterSum = uicontrol('Style','checkBox', 'String', 'save raster sum', 'Callback', @optSaveRasterSumCallback, 'Position', pos.optSaveRasterSum, 'value', mode.saveRasterSum);
    handles.optSaveExcelSum = uicontrol('Style','checkBox', 'String', 'save excel sum', 'Callback', @optSaveExcelSumCallback, 'Position', pos.optSaveExcelSum, 'value', mode.saveExcelSum);    
        
    % obsolete options begin
    handles.optNoOverwrite  = uicontrol('Style','checkBox', 'String', 'skip existing data', 'Callback', @optNoOverwriteCallback, 'Position', pos.optNoOverwrite, 'value', mode.noOverwrite, 'visible','off');
    handles.optExitDone = uicontrol('Style','checkBox', 'String', 'exit when done', 'Callback', @optExitDoneCallback, 'Position', pos.optExitDone, 'value', mode.exitDone, 'visible', 'off');
    handles.optGenFig = uicontrol('Style','checkBox', 'String', 'open figures', 'Callback', @optGenFigCallback, 'Position', pos.optGenFig, 'value', mode.generateFIG, 'visible', 'off');    
    % obsolete options end
    

    vpi = 3; hpi = 2; subplot(gui.panel(vpi,hpi));        
    
    handles.eventCurrText = uicontrol('Style', 'text', 'String', getCurrEventHeader(0), 'HorizontalAlignment', 'center', 'Position', pos.eventCurrText);
    handles.eventCurrLatencyT = uicontrol('Style', 'text', 'String', '  Latency to 1st: ', 'HorizontalAlignment', 'right', 'Position', pos.eventCurrLatencyT);
    handles.eventCurrLatencyV = uicontrol('Style', 'text', 'String', '', 'HorizontalAlignment', 'left', 'Position', pos.eventCurrLatencyV);
    handles.eventCurrDurationT = uicontrol('Style', 'text', 'String', '  Total Duration: ', 'HorizontalAlignment', 'right', 'Position', pos.eventCurrDurationT);
    handles.eventCurrDurationV = uicontrol('Style', 'text', 'String', '', 'HorizontalAlignment', 'left', 'Position', pos.eventCurrDurationV);    
    handles.eventCurrVolumeT = uicontrol('Style', 'text', 'String', '  Total Volume: ', 'HorizontalAlignment', 'right', 'Position', pos.eventCurrVolumeT);
    handles.eventCurrVolumeV = uicontrol('Style', 'text', 'String', '', 'HorizontalAlignment', 'left', 'Position', pos.eventCurrVolumeV);

    handles.eventSelectedText = uicontrol('Style', 'text', 'String', 'Meal Summary: Selected', 'HorizontalAlignment', 'center', 'Position', pos.eventSelectedText, 'visible', 'off');
%    handles.eventTotalT = uicontrol('Style', 'text', 'String', '  Total: ', 'HorizontalAlignment', 'right', 'Position', eventCurrTotalT);
%    handles.eventTotalV = uicontrol('Style', 'text', 'String', '', 'HorizontalAlignment', 'left', 'Position', eventCurrTotalV);
    handles.eventSelectedLatencyT = uicontrol('Style', 'text', 'String', '  StartTime: ', 'HorizontalAlignment', 'right', 'Position', pos.eventSelectedLatencyT, 'visible', 'off');
    handles.eventSelectedLatencyV = uicontrol('Style', 'text', 'String', '', 'HorizontalAlignment', 'left', 'Position', pos.eventSelectedLatencyV, 'visible', 'off');
    handles.eventSelectedDurationT = uicontrol('Style', 'text', 'String', '  Duration: ', 'HorizontalAlignment', 'right', 'Position', pos.eventSelectedDurationT, 'visible', 'off');
    handles.eventSelectedDurationV = uicontrol('Style', 'text', 'String', '', 'HorizontalAlignment', 'left', 'Position', pos.eventSelectedDurationV, 'visible', 'off');    
    handles.eventSelectedVolumeT = uicontrol('Style', 'text', 'String', '  Volume: ', 'HorizontalAlignment', 'right', 'Position', pos.eventSelectedVolumeT, 'visible', 'off');
    handles.eventSelectedVolumeV = uicontrol('Style', 'text', 'String', '', 'HorizontalAlignment', 'left', 'Position', pos.eventSelectedVolumeV, 'visible', 'off');
    
    handles.eventBulkText = uicontrol('Style', 'text', 'String', getBulkEventHeader(0,0), 'HorizontalAlignment', 'center', 'Position', pos.eventBulkText, 'visible', 'off');
    handles.eventBulkLatencyT = uicontrol('Style', 'text', 'String', '  Avg Latency: ', 'HorizontalAlignment', 'right', 'Position', pos.eventBulkLatencyT, 'visible', 'off');
    handles.eventBulkLatencyV = uicontrol('Style', 'text', 'String', '', 'HorizontalAlignment', 'left', 'Position', pos.eventBulkLatencyV, 'visible', 'off');
    handles.eventBulkDurationT = uicontrol('Style', 'text', 'String', '  Total Duration: ', 'HorizontalAlignment', 'right', 'Position', pos.eventBulkDurationT, 'visible', 'off');
    handles.eventBulkDurationV = uicontrol('Style', 'text', 'String', '', 'HorizontalAlignment', 'left', 'Position', pos.eventBulkDurationV, 'visible', 'off');    
    handles.eventBulkVolumeT = uicontrol('Style', 'text', 'String', '  Total Volume: ', 'HorizontalAlignment', 'right', 'Position', pos.eventBulkVolumeT, 'visible', 'off');
    handles.eventBulkVolumeV = uicontrol('Style', 'text', 'String', '', 'HorizontalAlignment', 'left', 'Position', pos.eventBulkVolumeV, 'visible', 'off');    
    
    handles.eventBulkText2 = uicontrol('Style', 'text', 'String', 'Meal Summary: Bulk (continued)', 'HorizontalAlignment', 'center', 'Position', pos.eventBulkText2, 'visible', 'off');
    handles.eventBulkInterEvtT = uicontrol('Style', 'text', 'String', '  Avg Inter Meal: ', 'HorizontalAlignment', 'right', 'Position', pos.eventBulkInterEvtT, 'visible', 'off');
    handles.eventBulkInterEvtV = uicontrol('Style', 'text', 'String', '', 'HorizontalAlignment', 'left', 'Position', pos.eventBulkInterEvtV, 'visible', 'off');
    handles.eventBulkEvtPercT = uicontrol('Style', 'text', 'String', '  Feeding %: ', 'HorizontalAlignment', 'right', 'Position', pos.eventBulkEvtPercT, 'visible', 'off');
    handles.eventBulkEvtPercV = uicontrol('Style', 'text', 'String', '', 'HorizontalAlignment', 'left', 'Position', pos.eventBulkEvtPercV, 'visible', 'off');    
    handles.eventBulkSpeedT = uicontrol('Style', 'text', 'String', '  Avg Meal Speed: ', 'HorizontalAlignment', 'right', 'Position', pos.eventBulkSpeedT, 'visible', 'off');
    handles.eventBulkSpeedV = uicontrol('Style', 'text', 'String', '', 'HorizontalAlignment', 'left', 'Position', pos.eventBulkSpeedV, 'visible', 'off');    
    
    
    % obsolete debug elements begin
    handles.debugValText = uicontrol('Style', 'text', 'String', 'Debug: Values', 'HorizontalAlignment', 'center', 'Position', pos.debugValText, 'visible', 'off');
    handles.debugValThreshText = uicontrol('Style', 'text', 'String', 'Threshold: ', 'HorizontalAlignment', 'left', 'Position', pos.debugValThreshText, 'visible', 'off');
	devfakt = getCommon('devfakt');
    handles.debugValThreshEdit = uicontrol('Style','edit', 'String', num2str(devfakt.event), 'Callback', @debugValThreshEditCallback, 'BackgroundColor','w', 'HorizontalAlignment', 'right', 'Position', pos.debugValThreshEdit, 'visible', 'off');
        
    debugMeth = uibuttongroup('Position', pos.debugMeth, 'visible', 'off');
	handles.debugMethText = uicontrol('Style', 'text',                       'String', 'Debug: Methods', 'HorizontalAlignment', 'center', 'Position', pos.debugMethText, 'visible', 'off');
    handles.debugMethOpt1 = uicontrol('Style', 'radio', 'parent', debugMeth, 'String', 'Threshold Learned Variance', 'HorizontalAlignment', 'left', 'Position', pos.debugMethOpt1, 'HandleVisibility','off', 'visible', 'off');
    handles.debugMethOpt2 = uicontrol('Style', 'radio', 'parent', debugMeth, 'String', 'Threshold Fixed Volume', 'HorizontalAlignment', 'left', 'Position', pos.debugMethOpt2, 'HandleVisibility','off', 'enable', 'off', 'visible', 'off');
    handles.debugMethOpt3 = uicontrol('Style', 'radio', 'parent', debugMeth, 'String', 'SVM Classifier', 'HorizontalAlignment', 'left', 'Position', pos.debugMethOpt3, 'HandleVisibility','off', 'enable', 'off', 'visible', 'off');
    set(debugMeth,'SelectionChangeFcn',@debugMethCallback);
    set(debugMeth,'SelectedObject', handles.debugMethOpt1);  % empty for no selection
    % obsolete debug elements end


%    setCommon('handles', handles);

    gui.handles = handles;
    gui.pos = pos;
    set(hfig, 'userdata', gui);    
    guidata(hfig, gui);
    
%    setCommon('gui', gui);

    %SetupTimer(hfig);
    %g = guidata(hfig);
    setCommon('init_complete', true);
end

    function hdf5Data = myHdf5ToStruct(fileName, p_loadData)
    % Usage: data = hdf5ToStruct('my_file.hdf5')
    %        data = hdf5ToStruct('my_file.hdf5', true)
    %
    % Reads hdf5 file and recursively converts it to a matlab structure. Each 
    % dataset and group become a field in the structure, sub-groups become
    % sub-fields, etc. 
    %
    % when p_loadData == 1 fields will contain their data from file
    % when p_loadData == 0 fields will contain their hdf5-fieldpath (for loading)
    % -------------------------------------------------------------------------
    if nargin < 2
        p_loadData = false;
    end

    info = hdf5info(fileName);
    group = info.GroupHierarchy;
    hdf5Data = groupToStruct(group, fileName);

        function groupStruct = groupToStruct(group, fileName)
            groupStruct = struct;
            if isfield(group,'Datasets')
                for i = 1:length(group.Datasets)    
                    name = group.Datasets(i).Name;
                    nameEnd = getNameEnd(name);
                    if p_loadData
                        groupStruct.(nameEnd) = hdf5read(fileName,name);     
                    else
                        groupStruct.(nameEnd) = name;
                    end
                end
            end
            if isfield(group, 'Groups')
                for i = 1:length(group.Groups)
                    name = group.Groups(i).Name;
                    nameEnd = getNameEnd(name);
                    groupStruct.(nameEnd) = groupToStruct(group.Groups(i),fileName);
                end
            end

            function nameEnd = getNameEnd(name)
                nameParts = regexp(name(2:end),'/','split');
                nameEnd = nameParts{end};
            end
        end
    end
    function hout = labelInsideAxes(axesHandle, ax, axesLabel, offset)
        
        if nargin < 4
            offset = [0 0];
            if nargin < 2
                ax = 'xy';
            end
        end
        if nargin < 3
            axesLabel = {'',''};
        end
        
        textHandles = get(axesHandle, 'userdata');
        if ~isempty(textHandles)           
            delete(textHandles.numbers);
            delete(textHandles.units);
        end
        
        hout = [];        
        hLabel = [];
        for xi = 1 : length(ax)
            axi = ax(xi);
            
            %get axis properties
            tick=get(axesHandle,[axi 'tick']);
            set(axesHandle,[axi 'TickLabelMode'],'auto')            
            ticklabel=get(axesHandle,[axi 'ticklabel']);
            set(axesHandle,[axi 'TickLabel'],[]);
            lim=get(axesHandle,[axi 'lim']);
            %reformat labels
            ticklabel=cellstr(ticklabel);
            
            %normalized x positions
            ticknorm=(tick-lim(1))/(lim(2)-lim(1));
%            ind=ticknorm ~=0 & ticknorm~=1; %don't print numbers at the corners.
            ind=ticknorm >= 0.015 & ticknorm <= 0.97; %don't print numbers at the corners.
            
            % add axes label
%            [ticklabel(ind), hLabel] = addAxesLabel(ticklabel(ind), axesLabel, axi);
            hLabelNew = addAxesLabel(ticklabel(ind), axesLabel, axi, ticknorm(ind));
            hLabel = [hLabel; hLabelNew];

            if axi == 'x'
                hNew=text(ticknorm(ind)+offset(1),repmat(.05,1,sum(ind))+offset(2),ticklabel(ind),...
                    'units','normalized','horizontalalignment','center', 'fontsize',get(axesHandle,'fontsize'),'parent',axesHandle);
            elseif axi == 'y'
                hNew=text(repmat(.015,1,sum(ind))+offset(1),ticknorm(ind)+offset(2),ticklabel(ind),...
                    'units','normalized','horizontalalignment','left','fontsize',get(axesHandle,'fontsize'),'parent',axesHandle);            
            end
%            set(axesHandle,[axi 'ticklabel'],[]);
            hout = [hout; hNew];            
        end
        textHandles.numbers = hout;
        textHandles.units = hLabel;
        set(axesHandle, 'userdata', textHandles);
        
        function hLabel = addAxesLabel(ticklabel, axesLabel, axi, ticknorm)
            if axi == 'x'
%                ticklabel{1} = [axesLabel{1} '\newline' ticklabel{1}]; % [ticklabel{1} ' ' axesLabel{1}]; % axesLabel{abs('x'-axi)+1}
%                ticklabel{1} = [ticklabel{1} ' ' axesLabel{1}]; % axesLabel{abs('x'-axi)+1}
                hLabel = text(mean(ticknorm(1:2)+offset(1)),.05+offset(2), axesLabel{1}, ...
                    'units','normalized','horizontalalignment','center', 'fontsize',get(axesHandle,'fontsize'),'parent',axesHandle,...
                    'color',[0.7 0.7 0.7]);
            elseif axi == 'y'
%                ticklabel{end} = [ticklabel{1} ' ' axesLabel{2}];% axesLabel{abs('x'-axi)+1}
%                hLabel = [];
                
                hLabel = text(0.015+offset(1)+0.011*length(ticklabel{end}), ticknorm(end)+offset(2), axesLabel{2}, ...
                    'units','normalized','horizontalalignment','left', 'fontsize',get(axesHandle,'fontsize'),'parent',axesHandle,...
                    'color',[0.7 0.7 0.7]);
            end
        end
    end
end
