function out = clhe(varargin)

%-----------------------------------------------------------------------------

matlab.images.internal.errorIfgpuArray(varargin{:});

[I, selectedRange, fullRange, numTiles, dimTile, clipLimit, numBins, ...
    noPadRect, distribution, alpha, int16ClassChange] = parseInputs(varargin{:});

tileMappings = makeTileMappings(I, numTiles, dimTile, numBins, clipLimit, ...
    selectedRange, fullRange, distribution, alpha);

%Synthesize the output image based on the individual tile mappings.
out = makeClaheImage(I, tileMappings, numTiles, selectedRange, numBins,...
    dimTile);

if int16ClassChange
    % Change uint16 back to int16 so output has same class as input.
    out = images.internal.builtins.uint16toint16(out);
end

if ~isempty(noPadRect) %do we need to remove padding?
    out = out(noPadRect.ulRow:noPadRect.lrRow, ...
        noPadRect.ulCol:noPadRect.lrCol);
end

%-----------------------------------------------------------------------------

function tileMappings = ...
    makeTileMappings(I, numTiles, dimTile, numBins, clipLimit,...
    selectedRange, fullRange, distribution, alpha)

numPixInTile = prod(dimTile);

tileMappings = cell(numTiles);

% extract and process each tile
imgCol = 1;
for col=1:numTiles(2),
    imgRow = 1;
    for row=1:numTiles(1),

        tile = I(imgRow:imgRow+dimTile(1)-1,imgCol:imgCol+dimTile(2)-1);

        % for speed, call MEX file directly thus avoiding costly
        % input parsing of imhist
        tileHist = images.internal.builtins.imhistc(tile, numBins, 1, fullRange(2));

        tileHist = clipHistogram(tileHist, clipLimit, numBins);

        tileMapping = makeMapping(tileHist, selectedRange, fullRange, ...
            numPixInTile, distribution, alpha);

        % assemble individual tile mappings by storing them in a cell array;
        tileMappings{row,col} = tileMapping;

        imgRow = imgRow + dimTile(1);
    end
    imgCol = imgCol + dimTile(2); % move to the next column of tiles
end

%-----------------------------------------------------------------------------
% Calculate the equalized lookup table (mapping) based on cumulating the input
% histogram.  Note: lookup table is rescaled in the selectedRange [Min..Max].

function mapping = makeMapping(imgHist, selectedRange, fullRange, ...
    numPixInTile, distribution, alpha)

histSum = cumsum(imgHist);
valSpread  = selectedRange(2) - selectedRange(1);

switch distribution
    case 'uniform',
        scale =  valSpread/numPixInTile;
        mapping = min(selectedRange(1) + histSum*scale,...
            selectedRange(2)); %limit to max

    case 'rayleigh', % suitable for underwater imagery
        % pdf = (t./alpha^2).*exp(-t.^2/(2*alpha^2))*U(t)
        % cdf = 1-exp(-t.^2./(2*alpha^2))
        hconst = 2*alpha^2;
        vmax = 1 - exp(-1/hconst);
        val = vmax*(histSum/numPixInTile);
        val(val>=1) = 1-eps; % avoid log(0)
        temp = sqrt(-hconst*log(1-val));
        mapping = min(selectedRange(1)+temp*valSpread,...
            selectedRange(2)); %limit to max

    case 'exponential',
        % pdf = alpha*exp(-alpha*t)*U(t)
        % cdf = 1-exp(-alpha*t)
        vmax = 1 - exp(-alpha);
        val = (vmax*histSum/numPixInTile);
        val(val>=1) = 1-eps;
        temp = -1/alpha*log(1-val);
        mapping = min(selectedRange(1)+temp*valSpread,selectedRange(2));

    otherwise,
        error(message('images:adapthisteq:distributionType')) %should never get here

end

%rescale the result to be between 0 and 1 for later use by the GRAYXFORMMEX
%private mex function
mapping = mapping/fullRange(2);

%-----------------------------------------------------------------------------
% This function clips the histogram according to the clipLimit and
% redistributes clipped pixels across bins below the clipLimit

function imgHist = clipHistogram(imgHist, clipLimit, numBins)

% total number of pixels overflowing clip limit in each bin
totalExcess = sum(max(imgHist - clipLimit,0));

% clip the histogram and redistribute the excess pixels in each bin
avgBinIncr = floor(totalExcess/numBins);
upperLimit = clipLimit - avgBinIncr; % bins larger than this will be
% set to clipLimit

% this loop should speed up the operation by putting multiple pixels
% into the "obvious" places first
for k=1:numBins
    if imgHist(k) > clipLimit
        imgHist(k) = clipLimit;
    else
        if imgHist(k) > upperLimit % high bin count
            totalExcess = totalExcess - (clipLimit - imgHist(k));
            imgHist(k) = clipLimit;
        else
            totalExcess = totalExcess - avgBinIncr;
            imgHist(k) = imgHist(k) + avgBinIncr;
        end
    end
end

% this loops redistributes the remaining pixels, one pixel at a time
k = 1;
while (totalExcess ~= 0)
    %keep increasing the step as fewer and fewer pixels remain for
    %the redistribution (spread them evenly)
    stepSize = max(floor(numBins/totalExcess),1);
    for m=k:stepSize:numBins
        if imgHist(m) < clipLimit
            imgHist(m) = imgHist(m)+1;
            totalExcess = totalExcess - 1; %reduce excess
            if totalExcess == 0
                break;
            end
        end
    end

    k = k+1; %prevent from always placing the pixels in bin #1
    if k > numBins % start over if numBins was reached
        k = 1;
    end
end

%-----------------------------------------------------------------------------
% This function interpolates between neighboring tile mappings to produce a
% new mapping in order to remove artificially induced tile borders.
% Otherwise, these borders would become quite visible.  The resulting
% mapping is applied to the input image thus producing a CLAHE processed
% image.

function claheI = makeClaheImage(I, tileMappings, numTiles, selectedRange,...
    numBins, dimTile)

%initialize the output image to zeros (preserve the class of the input image)
claheI = I;
claheI(:) = 0;

%compute the LUT for looking up original image values in the tile mappings,
%which we created earlier
if ~(isa(I,'double') || isa(I,'single'))
    k = selectedRange(1)+1 : selectedRange(2)+1;
    aLut = zeros(length(k),1);
    aLut(k) = (k-1)-selectedRange(1);
    aLut = aLut/(selectedRange(2)-selectedRange(1));
else
    % remap from 0..1 to 0..numBins-1
    if numBins ~= 1
        binStep = 1/(numBins-1);
        start = ceil(selectedRange(1)/binStep);
        stop  = floor(selectedRange(2)/binStep);
        k = start+1:stop+1;
        aLut(k) = 0:1/(length(k)-1):1;
    else
        aLut(1) = 0; %in case someone specifies numBins = 1, which is just silly
    end
end

imgTileRow=1;
for k=1:numTiles(1)+1

            imgTileNumRows = dimTile(1);
            mapTileRows = [k, k]; %[upperRow lowerRow]


    % loop over columns of the tileMappings cell array
    imgTileCol=1;
    for l=1:numTiles(2)+1

                imgTileNumCols = dimTile(2);
                mapTileCols = [l, l]; % right left


        % Extract four tile mappings
        ulMapTile = tileMappings{mapTileRows(1), mapTileCols(1)};
        %      urMapTile = tileMappings{mapTileRows(1), mapTileCols(2)};
        %     blMapTile = tileMappings{mapTileRows(2), mapTileCols(1)};
        %     brMapTile = tileMappings{mapTileRows(2), mapTileCols(2)};

        normFactor = imgTileNumRows*imgTileNumCols; %normalization factor
        imgTileIdx = {imgTileRow:imgTileRow+imgTileNumRows-1, ...
            imgTileCol:imgTileCol+imgTileNumCols-1};

        imgPixVals = images.internal.builtins.grayxform(I(imgTileIdx{1},imgTileIdx{2}), aLut);

        % calculate the weights used for linear interpolation between the
        % four mappings
        %     rowW = repmat((0:imgTileNumRows-1)',1,imgTileNumCols);
        %     colW = repmat(0:imgTileNumCols-1,imgTileNumRows,1);
        %     rowRevW = repmat((imgTileNumRows:-1:1)',1,imgTileNumCols);
        %     colRevW = repmat(imgTileNumCols:-1:1,imgTileNumRows,1);

        claheI(imgTileIdx{1}, imgTileIdx{2}) = ...
            double(images.internal.builtins.grayxform(imgPixVals, ulMapTile));

        imgTileCol = imgTileCol + imgTileNumCols;
    end %over tile cols
    imgTileRow = imgTileRow + imgTileNumRows;
end %over tile rows

%-----------------------------------------------------------------------------

function [I, selectedRange, fullRange, numTiles, dimTile, clipLimit,...
    numBins, noPadRect, distribution, alpha, ...
    int16ClassChange] = parseInputs(varargin)

narginchk(1,13);

I = varargin{1};
validateattributes(I, {'uint8', 'uint16', 'double', 'int16', 'single'}, ...
    {'real', '2d', 'nonsparse', 'nonempty'}, ...
    mfilename, 'I', 1);

% convert int16 to uint16
if isa(I,'int16')
    I = images.internal.builtins.int16touint16(I);
    int16ClassChange = true;
else
    int16ClassChange = false;
end

if any(size(I) < 2)
    error(message('images:adapthisteq:inputImageTooSmall'))
end

%Other options
%%%%%%%%%%%%%%

%Set the defaults
distribution = 'uniform';
alpha   = 0.4;

if isa(I, 'double') || isa(I,'single')
    fullRange = [0 1];
else
    fullRange(1) = I(1);         %copy class of the input image
    fullRange(1:2) = [-Inf Inf]; %will be clipped to min and max
    fullRange = double(fullRange);
end

selectedRange   = fullRange;

%Set the default to 256 bins regardless of the data type;
%the user can override this value at any time
numBins = 256;
normClipLimit = 0.01;
numTiles = [8 8];

checkAlpha = false;

validStrings = {'NumTiles','ClipLimit','NBins','Distribution',...
    'Alpha','Range'};

if nargin > 1
    done = false;

    idx = 2;
    while ~done
        input = varargin{idx};
        inputStr = validatestring(input, validStrings,mfilename,'PARAM',idx);

        idx = idx+1; %advance index to point to the VAL portion of the input

        if idx > nargin
            error(message('images:adapthisteq:missingValue', inputStr))
        end

        switch inputStr

            case 'NumTiles'
                numTiles = varargin{idx};
                validateattributes(numTiles, {'double'}, {'real', 'vector', ...
                    'integer', 'finite','positive'},...
                    mfilename, inputStr, idx);

                if (any(size(numTiles) ~= [1,2]))
                    error(message('images:adapthisteq:invalidNumTilesVector', inputStr))
                end

                if any(numTiles < 2)
                    error(message('images:adapthisteq:invalidNumTilesValue', inputStr))
                end

            case 'ClipLimit'
                normClipLimit = varargin{idx};
                validateattributes(normClipLimit, {'double'}, ...
                    {'scalar','real','nonnegative'},...
                    mfilename, inputStr, idx);

                if normClipLimit > 1
                    error(message('images:adapthisteq:invalidClipLimit', inputStr))
                end

            case 'NBins'
                numBins = varargin{idx};
                validateattributes(numBins, {'double'}, {'scalar','real','integer',...
                    'positive'}, mfilename, 'NBins', idx);

            case 'Distribution'
                validDist = {'rayleigh','exponential','uniform'};
                distribution = validatestring(varargin{idx}, validDist, mfilename,...
                    'Distribution', idx);

            case 'Alpha'
                alpha = varargin{idx};
                validateattributes(alpha, {'double'},{'scalar','real',...
                    'nonnan','positive','finite'},...
                    mfilename, 'Alpha',idx);
                checkAlpha = true;

            case 'Range'
                validRangeStrings = {'original','full'};
                rangeStr = validatestring(varargin{idx}, validRangeStrings,mfilename,...
                    'Range',idx);

                if strmatch(rangeStr,'original')
                    selectedRange = double([min(I(:)), max(I(:))]);
                end

            otherwise
                error(message('images:adapthisteq:inputString')) %should never get here
        end

        if idx >= nargin
            done = true;
        end

        idx=idx+1;
    end
end


% Pre-process the inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%

dimI = size(I);
dimTile = dimI ./ numTiles;

%check if tile size is reasonable
if any(dimTile < 1)
    error(message('images:adapthisteq:inputImageTooSmallToSplit', num2str( numTiles )))
end

if checkAlpha
    if strcmp(distribution,'uniform')
        error(message('images:adapthisteq:alphaShouldNotBeSpecified', distribution))
    end
end

%check if the image needs to be padded; pad if necessary;
%padding occurs if any dimension of a single tile is an odd number
%and/or when image dimensions are not divisible by the selected
%number of tiles
rowDiv  = mod(dimI(1),numTiles(1)) == 0;
colDiv  = mod(dimI(2),numTiles(2)) == 0;

if rowDiv && colDiv
    rowEven = mod(dimTile(1),2) == 0;
    colEven = mod(dimTile(2),2) == 0;
end

noPadRect = [];
if  ~(rowDiv && colDiv && rowEven && colEven)
    padRow = 0;
    padCol = 0;

    if ~rowDiv
        rowTileDim = floor(dimI(1)/numTiles(1)) + 1;
        padRow = rowTileDim*numTiles(1) - dimI(1);
    else
        rowTileDim = dimI(1)/numTiles(1);
    end

    if ~colDiv
        colTileDim = floor(dimI(2)/numTiles(2)) + 1;
        padCol = colTileDim*numTiles(2) - dimI(2);
    else
        colTileDim = dimI(2)/numTiles(2);
    end

    %check if tile dimensions are even numbers
    rowEven = mod(rowTileDim,2) == 0;
    colEven = mod(colTileDim,2) == 0;

    if ~rowEven
        padRow = padRow+numTiles(1);
    end

    if ~colEven
        padCol = padCol+numTiles(2);
    end

    padRowPre  = floor(padRow/2);
    padRowPost = ceil(padRow/2);
    padColPre  = floor(padCol/2);
    padColPost = ceil(padCol/2);

    I = padarray(I,[padRowPre  padColPre ],'symmetric','pre');
    I = padarray(I,[padRowPost padColPost],'symmetric','post');

    %UL corner (Row, Col), LR corner (Row, Col)
    noPadRect.ulRow = padRowPre+1;
    noPadRect.ulCol = padColPre+1;
    noPadRect.lrRow = padRowPre+dimI(1);
    noPadRect.lrCol = padColPre+dimI(2);
end

%redefine this variable to include the padding
dimI = size(I);

%size of the single tile
dimTile = dimI ./ numTiles;

numPixInTile = prod(dimTile);
minClipLimit = ceil(numPixInTile/numBins);
clipLimit = minClipLimit + round(normClipLimit*(numPixInTile-minClipLimit));

%-----------------------------------------------------------------------------
