function [RDKcoordinates,dots] = hb_getRDK( dots )
% Usage : RDKframes = hb_getRDK(dots)
% 
%  EXAMPLE USE
% % dots = [];
% % dots.params.fieldSize = [500 500];
% % dots.params.fieldShape = 1; % 0~1 : Spare factor, 1~ : Square-like field 
% % dots.params.nDots = 100; % number of dots
% % dots.params.coherence = .5; % 0~1 : Proportion of coherent dot
% % dots.params.speed = 100; % Pixel per second
% % dots.params.duration = 2; %
% % dots.params.ori = 0; % 0~360 : Orientation
% % dots.params.frameRate = 60;
% % dots.params.dotLifeOption = 'N'; % 'N' for normal dist, 'U' for uniform dist
% % dots.params.dotLife = dots.params.duration * dots.params.frameRate; % frames, Mean of dots' life
% % nFrames = dots.params.frameRate * dots.params.duration;
% % [RDKcoordinates,dots] = hb_getRDK(dots);
% % for frameIdx = 1:nFrames
% %     plot( RDKcoordinates(:, 1, frameIdx),RDKcoordinates(:, 2, frameIdx), 'kp' );
% %     axis([-1 1 -1 1]*.75*dots.circField.size(1));
% %     axis xy;
% %     drawnow;
% % end
% Written by Hio-Been Han, hiobeen.han@seoutech.ac.kr, 2016-08-01
% Simple alternative for VCRDM Toolbox from Shadlen Lab (www.shadlenlab.columbia.edu) 
% 쉐들렌 코드가 구려서 직접 만들었음

warning off;
if nargin < 1
    disp(['No input argument detected. Make with default setting'])
    dots = [];
    dots.params.fieldSize = [500 500]; % Unit : Pixel
    dots.params.fieldShape = 1; % 0~1 : Spare factor, 1~ : Square-like field
    dots.params.nDots = 100; % number of dots
    dots.params.coherence = .5; % 0~1 : Proportion of coherent dot
    
    dots.params.speed = 100; % Unit : Pixel per second
    dots.params.duration = 2; %
    dots.params.ori = 0; % 0~360 : Orientation
    dots.params.frameRate = 60;
    dots.params.dotLifeOption = 'N'; % 'N' for normal dist, 'U' for uniform dist
    dots.params.dotLife = [50, 5]; % frames, [ Mean, SD ] of dots' life, SD is optional for 'N'
    nFrames = dots.params.frameRate * dots.params.duration;
end

%% (1) Create Field
dots.circField = [];
dots.circField.size= dots.params.fieldSize; 
dots.circField.oval = hb_cropOval( ones(dots.circField.size), dots.params.fieldShape, 0);
dots.circField.coords = [];
[dots.circField.coords(:,1),dots.circField.coords(:,2)]= ind2sub( size(dots.circField.oval), find(dots.circField.oval==1));
dots.circField.randperm = hb_randperm( size(dots.circField.coords,1) );

if dots.params.dimension == 3
    temp = [];
    temp.gaussField = hb_gaussianBlob( dots.params.speedField.size, dots.params.speedField.sd );
    temp.gaussField_mincut = temp.gaussField - min(min(temp.gaussField));
    temp.gaussField_maxcut = temp.gaussField_mincut / max(max(temp.gaussField_mincut));    
    dots.circField.gaussField =temp.gaussField_maxcut *dots.params.speedField.maxvalue;
    clear temp;    
    if dots.params.fieldShape == 2
        dots.circField.gaussField = repmat( dots.circField.gaussField( :, round(size(dots.circField.gaussField, 1)*.5)),...
            [1, size(dots.circField.gaussField, 2)]);
    end
    
    plotOption = 0;
    if plotOption
        meshz(dots.circField.gaussField); %colormap(gray);
    end
    
end
%% (2) Set movement parameters
dots.params.angles = (hb_Shuffle([ zeros([1,floor(dots.params.nDots*dots.params.coherence)]),...
    ones([1,dots.params.nDots-floor(dots.params.nDots*dots.params.coherence)]) ])...
    *360) .* rand([1,dots.params.nDots]);
dots.params.angles(find(dots.params.angles==0)) = dots.params.ori;

nFrames = floor(dots.params.frameRate * dots.params.duration);

%% (3) Initial prototype of dot configuration
dots.dotInfo = [];
for dotIdx = 1:dots.params.nDots
    dots.dotInfo(dotIdx).dotIdx = dotIdx;
    dots.dotInfo(dotIdx).angle = dots.params.angles(dotIdx);
    dots.dotInfo(dotIdx).dx = dots.params.speed*sin(dots.params.angles(dotIdx)*pi/180)/(dots.params.frameRate-1);
    dots.dotInfo(dotIdx).dy = dots.params.speed*cos(dots.params.angles(dotIdx)*pi/180)/(dots.params.frameRate-1);
    dots.dotInfo(dotIdx).XY = [dots.circField.coords(dots.circField.randperm(dotIdx),1)-(dots.circField.size(1)*.5),...
    dots.circField.coords(dots.circField.randperm(dotIdx),2)-(dots.circField.size(2)*.5)];
    switch(dots.params.dotLifeOption)
        case('U') % Uniform
            dots.dotInfo(dotIdx).remained_life = floor(rand()*dots.params.dotLife(1)*2)+1;
        case('N') % Normal
            dots.dotInfo(dotIdx).remained_life = floor(normrnd(dots.params.dotLife(1),dots.params.dotLife(2)))+1;
    end    
end

%% (4) Generate struct type dotFrames 
dots.dotFrames = [];
dots.dotFrames(1).dotInfo = dots.dotInfo;
dots.circField.randperm_counter = dots.params.nDots;
for frameIdx = 2:nFrames
    dots.dotFrames(frameIdx).dotInfo = dots.dotInfo;
    for dotIdx = 1:dots.params.nDots        
        refresh = 0;
        % Incremental change (life)
        life = dots.dotFrames(frameIdx-1).dotInfo(dotIdx).remained_life -1;
        if life > 0
            dots.dotFrames(frameIdx).dotInfo(dotIdx).remained_life = life;
        else
            % Get new dot XY
            newdotIdx = dots.circField.randperm(dots.circField.randperm_counter);
            dots.circField.randperm_counter=dots.circField.randperm_counter+1;
                XY =...
                    [dots.circField.coords(dots.circField.randperm(newdotIdx),1)-(dots.circField.size(1)*.5),...
                    dots.circField.coords(dots.circField.randperm(newdotIdx),2)-(dots.circField.size(2)*.5)];
            refresh = 1;
        end
            
        if ~refresh
            % Incremental change (pos)
            if dots.params.dimension == 2
                dx = dots.dotFrames(frameIdx).dotInfo(dotIdx).dx;
                dy = dots.dotFrames(frameIdx).dotInfo(dotIdx).dy;                
                XY = dots.dotFrames(frameIdx-1).dotInfo(dotIdx).XY+ [dx, dy];
            else
                
                XYidx = round(( dots.dotFrames(frameIdx-1).dotInfo(dotIdx).XY ) + dots.circField.size*.5);
                if XYidx(1) < 1; XYidx(1)=1; end; 
                if XYidx(1) > dots.circField.size(1); XYidx(1)=dots.circField.size(1); end; 
                if XYidx(2) < 1; XYidx(2)=1; end; 
                if XYidx(2) > dots.circField.size(2); XYidx(2)=dots.circField.size(2); end; 
                dx = dots.dotFrames(frameIdx).dotInfo(dotIdx).dx * ...
                    dots.circField.gaussField( XYidx(1),XYidx(2) );                
                dy = dots.dotFrames(frameIdx).dotInfo(dotIdx).dy * ...
                    dots.circField.gaussField(  XYidx(1),XYidx(2) );
                XY = dots.dotFrames(frameIdx-1).dotInfo(dotIdx).XY+ [dx, dy];
            end
            
            try % 트라이캐치 나중에 바꾸기
                inCircle = dots.circField.oval(round(XY(1)+(dots.circField.size(1)*.5)),round(XY(2)+(dots.circField.size(2)*.5)));
            catch
                inCircle = false;
            end
            
            if ~inCircle
                % Get new dot XY
                newdotIdx = dots.circField.randperm(dots.circField.randperm_counter);
                dots.circField.randperm_counter=dots.circField.randperm_counter+1;
                XY =...
                    [dots.circField.coords(dots.circField.randperm(newdotIdx),1)-(dots.circField.size(1)*.5),...
                    dots.circField.coords(dots.circField.randperm(newdotIdx),2)-(dots.circField.size(2)*.5)];
            end
        end
        dots.dotFrames(frameIdx).dotInfo(dotIdx).XY = XY;
            
    end
end

%% (5) Convert into 3-D matrix 
dots.RDKcoordinates = [];
for frameIdx = 1:nFrames
    XYs = [];
    for dotIdx = 1:dots.params.nDots
        XY = dots.dotFrames(frameIdx).dotInfo(dotIdx).XY;
        XYs = [XYs;XY];
    end
    dots.RDKcoordinates = cat(3, dots.RDKcoordinates, XYs);
end
RDKcoordinates = dots.RDKcoordinates;
return


%% Define CropOval
function resultImg = hb_cropOval(inputImg, spareFactor, bgrcolor)
if nargin < 3
    bgrcolor = 0;
end
hwidth = size(inputImg, 2) / 2.0;
hheight = size(inputImg, 1) / 2.0;

spareWidth = hwidth * spareFactor;
spareHeight = hheight * spareFactor;

[ww, hh] = meshgrid(1:hwidth, 1:hheight);

% simple ellipse equation gets us part three of your mask
mask_rightBottom = (((ww.^2)/spareWidth^2+(hh.^2)/spareHeight^2)<=1); 
mask_rightTop = flipud(mask_rightBottom);
mask_leftBottom = fliplr(mask_rightBottom);
mask_leftTop = flipud(mask_leftBottom);

mask_integrated = [mask_leftTop, mask_rightTop; ...
    mask_leftBottom, mask_rightBottom];

resultImg = inputImg;
[~,~,nDim] = size(resultImg);
if nDim == 1
    resultImg(mask_integrated(:,:)==0) = bgrcolor;
else
    multichannel_mask = repmat(mask_integrated,[1 1 nDim]);
    resultImg(multichannel_mask==0) = bgrcolor;
end

return
function shuffled_v = hb_Shuffle(v)
shuffled_v = v([hb_randperm(length(v))]);
return
function perm = hb_randperm(N)
%  USAGE -> perm= hb_randperm(N)
%  perm = hb_randperm(N) returns a vector containing a random permutation of the
%    integers 1:N.  For example, randperm(6) might be [2 4 5 6 1 3].
% 
[~, perm]=sort(rand([N,1]));
return
function f3 = hb_gaussianBlob(N,sigma)

if nargin < 1
    N = 30^2;
    sigma = sqrt(N) ^ 1.5;
end

% Generate basic gaussian blob
[x, y] = meshgrid(floor(-N/2):floor(N/2)-1, floor(-N/2):floor(N/2)-1);
f0 = exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
f1 = f0./sum(f0(:));

% Range re-scaling
f2 = (f1 - min(min(f1)));
f3 = 127.5  + (127.5 * (f2 / max(max(f2))));

return