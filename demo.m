% addpath hb_PTBsubfunctions/

dots = [];
dots.params.fieldSize = [200 200];
dots.params.fieldShape = 1;  % 1 : Circle, 2 : Square
dots.params.nDots = 100;      % number of dots
dots.params.coherence = .95 ;  % 0~1 : Proportion of coherent dot
dots.params.speed = 250;      % Pixel per second
dots.params.duration = 3;      % Second
dots.params.ori = 45;            % 0~360 : Orientation
dots.params.frameRate = 60;  % Hz
dots.params.dotLifeOption = 'N'; % 'N' for random sampling from normal dist., 'U' for uniform dist.
dots.params.dotLife = [25, 5]; % frames, [ Mean, SD ] of dots' life, SD is optional for 'N'

% 3D RDK parameters
dots.params.dimension = 3; % 2D- or 3D- RDK
if dots.params.dimension == 3
    % Gaussian kernel option
    dots.params.speedField.size = dots.params.fieldSize(1);
    dots.params.speedField.sd = dots.params.fieldSize(1)*.25;
    dots.params.speedField.maxvalue = 1;
end

% Get Dots
[RDKcoordinates1,dots] = hb_getRDK(dots);
dots.params.ori = dots.params.ori + 180; % Mirror-reversed orientation
[RDKcoordinates2,dots] = hb_getRDK(dots);
% title('Speed Factor (2D Gaussian)');

%% Demo 1. Visualization through plot
figure(1);
nFrames = dots.params.frameRate * dots.params.duration;
markerSize = 5;
for frameIdx = 1:nFrames
    plot( RDKcoordinates1(:, 1, frameIdx),RDKcoordinates1(:, 2, frameIdx), 'ro', 'MarkerSize', markerSize ); % Red-o
    hold on;
    plot( RDKcoordinates2(:, 1, frameIdx),RDKcoordinates2(:, 2, frameIdx), 'bo', 'MarkerSize', markerSize ); % Blue-o
    hold off;
    axis([-1 1 -1 1]*.75*dots.circField.size(1));
    xlabel('X coordinate'); ylabel('Y coordinate');
    axis xy;
    drawnow;
end

%% Demo 2. Visualization through PTB
nFrames = dots.params.frameRate * dots.params.duration;
Screen('Preference', 'SkipSyncTests',2);
PsychDefaultSetup(2);
scrsize = [0 0 dots.params.fieldSize(1) dots.params.fieldSize(2) ]*2;
[win,rect]= PsychImaging('OpenWindow', 0, .5, scrsize);
Screen('BlendFunction', win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
% dots.params.frameRate = round(1/Screen('GetFlipInterval', win));
dotsize= 3;  dottype = 0;
dotcolor1=[255 255 255]*1;
dotcolor2=[255 255 255]*0;
flipt= Screen('Flip', win);
for frameIdx = 1:nFrames
    Screen('DrawDots', win, RDKcoordinates1(:, :, frameIdx)',...
        dotsize, dotcolor1, rect(3:4)*.5, dottype);
    if dots.params.dimension == 3
        Screen('DrawDots', win, RDKcoordinates2(:, :, frameIdx)',...
            dotsize, dotcolor2, rect(3:4)*.5, dottype);
    end
    Screen('Flip', win, flipt+((1/dots.params.frameRate)*.7));
end
%
%Screen('CloseAll');



