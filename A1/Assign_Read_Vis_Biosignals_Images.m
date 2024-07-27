%% CMPT340 Assignment - Reading & Visualizing Bio-signals
% ============================ Goal =========================:
%   In this activity, you will load and visualize different data  
%   types acquired from different modalities. 
%
% ===================== Download Data   ==============================
%
% First, download and unzip the BioData from Canvas: 
%
%
% The data folder contains 5 files:
%
%       (1) emg_myopathy.txt (EMG)
%
%       (2) bMode.mat (B-mode Ultrasound)
%
%       (3) vol4d.mat (fMRI)
%
%       (4) XYStripes20x20Noise_9_L_3.mat (tensor field)
%       
%       (5) clathrin_ROI_cmpt340.csv (SMLM data)
%                 
%% ============================ Task 1 ===============================
%% (1) EMG Data 
% Data Description: EMG is short for electromyogram, which measures electrical activity in muscles
% Requirement: 
%      Load the data from emg_myopathy.txt and extract the two columns of values: 
%%             i) First column -- time
%%             ii) Second column -- the amplitude of an EMG signal. 

close all; clear; clc;
% Load data into MATLAB:
% Write your own code here:
data = importdata("emg_myopathy.txt");
time = data(:,1);
amp = data(:,2);

% Q1: What are the dimensions of the loaded variable (eg. 5 x 3 x 10)? 
% Write your own code here: 
dimensions = size(data);

disp(['Variable dimensions: ', num2str(dimensions)]); % display your answer

% How many EMG amplitude values does the signal contain?
% Write your own code here:
emg_values = size(amp);

disp(['Number of amplitudes: ', num2str(emg_values(:, 1))]); % display your answer

% Plot the amplitude vs time using the values in the two columns. 
% Properly label both axes. Add a suitable title to your figure.
figure(1)
plot(time, amp);
xlabel("Time");
ylabel("Amplitude");
title('Amplitude VS Time of an EMG signal');

% There is a period of time with a clearly reduced EMG amplitude. 
% Examine the previous plot, use the "Zoom In" icon on the toolbar 
% to close up on this region of the plot. 
% Visually (roughly) identify the approximate starting and ending 
% time for this period.
% Display these values :
% Write your own code here
start_time = 16.7638;
end_time = 16.9965;

disp(['Start : ', num2str(start_time), ' End: ', num2str(end_time)]); % update the display 

% What are the time and amplitude values of the 100th sample of the data? 
% Write your own code here
time_100 = time(100);
amp_100 = amp(100);

% update the display
disp(['Time_100th = ', num2str(time_100)]);
disp(['Amp_100th = ', num2str(amp_100)]);

% Plot the EMG signal from the 100th sample (in (1d)) up until the 1100th
% sample. Give a title to your figure. 
figure(2)
x_axis = time(100:1100);
y_axis = amp(100:1100);
plot(x_axis, y_axis)
xlabel("Time");
ylabel("Amplitude");
title('Amplitude VS Time for an EMG signal (only 100-1100 entries)');

% ============================ Task 2 ===============================
%% (2) B-mode Ultrasound sequence 
% Data Description: This data represent an B-mode ultrasound video sequence (B stands for Brightness)
% Requirement: Load the data in bMode.mat into MATLAB, explore the video frames and create a GIF
%  
close all; clear; clc;
% Load data into MATLAB:
% Write your own code here:
bmode_struct = load("bMode.mat");
bmode = bmode_struct.bMode;

% What are the dimensions of the loaded variable (eg. 5 x 3 x 10)? 
% You should output the variable dim_bMode as a vector of size 1x3.. 
% Write your own code here
bmode_dimensions = size(bmode);

disp(['bMode dimensions: ', num2str(bmode_dimensions)]); % display your answers

% How many frames (2D images) do we have in the video? 
% Write your own code here
frames = bmode_dimensions(:, 3);

disp(['Number of frames: ', num2str(frames)]); % display your answer

% What is the size (rows and columns) of each frame?
% Write your own code here
[rows, cols, ~] = size(bmode);

disp(['Nb_rows: ', num2str(rows)]);
disp(['Nb_cols: ', num2str(cols)]);

% Extract and display the 9-th frame of the sequence in a figure. 
% Apply the command that ensures that pixels are not stretched 
% (i.e. maintain aspect ratio of the image). Apply the command 
% that uses a grey colour-map to display the image.
% Write your own code here
figure(1)
ninth_frame = bmode(:,:,9);
imagesc(ninth_frame);
axis image;
colormap gray;
title('9th frame sequence');

% What is the intensity value at row 30 and column 40, at frame 20?
% Write your own code here 
intensity_value = bmode(30, 40, 20);
disp(['Intensity_(30,40,20): ', num2str(intensity_value)]); % update the display to show your answer

% Extract a cropped rectangular region of frame 15. 
% The cropped part should extend from row 10 to row 30, and from column 45 
% to column 60 of the image.  Display this cropped region 
% (also called region of interest or ROI).
% Write your own code here
% 
figure(2) % update the colormap and figure title
region = bmode(10:30, 45:60, 15);
imagesc(region);
axis image;
colormap gray;
title("ROI (Region of Interest)");

% Create an animated GIF that shows the ultrasound sequence
% (like a video clip). 
% Save this GIF image using the name US_seq.gif . 
% You can make use of the MATLAB code used during the lectures
% Write your own code here 
fId = figure(3);

for a=1:frames
   imagesc(bmode(:,:, a));
   axis image; 
   shading flat;
   pause(0.01);
   % camorbit(5,0,'data',[0 1 1])
   drawnow

   % This part records the GIF.
   frame = getframe(1);
   im = frame2im(frame);
   [A,map] = rgb2ind(im,256);
   filename = fullfile('US_seq.gif');
   if a == 1
       imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.25);
   else
       imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.25);
   end
   % end GIF stuff.

end
close(fId);

%% ============================ Task 3 ===============================
%% (3) Functional MRI
% Data Description: This is a 3D+time f-MRI data (functional magnetic resonance image). 
% fMRI captures the changes in oxygenation with time (4th dimension) at every (x,y,z) location of a volume.
% Requirement: Load the data in vol4d.mat into MATLAB, explore how oxygenation level changes with time

close all; clear; clc;
% Load data into MATLAB:
% Write your own code here:
vol4d_struct = load("vol4d.mat");
vol4d = vol4d_struct.vol4d;

% What are the dimensions of the loaded variable (eg. 5 x 3 x 10)? 
% Write your own code here:
vol4d_dimensions = size(vol4d);


disp(['Size fMRI: ', num2str(vol4d_dimensions)]); % display your answer

% Noting that the first 3 dimensions represent the spatial x, y, and z 
% dimensions, identify how many rows, columns, slices, and time steps 
% does this fMRI data have?
% Write your own code here:
[rows, cols, slices, time] = size(vol4d);

% update the display:
disp(['Size fMRI (rows): ', num2str(rows)]); % display # of rows in each slice
disp(['Size fMRI (cols): ', num2str(cols)]); % display # of columns in each slice
disp(['Size fMRI (slices): ', num2str(slices)]); % display # of slices
disp(['Size fMRI (time-steps): ', num2str(time)]); % display time

% Plot a curve that shows how the oxygenation level changes 
% with time at voxel location (x=32, y=34 , z=10).
% Define axis and figure title. 
figure(1);
vox = vol4d(32, 34, 10, :);
vox = squeeze(vox);
plot(1:size(vox,1),vox,'-o');
title('The oxygenation levels changing over time');
xlabel('time');
ylabel('oxegenation level');
% 
% % Extract and display the slice in the X-Y plane, at slice x
% % location z=15, of the 10th time sample.
% 
figure(2);
xyslice = vol4d(:, :, :, 10);
slice(xyslice, [], [], 15)
title("X-Y plane slice at z=15, t=10")
xlabel('x axis');
ylabel('y axis');
zlabel('z axis');
% 
% 
% % Extract and visualize the slice in the X-Z plane, 
% % at slice location y=10, of the 15th time sample.
% 
figure(3)
xzslice = vol4d(:, :, :, 15);
slice(xzslice,[], 10, [])
title("X-Z plane slice at y=10, t=15")
xlabel('x axis');
ylabel('y axis');
zlabel('z axis');

% Create an animated GIF that shows how the slice in the previous question
% changes with time (over the whole time steps). 
% Save this GIF image using the name fMRI.gif. 
% You can make use of the MATLAB code used during the lectures. 

% Write your own code here:
fId = figure(4);
for a=1:time
    vol = vol4d(:, :, :, a);
    slice(vol, [], 10, [])
    axis tight;
    shading flat;
    colorbar

    title(sprintf('3D volume:%i',a));
    %snapnow;
    frame = getframe(1);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);
    filename = fullfile('fMRI.gif');
    if a == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
    end

    pause(0.5)
end
close(fId);

%% ============================ Task 4 ===============================
%% (5) Diffusion Tensor MRI 
% Data Description: This is a single 2D slice of a Diffusion Tensor MRI data, 
% which is a 2D image with a 3x3 matrix at every pixel.
% Requirement: Load the data XYStripes20x20Noise_9_L_3.mat into MATLAB, calculate the eigenvalues.

close all; clear; clc;
% Load data into MATLAB:
% Write your own code here:
xystripe_struct = load("XYStripes20x20Noise_9_L_3.mat");
TensorIM = xystripe_struct.TensorIM;

% What is the dimensions of the loaded variable (eg. 5 x 3 x 10)?
% Write your own code here:
TensorIM_dimensions = size(TensorIM);

disp(['Size fMRI: ', num2str(TensorIM_dimensions)]); %display your answer 

% How many rows and columns of pixels does this 2D image have?
% Write your own code here: 
[rows, cols, ~, ~] = size(TensorIM);

disp(['Nb rows:', num2str(rows) ', Nb cols:', num2str(cols)]); % update to display the correct values

% Extract and display (as an array of numbers) the 3x3 
% matrix at pixel location: (row 10, column 15). Use variable name *arr3x3*

% your code here, update the display to print your variables
arr3x3 = squeeze(TensorIM(10, 15, :, :));

disp('Array at x=10, y=15 : ')
disp(arr3x3);

% Use the "eig" function in MATLAB to calculate and display the 3 eigen
% vectors and the corresponding 3 eigen values of the 3x3 matrix arr3x3;
% Use the variable name *arr3x3_eig*.

% your code here, add a display at the end of your code to print out your
% answers
[arr3x3_eig, arr3x3_val] = eig(arr3x3);

disp('Eigenvectors : '); 
disp(arr3x3_eig);
disp('Eigenvalues : '); 
disp(arr3x3_val);

% Create and display a 2D image *DT_2D* with a scalar number at every pixel. 
% The scalar number should be the entry (row=1, column=2) of the 
% 3x3 matrix at every pixel.

% Write your own code here:
figure(1)
image = zeros(size(TensorIM,1));
for i = 1:size(TensorIM, 1)
    for j = 1:size(TensorIM, 2)
        image(i, j) = TensorIM(i, j, 1, 2);
    end
end
% disp(image);
imagesc(image);
axis image;
colorbar;


% Create and display a 2D image with a scalar number at every pixel. 
% The scalar number should be the largest eigenvalue of the 3x3 matrix at 
% each pixel. 
% The challenge is to do this WITHOUT using loops (without for or while).

% your code here
% figure(2)
% image2 = zeros(size(TensorIM, 1));
% for i = 1:size(TensorIM, 1)
%     for j = 1:size(TensorIM, 2)
%         mat_3x3 = squeeze(TensorIM(i, j, :, :));
%         [~, val] = eig(mat_3x3);
%         max_eig_val = max(max(val));
%         disp(max_eig_val)
%         image2(i, j) = max_eig_val;
%     end
% end
% imagesc(image2);
% colorbar;
figure(2)
submatrix_dimensions = [3, 3];
matrix_dimensions = [20, 20];
% at every pixel, extract the 3x3 matrix
submatrix = mat2cell(TensorIM, ones(1, matrix_dimensions(1)), ones(1, matrix_dimensions(2)), submatrix_dimensions(1), submatrix_dimensions(2));
% find the eigenvalue
eig_val = cellfun(@(x) eig(reshape(x, 3, 3)), submatrix, 'UniformOutput', false);
% extract max eigenvalue at each pixel
max_eig_val = cellfun(@(x) max(x), eig_val);
% show the image
image2 = max_eig_val;
imagesc(image2);
axis image;
colorbar;


%% ============================ Task 5 ===============================
%% (6) SMLM -- Single-Molecule Localization Microscopy
% Data Description: 3D point cloud data illustrating clathrin protein locations. 
% The 3 channels in the 2nd column represents 3D coordinates.

close all; clear; clc;
% Load data into MATLAB:
% Write your own code here:
dat = importdata('clathrin_ROI_cmpt340.csv');
rawData = dat.data;

% What is the dimensions of the loaded variable (eg. 5 x 3 x 10)? 
% your code here
rawData_dimensions = size(rawData);

disp(['Size SMLM: ', num2str(rawData_dimensions)]); %display your answer 

% Select 1001 points from the imported data. You can clip from the middle, e.g., use format such as rawData(ind-xxx:ind:xxx)
% Write your code here
points3d = rawData(round(end/2)-500:round(end/2)+500, :);

%% Plot the raw data
figure(1) % you can use proper 'position' property to make the figure look nice
% Write your code here
% Use colormap jet, black background
% Label X, Y, Z dimensions properly, display Unit = 'nm'
% rotate camera to visualize your data
scatter3(points3d(:, 1), points3d(:, 2), points3d(:, 3), 20, 'b', 'filled')
set(gca,'Color','k');
axis equal;
xlabel('X nm');
ylabel('Y nm');
zlabel('Z nm');
title('Scatter plot');
grid on;
colormap(jet);
rotate3d;

%% Construct graph (this method only works for small data) 
% this is the naive way to construct the graph using the distance matrix (for small data)
PT = 800; % set proximity threshold -- proximity threshold means how you decide to clip the data according to the distances among neighboring points
% Calculate pairwise distance among all the points
distData = pdist(points3d);
% Use adjacency matrix representation
distMatrix = squareform(distData);
% Threshold the adjacency distance matrix
distMatrix_thresh = distMatrix.*(distMatrix <= PT);
% Create a graph from the distance matrix, save it as variable g
g = graph(distMatrix_thresh, array2table(points3d));

%% What if you selected more than 1001 points from the data 
% When the size of the data is large, creating a graph from 
% the (thresholded) distance matrix can produce out-of-memory errors. 
% An alternative is to call the graph by supplying the source, the target, 
% and the weights of the edges in the graph, i.e.,
%      g = graph(source, target, weight);
% 
% Write the code (which involves a loop over the number of nodes) 
% that can create source, target and weight from d=pdist(rawDat);
% Write your code here
dim_points3d = size(points3d, 1);
source = [];
target = [];
weight = [];
for i = 1:dim_points3d
    for j = i+1:dim_points3d
        buffer = j-i;
        index_corrector = i-1;
        threshold = (dim_points3d-i/2)*index_corrector + buffer;
        if distData(threshold) <= PT
            source = [source, i];
            target = [target, j];
            weight = [weight, distData(threshold)];
        end
    end
end
g = graph(source, target, weight);

% plot the constructed graph preserving the spatial info
figure(2)
p= plot(g, 'XData', points3d(:, 1),'YData', points3d(:, 2),'ZData', points3d(:, 3),'NodeLabel',{}); % Fill in the code here
axis equal;
% set marker to 'o', use a marker size=2, red color and white edge color
% set X, Y, Z axis labels
% Write your code here
set(gca,'Color','k')
p.Marker = 'o';
p.MarkerSize = 2;
p.NodeColor = 'r';
p.EdgeColor = 'w';
xlabel('X nm');
ylabel('Y nm');
zlabel('Z nm');
grid on;
view(3);

% %% find the node degree, save it into variable deg
% % Write your code here
deg = degree(g);
% 
% %% plot the points colored by their degree
figure(3)

% Write your code here
% use colormap jet, black background, set X, Y, Z axis labels using unit='nm'
% Write your code here
scatter3(points3d(:, 1), points3d(:, 2), points3d(:, 3), 20, deg, 'filled')
set(gca,'Color','k')
axis equal;
colormap(jet);
xlabel('X nm');
ylabel('Y nm');
zlabel('Z nm');
grid on;
colorbar;
rotate3d;

% %%%%%%%%%%%%%%
figure(4)
% % For the original graph (created using the adjacency matrix) change graph visualization layout to 'force' 
% % Write your code here
plot = plot(g, 'Layout', 'Force');