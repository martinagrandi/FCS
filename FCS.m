%% EX 1

% Here you give as input the image with a lot of dots. You set the pixel
% size and the number of measurement yoou want to do (iterations, here 10). The code
% automatically will ask you 10 times to select two points. Then it saves
% the resolution (half the distance of the 2 points) in an array and
% performs some data analysis. What you have at the end is the measure of
% the mean value and the error associated to it. =)

file_tiff1 = 'ex1_px100nm.tif'; % write the name of the image you want to use
px_size = 100; %nm
iterations = 10;
res = resolution(file_tiff1, px_size, iterations);
M_res = mean(res);
st_dev_res = sqrt(sum((res-M_res).^2)/(length(res)-1));
err_mean_res = st_dev_res/ sqrt(length(res));
fprintf('EXPERIMENT 1: The resolution is (%s pm %d)nm .\n',M_res,err_mean_res);


%% EX 2 - DATA EXTRACTION

% Here it is a bit less smooth. You should run several times this section
% of the code, saving by hand the value of c that will appear and the
% lower bound of the confidence bounds. When ou have done that, go to next
% section =)

file_tiff = 'ex2_im4_50nm.tif'; % write the name of the image you want to use

%Cut the image
cut_im = cutImageTIFF(file_tiff);
[~, name, esten] = fileparts(file_tiff);
name_cut_file = fullfile([name '_tagliata' esten]); %builds the name of the new file

%Finds the value of the pixels in the line
values = readPixelTIFF(name_cut_file);

% Gaussian fitting
nr = length(values); %nr of the pixel considered
x = 1:1:nr;
f = fit(x',values','gauss1') % fit to a Gaussian curve
figure;
plot(f,x,values)
hold off

%% EXPERIMENT 2 - DATA ANALYSIS

% When you have saved the values of c, we convert it into real distances
% (*px_size2) and we perform some data analysis.


px_size2 = 50;
c2 = [8.8486 10.3403 13.5361 12.5305 14.8732]*px_size2; %nm
dc2 = [0.7558 0.90355 0.64105 0.63445 0.9225]*px_size2; %nm

% This is the mean value of the FWHM with the associated error
M2 = mean(c2);
st_dev2 = sqrt(sum((dc2.^2 .* (c2 - M2).^2)) / sum(dc2.^2));
err_mean2 = st_dev2/ sqrt(length(c2));
fprintf('EXPERIMENT 2a: The mean is (%s pm %d)nm \n', M2,err_mean2);


% We compute the FWHM_res with the formula FWHM^2_meas - FWHM^2_sphere = FWHM^2_res
FWHM_res = sqrt(c2.^2 - 100^2);
delta = c2./FWHM_res.*dc2;
M3 = mean(FWHM_res);
st_dev3 = sqrt(sum((delta.^2 .* (FWHM_res - M3).^2)) / sum(delta.^2));
err_mean3 = st_dev3/ sqrt(length(FWHM_res));
fprintf('EXPERIMENT 2b: The FWHMW_resolution is (%s pm %d)nm \n ', M3,err_mean3);




function cut_im = cutImageTIFF(file_tiff)
    % Read the TIFF file
    immage = imread(file_tiff);

    % Visualize the image
    imshow(immage);
    title('Select the area you want to cut');

    % Allow the user to select a rectangular shape
    rect = round(getrect);

    % Save the coordiantes of the selected rectangle
    x1 = rect(1);
    y1 = rect(2);
    width = rect(3);
    height = rect(4);

    % Cut the image
    cut_im = immage(y1:y1+height-1, x1:x1+width-1, :);

    % Save the cut image in the current folder
    [~, name, esten] = fileparts(file_tiff);
    name_cut_file = fullfile(pwd, [name '_tagliata' esten]);
    imwrite(cut_im, name_cut_file);
    
    disp(['Image saved as "', name_cut_file, '".']);
end

function pixel_values = readPixelTIFF(file_tiff)
    % Read TIFF file
    image = imread(file_tiff);
    
    % Visualize the image and select 2 points inside it
    imshow(image);
    title('Select two points inside the image');
    [x, y] = ginput(2);
    coord1 = [x(1), y(1)];
    coord2 = [x(2), y(2)];

    % Compute the distance between the points
    distance = sqrt((coord2(1) - coord1(1))^2 + (coord2(2) - coord1(2))^2);

    % Compute the nr of pixel on the line
    num_pixel = ceil(distance) + 1;

    % Compute the coordinates of pixels along the line
    x_coords = round(linspace(coord1(1), coord2(1), num_pixel));
    y_coords = round(linspace(coord1(2), coord2(2), num_pixel));

    % Ottieni i valori dei pixel lungo il segmento
    pixel_values = image(sub2ind(size(image), y_coords, x_coords));
end

function res = resolution(file_tiff,px_size,iteration)
    % Read TIFF file
    image = imread(file_tiff);
    res = zeros(1,iteration);
    for i = 1:iteration
        % Visualize the image and select 2 points inside it
        imshow(image);
        title('Select two points inside the image');
        [x, y] = ginput(2);
        coord1 = [x(1), y(1)];
        coord2 = [x(2), y(2)];

        % Compute the resolution between the points
        res(i) = 0.5*px_size*sqrt((coord2(1) - coord1(1))^2 + (coord2(2) - coord1(2))^2);
    end
end
