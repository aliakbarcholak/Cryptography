clc; clear; close all;

% Set the folder path containing images
folder_path = 'H:\DataBase\BMP'; % Change this to your actual folder path
image_files = dir(fullfile(folder_path, '*.BMP')); % Get list of JPG images

% Check if images exist
if isempty(image_files)
    error('No images found in the folder. Make sure the folder path is correct and contains images.');
end

% Initialize MSE storage
num_images = length(image_files);
MSE_encrypted = zeros(1, num_images);
MSE_decrypted = zeros(1, num_images);

% Loop through all images
for i = 1:num_images
    % Read image
    img_path = fullfile(folder_path, image_files(i).name);
    img = imread(img_path);
    img = imresize(img, [128, 128]); % Resize to a fixed size

    % Convert to grayscale if needed
    img = im2gray(img); % Works for both RGB and grayscale images
    img = im2double(img); % Convert to double for Fourier processing

    % Apply DRPE Encryption
    [encrypted_img, phase_key1, phase_key2] = drpe_encrypt(img);

    % Decrypt image by applying inverse DRPE
    decrypted_img = drpe_decrypt(encrypted_img, phase_key1, phase_key2);

    % Compute MSE (Original vs Encrypted)
    MSE_encrypted(i) = immse(img, abs(encrypted_img));

    % Compute MSE (Original vs Decrypted)
    MSE_decrypted(i) = immse(img, abs(decrypted_img));

    % Display results
    fprintf('Image %d: %s\n', i, image_files(i).name);
    fprintf('  MSE (Original vs Encrypted): %.4f\n', MSE_encrypted(i));
    fprintf('  MSE (Original vs Decrypted): %.4f\n\n', MSE_decrypted(i));

    % Display images
    figure;
    subplot(1, 3, 1); imshow(img, []); title('Original Image');
    subplot(1, 3, 2); imshow(log(1 + abs(encrypted_img)), []); title('Encrypted Image (Fourier Domain)');
    subplot(1, 3, 3); imshow(decrypted_img, []); title('Decrypted Image');
end

% Display overall results
disp('MSE between Original and Encrypted Images:');
disp(MSE_encrypted);
disp('MSE between Original and Decrypted Images:');
disp(MSE_decrypted);

%% ---- Helper Functions ----

% DRPE Encryption Function
function [encrypted_img, phase_key1, phase_key2] = drpe_encrypt(img)
    [rows, cols] = size(img);

    % Generate two random phase masks
    phase_key1 = exp(1i * 2 * pi * rand(rows, cols));
    phase_key2 = exp(1i * 2 * pi * rand(rows, cols));

    % Apply DRPE: Fourier Transform -> Multiply by Phase Key 1 -> Fourier Transform -> Multiply by Phase Key 2
    fft_img = fft2(img);
    encrypted_img = fft2(fft_img .* phase_key1) .* phase_key2;
end

% DRPE Decryption Function
function decrypted_img = drpe_decrypt(encrypted_img, phase_key1, phase_key2)
    % Apply inverse DRPE: Inverse Fourier Transform -> Remove Phase Key 2 -> Inverse Fourier Transform -> Remove Phase Key 1
    ifft_img = ifft2(encrypted_img ./ phase_key2);
    decrypted_img = abs(ifft2(ifft_img ./ phase_key1)); % Take absolute value to recover real part
end
