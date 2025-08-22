clc; clear; close all;

% Set folder path containing images (CHANGE THIS TO YOUR FOLDER)
folder_path = 'H:\DataBase\TIFF'; % Change this to your actual folder path
image_files = dir(fullfile(folder_path, '*.TIFF')); % Get list of JPG images

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
    img = im2uint8(img); % Convert to uint8 format

    % Apply Visual Cryptography
    [share1, share2] = visual_cryptography(img);

    % Decrypt image by combining the shares
    decrypted_img = bitxor(share1, share2);

    % Ensure all images are of the same size and data type
    share1 = imresize(share1, size(img(:,:,1))); % Resize share1
    share2 = imresize(share2, size(img(:,:,1))); % Resize share2
    decrypted_img = imresize(decrypted_img, size(img(:,:,1)));

    % Compute MSE (Original vs Encrypted)
    MSE_encrypted(i) = immse(double(img), double(share1));

    % Compute MSE (Original vs Decrypted)
    MSE_decrypted(i) = immse(double(img), double(decrypted_img));

    % Display results
    fprintf('Image %d: %s\n', i, image_files(i).name);
    fprintf('  MSE (Original vs Encrypted): %.4f\n', MSE_encrypted(i));
    fprintf('  MSE (Original vs Decrypted): %.4f\n\n', MSE_decrypted(i));

    % Display images
    figure;
    subplot(1, 4, 1); imshow(img); title('Original Image');
    subplot(1, 4, 2); imshow(share1); title('Share 1 (Encrypted)');
    subplot(1, 4, 3); imshow(share2); title('Share 2 (Encrypted)');
    subplot(1, 4, 4); imshow(decrypted_img); title('Decrypted Image');
end

% Display overall results
disp('MSE between Original and Encrypted Images:');
disp(MSE_encrypted);
disp('MSE between Original and Decrypted Images:');
disp(MSE_decrypted);

%% ---- Helper Functions ----

% Function for Visual Cryptography (VC) for Color Images
function [share1, share2] = visual_cryptography(img)
    [rows, cols, channels] = size(img);
    share1 = uint8(randi([0, 255], rows, cols, channels)); % Generate random share
    share2 = bitxor(img, share1); % XOR with original image to get the second share
end
