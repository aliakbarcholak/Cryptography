clc; clear; close all;

% Set folder path for images
folder_path = 'H:\DataBase\TIFF'; % Change to the actual folder path
image_files = dir(fullfile(folder_path, '*.TIFF')); % Get all JPG images

% Check if images exist
if isempty(image_files)
    error('No images found in the folder. Ensure the correct folder path and images are available.');
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
    img = imresize(img, [256, 256]); % Resize to fixed size
    
    % Convert to grayscale if needed
    img = im2gray(img); % Handles both RGB and grayscale images
    img = im2double(img); % Convert to double precision for processing

    % Encrypt the image using SCAN (Lossless Image Compression and Encryption)
    encrypted_img = scan_encryption(img);

    % Decrypt the image using SCAN inverse algorithm
    decrypted_img = scan_decryption(encrypted_img);

    % Compute MSE (Original vs Encrypted)
    MSE_encrypted(i) = immse(img, encrypted_img);

    % Compute MSE (Original vs Decrypted)
    MSE_decrypted(i) = immse(img, decrypted_img);

    % Display results
    fprintf('Image %d: %s\n', i, image_files(i).name);
    fprintf('  MSE (Original vs Encrypted): %.4f\n', MSE_encrypted(i));
    fprintf('  MSE (Original vs Decrypted): %.4f\n\n', MSE_decrypted(i));

    % Display images
    figure;
    subplot(1, 3, 1); imshow(img, []); title('Original Image');
    subplot(1, 3, 2); imshow(encrypted_img, []); title('Encrypted Image');
    subplot(1, 3, 3); imshow(decrypted_img, []); title('Decrypted Image');
end

% Display overall results
disp('MSE between Original and Encrypted Images:');
disp(MSE_encrypted);
disp('MSE between Original and Decrypted Images:');
disp(MSE_decrypted);

%% ---- Helper Functions ----

% SCAN-Lossless Image Compression and Encryption Algorithm
function encrypted_img = scan_encryption(img)
    % Simulate SCAN-based encryption
    % This is a simplified version for demonstration purposes
    % Step 1: Perform zigzag scanning (simplified)
    [rows, cols] = size(img);
    zigzag = [];
    for i = 1:rows
        for j = 1:cols
            zigzag = [zigzag img(i,j)];
        end
    end
    
    % Step 2: Encrypt by reversing the scanned array
    encrypted_img = flip(zigzag);
    encrypted_img = reshape(encrypted_img, [rows, cols]);
end

% SCAN Decryption (inverse of the SCAN encryption)
function decrypted_img = scan_decryption(encrypted_img)
    % Step 1: Reverse the encrypted image array (undo flip)
    [rows, cols] = size(encrypted_img);
    zigzag = reshape(encrypted_img, 1, rows*cols);
    decrypted_zigzag = flip(zigzag);
    
    % Step 2: Reconstruct the image from the decrypted zigzag pattern
    decrypted_img = reshape(decrypted_zigzag, [rows, cols]);
end
