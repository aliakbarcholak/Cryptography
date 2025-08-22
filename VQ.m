clc; clear; close all;

% Set folder path for images (CHANGE THIS TO YOUR FOLDER)
folder_path = 'H:\small-DataBase\JPEG'; % Modify with your actual folder path
image_files = dir(fullfile(folder_path, '*.JPEG')); % Get all JPG images

% Check if images exist
if isempty(image_files)
    error('No images found in the folder. Ensure the correct folder path and images are available.');
end

% Initialize MSE storage
num_images = length(image_files);
MSE_encrypted = zeros(1, num_images);
MSE_decrypted = zeros(1, num_images);

% Set VQ parameters
block_size = 4; % Size of blocks for vector quantization
num_clusters = 256; % Number of clusters (codebook size)

% Loop through all images
for i = 1:num_images
    % Read image
    img_path = fullfile(folder_path, image_files(i).name);
    img = imread(img_path);
    img = imresize(img, [128, 128]); % Resize to fixed size
    
    % Convert to grayscale if needed
    img = im2gray(img); % Handles both RGB and grayscale images
    img = im2double(img); % Convert to double precision for processing

    % Apply VQ encryption
    [encrypted_img, codebook] = vq_encrypt(img, block_size, num_clusters);

    % Decrypt the image using the same codebook
    decrypted_img = vq_decrypt(encrypted_img, codebook, block_size, size(img));

    % Resize encrypted image to original dimensions before MSE calculation
    encrypted_resized = imresize(encrypted_img, size(img));  

    % Compute MSE (Original vs Encrypted)
    MSE_encrypted(i) = immse(img, encrypted_resized);

    % Compute MSE (Original vs Decrypted)
    MSE_decrypted(i) = immse(img, decrypted_img);

    % Display results
    fprintf('Image %d: %s\n', i, image_files(i).name);
    fprintf('  MSE (Original vs Encrypted): %.4f\n', MSE_encrypted(i));
    fprintf('  MSE (Original vs Decrypted): %.4f\n\n', MSE_decrypted(i));

    % Display images
    figure;
    subplot(1, 3, 1); imshow(img, []); title('Original Image');
    subplot(1, 3, 2); imshow(encrypted_resized, []); title('Encrypted Image (Resized)');
    subplot(1, 3, 3); imshow(decrypted_img, []); title('Decrypted Image');
end

% Display overall results
disp('MSE between Original and Encrypted Images:');
disp(MSE_encrypted);
disp('MSE between Original and Decrypted Images:');
disp(MSE_decrypted);

%% ---- Helper Functions ----

% VQ Encryption Function
function [encrypted_img, codebook] = vq_encrypt(img, block_size, num_clusters)
    [rows, cols] = size(img);

    % Reshape image into blocks
    blocks = im2col(img, [block_size, block_size], 'distinct');

    % Apply k-means clustering to form the codebook
    [idx, codebook] = kmeans(blocks', num_clusters);

    % Replace blocks with cluster indices for encryption
    encrypted_img = reshape(idx, [rows / block_size, cols / block_size]);
end

% VQ Decryption Function
function decrypted_img = vq_decrypt(encrypted_img, codebook, block_size, orig_size)
    [rows, cols] = size(encrypted_img);
    
    % Map indices back to codebook values
    blocks = codebook(encrypted_img(:), :);
    
    % Reshape back to image and resize to original size
    decrypted_img = col2im(blocks', [block_size, block_size], orig_size, 'distinct');
end
