clc; clear; close all;

% Set folder path for images
folder_path = 'F:\DataBase\JPEG'; % Change to the actual folder path
image_files = dir(fullfile(folder_path, '*.jpg')); % Get all JPG images

% Check if images exist
if isempty(image_files)
    error('No images found in the folder. Ensure the correct folder path and images are available.');
end

% Initialize MSE storage
num_images = length(image_files);
MSE_encrypted = zeros(1, num_images);
MSE_decrypted = zeros(1, num_images);

% Chaotic map parameters (Logistic Map)
r = 3.99; % Chaotic parameter (0 < r ≤ 4)
x0 = 0.5; % Initial condition

% Loop through all images
for i = 1:num_images
    % Read image
    img_path = fullfile(folder_path, image_files(i).name);
    img = imread(img_path);
    img = imresize(img, [256, 256]); % Resize to fixed size
    
    % Convert to grayscale if needed
    img = im2gray(img); % Handles both RGB and grayscale images
    img = im2double(img); % Convert to double precision for processing

    % Generate chaotic sequence
    chaotic_seq = logistic_map(r, x0, numel(img));

    % Encrypt image using chaotic XOR operation
    encrypted_img = xor(img > 0.5, chaotic_seq > 0.5);

    % Decrypt image (since XOR is symmetric)
    decrypted_img = xor(encrypted_img, chaotic_seq > 0.5);

    % Convert logical encrypted and decrypted images back to double for MSE calculation
    encrypted_img = double(encrypted_img);
    decrypted_img = double(decrypted_img);

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

% Logistic Map Function for Chaotic Sequence
function chaotic_seq = logistic_map(r, x0, N)
    chaotic_seq = zeros(1, N);
    x = x0;
    for k = 1:N
        x = r * x * (1 - x);
        chaotic_seq(k) = x;
    end
    chaotic_seq = reshape(chaotic_seq, sqrt(N), sqrt(N));
end
