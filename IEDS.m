clc; clear; close all;

% Set folder path for images
folder_path = 'C:\Users\admin\Desktop\DATASET'; % Change to the actual folder path
image_files = dir(fullfile(folder_path, '*.PNG')); % Get all JPG images

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
    img = imresize(img, [1024, 1024]); % Resize to fixed size
    
    % Convert to grayscale if needed
    img = im2gray(img); % Handles both RGB and grayscale images
    img = im2double(img); % Convert to double precision for processing

    % Encrypt the image using Image Encryption using Digital Signatures (IEDS)
    [encrypted_img, signature] = ieds_encryption(img);

    % Decrypt the image
    decrypted_img = ieds_decryption(encrypted_img, signature);

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

% Image Encryption using Digital Signatures (IEDS) Algorithm
function [encrypted_img, signature] = ieds_encryption(img)
    % Generate a digital signature by hashing the image
    signature = hash_image(img);
    
    % Simulate encryption: XOR the image with the signature (for demonstration)
    encrypted_img = xor(img > 0.5, signature > 0.5); % XOR encryption
    
    % Convert to double for consistency
    encrypted_img = double(encrypted_img);
end

% Image Decryption using Digital Signatures (inverse of IEDS)
function decrypted_img = ieds_decryption(encrypted_img, signature)
    % Decrypt the image: XOR the encrypted image with the signature
    decrypted_img = xor(encrypted_img > 0.5, signature > 0.5); % XOR decryption
    
    % Convert to double for consistency
    decrypted_img = double(decrypted_img);
end

% Hashing function to generate a digital signature for an image
function signature = hash_image(img)
    % Create a hash based on the image content (example using a simple sum)
    % You can replace this with more sophisticated cryptographic hashing methods
    img_sum = sum(img(:));  % Sum of all pixel values (simple approach)
    
    % Convert the sum to a binary signature (for simplicity)
    signature = logical(mod(img_sum, 2)); % Simple binary signature (0 or 1)
    signature = repmat(signature, size(img)); % Expand the signature to match the image size
end
