clc; clear; close all;

% Set the folder path containing images
folder_path = 'H:\DataBase\BMP'; % !!! IMPORTANT: Using the original path from this file
image_files = dir(fullfile(folder_path, '*.BMP')); % Get list of BMP images

% Check if images exist
if isempty(image_files)
    error('No images found in the folder. Make sure the folder path is correct and contains BMP images.');
end

% Loop through all images (set to 1 for a single demonstration run)
% For a full run on all images, change '1' to 'length(image_files)'
for i = 1:1
    % --- Image Preparation ---
    img_path = fullfile(folder_path, image_files(i).name);
    img_orig = imread(img_path);

    % Resize for faster processing and consistency. Original was 128x128.
    img_orig = imresize(img_orig, [256, 256]);

    % Convert to grayscale and then to double format [0, 1] for processing
    img_gray = im2gray(img_orig);
    img_double = im2double(img_gray);

    fprintf('============================================================\n');
    fprintf('Processing Image: %s\n', image_files(i).name);
    fprintf('============================================================\n\n');

    % --- Encryption and Decryption ---
    [encrypted_img, phase_key1, phase_key2] = drpe_encrypt(img_double);
    decrypted_img = drpe_decrypt(encrypted_img, phase_key1, phase_key2);

    % --- Metric Calculation ---
    % For standardized metrics like entropy, PSNR, and SSIM, we convert
    % the images to the standard uint8 format [0, 255].
    original_uint8 = im2uint8(img_double);
    encrypted_uint8 = im2uint8(abs(encrypted_img));
    decrypted_uint8 = im2uint8(decrypted_img);

    % 1. Entropy: Measures the randomness of the encrypted image.
    entropy_val = entropy(encrypted_uint8);

    % 2. PSNR (Peak Signal-to-Noise Ratio): Measures image quality.
    psnr_encrypted = psnr(encrypted_uint8, original_uint8);
    psnr_decrypted = psnr(decrypted_uint8, original_uint8);

    % 3. SSIM (Structural Similarity Index): Measures structural similarity.
    ssim_encrypted = ssim(encrypted_uint8, original_uint8);
    ssim_decrypted = ssim(decrypted_uint8, original_uint8);

    % 4. & 5. NPCR and UACI (Differential Attack Analysis)
    % Create a slightly modified original image by changing one pixel.
    img_double_mod = img_double;
    % Flip the least significant bit of one pixel to create a minimal change.
    img_double_mod(1,1) = img_double_mod(1,1) + 1/255;
    if img_double_mod(1,1) > 1.0
        img_double_mod(1,1) = img_double_mod(1,1) - 2/255; % Handle wrap-around
    end

    % Encrypt the modified image using the SAME keys generated previously.
    fft_img_mod = fft2(img_double_mod);
    encrypted_img_mod = fft2(fft_img_mod .* phase_key1) .* phase_key2;
    encrypted_uint8_mod = im2uint8(abs(encrypted_img_mod));

    % Calculate NPCR and UACI by comparing the two encrypted images.
    [npcr_val, uaci_val] = calculate_npcr_uaci(encrypted_uint8, encrypted_uint8_mod);

    % --- Display All Results in Command Window ---

    fprintf('--- 1. Encryption Performance ---\n');
    fprintf('Entropy of Encrypted Image             : %.4f (Ideal > 7.8)\n', entropy_val);
    fprintf('PSNR (Original vs. Encrypted)        : %.4f dB (Ideal: Low, e.g., < 10)\n', psnr_encrypted);
    fprintf('SSIM (Original vs. Encrypted)        : %.4f (Ideal: Low, e.g., < 0.1)\n', ssim_encrypted);
    fprintf('\n');

    fprintf('--- 2. Differential Attack Resistance ---\n');
    fprintf('NPCR (Number of Pixels Change Rate)    : %.4f %% (Ideal > 99.6%%)\n', npcr_val);
    fprintf('UACI (Unified Average Changing Intensity): %.4f %% (Ideal ~ 33.46%%)\n', uaci_val);
    fprintf('\n');

    fprintf('--- 3. Decryption Quality ---\n');
    fprintf('PSNR (Original vs. Decrypted)        : %.4f dB (Ideal: High, e.g., > 30)\n', psnr_decrypted);
    fprintf('SSIM (Original vs. Decrypted)        : %.4f (Ideal: High, close to 1.0)\n', ssim_decrypted);
    fprintf('\n');

    % --- Display Images ---
    figure('Name', ['Encryption Analysis for ' image_files(i).name]);
    subplot(1, 3, 1); imshow(img_double, []); title('Original Image');
    subplot(1, 3, 2); imshow(encrypted_uint8, []); title('Encrypted Image');
    subplot(1, 3, 3); imshow(decrypted_img, []); title('Decrypted Image');
end

% --- Explanation of Metrics ---
fprintf('============================================================\n');
fprintf('How to Interpret the Security Metrics:\n');
fprintf('============================================================\n\n');
fprintf(['1. Entropy: Measures the randomness of the encrypted image. For an 8-bit image,\n' ...
    '   the maximum possible entropy is 8. A value > 7.8 indicates high randomness,\n' ...
    '   making the encrypted image appear as noise and resisting statistical attacks.\n\n']);
fprintf(['2. PSNR & SSIM (Encryption): These metrics compare the encrypted image to the original.\n' ...
    '   - A LOW PSNR (< 10 dB) and a LOW SSIM (< 0.1) are desired. They prove that the\n' ...
    '     encrypted image is significantly different from the original, and its structure\n' ...
    '     is completely obscured.\n\n']);
fprintf(['3. NPCR & UACI (Differential Resistance): These test the "avalanche effect". We change\n' ...
    '   one pixel in the original image and see how much the encrypted image changes.\n' ...
    '   - NPCR > 99.6%% (ideal) means almost all pixels changed.\n' ...
    '   - UACI ~ 33.46%% (ideal) means the changed pixels have significantly different values.\n' ...
    '   High NPCR and ideal UACI prove strong resistance to differential attacks.\n\n']);
fprintf(['4. PSNR & SSIM (Decryption): These metrics compare the decrypted image to the original.\n' ...
    '   - A HIGH PSNR (> 30 dB) and a HIGH SSIM (close to 1.0) are desired. They prove\n' ...
    '     that the original image was recovered with excellent quality and minimal data loss.\n\n']);


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

% NPCR and UACI Calculation Function
function [npcr, uaci] = calculate_npcr_uaci(C1, C2)
    % This function calculates the Number of Pixels Change Rate (NPCR) and
    % Unified Average Changing Intensity (UACI) between two uint8 images.
    [M, N] = size(C1);

    % --- NPCR Calculation ---
    % Count the number of pixels that are different between the two images.
    num_diff_pixels = sum(C1(:) ~= C2(:));
    % Express as a percentage of the total number of pixels.
    npcr = (num_diff_pixels / (M * N)) * 100;

    % --- UACI Calculation ---
    % Calculate the normalized absolute difference between the two images.
    % Convert to double for subtraction to avoid integer overflow issues.
    sum_abs_diff = sum(abs(double(C1(:)) - double(C2(:))));
    % Normalize by the max possible sum of differences (M*N*255) and express as a percentage.
    uaci = (sum_abs_diff / (M * N * 255)) * 100;
end
