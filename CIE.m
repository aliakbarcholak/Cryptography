%% CIE_SM_ALI — Secure Chaotic Image Encryption (Permutation + Two-Pass Diffusion)
% الأهداف/المؤشرات:
% - Entropy للصورة المشفّرة غالبًا ≥ 7.8
% - PSNR > 30 dB و SSIM ≥ 0.9 بعد فك التشفير (استرجاع صحيح)
% - مقاومة تفاضلية قوية: NPCR > 99% و UACI ≈ 33%
% - توثيق: PSNR_Encrypted و SSIM_Encrypted و MSE_Encrypted
%
% المكوّنات:
% - Permutation: Lorenz (مؤشر/فهرسة سطر/عمود)
% - Diffusion: Logistic Keystream + تمريرين (CBC أمامي + CBC رجعي) مع IVs مستقلة
% - بذور ومفاتيح مستخرجة من SHA-256 لضمان الثبات
%
% ملاحظات تشغيل:
% - عدّل المسار والامتداد حسب بياناتك.
% - الكود يحوّل الصور إلى رمادية (للتبسيط). يمكن توسعته للـRGB لاحقًا.

clc; clear; close all;

%% --------------------------- إعدادات عامة ---------------------------
folder_path = 'F:\DataBase\PNG';      % <-- عدّل المسار
ext_list    = {'*.png'};  % يمكنك تعديل/تقليل الأنواع

% جمع الملفات
image_files = [];
for e = 1:numel(ext_list)
    image_files = [image_files; dir(fullfile(folder_path, ext_list{e}))]; %#ok<AGROW>
end
if isempty(image_files)
    error('لا توجد صور في المجلد المحدد. تأكد من المسار ونوع الامتداد.');
end

% إعادة تحجيم اختياري (ضع [] لتعطيله)
targetSize = [1024, 1024];

% مفتاح/سر المستخدم
userSecret = 'myStrong!Secret#Key2025';

% Lorenz (للمفاتيح/الفهارس)
lorenzParams.sigma = 10; 
lorenzParams.rho   = 28; 
lorenzParams.beta  = 8/3;
lorenzParams.h     = 0.01;
lorenzParams.burn  = 1000;

% Logistic (لمفتاح الانتشار)
logisticParams.r    = 3.999;     % أعلى قليلاً لزيادة العشوائية
logisticParams.burn = 2000;      % حرق أطول لتقليل الارتباطات

% حاويات النتائج
num_images    = numel(image_files);
Entropy_enc   = zeros(1, num_images);
PSNR_dec      = zeros(1, num_images);
SSIM_dec      = zeros(1, num_images);
NPCR_scores   = zeros(1, num_images);
UACI_scores   = zeros(1, num_images);
PSNR_enc      = zeros(1, num_images);   % PSNR بين الأصل والمشفّر
SSIM_enc      = zeros(1, num_images);   % SSIM بين الأصل والمشفّر
MSE_enc       = zeros(1, num_images);   % NEW: MSE (Orig vs Enc)
MSE_dec       = zeros(1, num_images);   % NEW: MSE (Dec  vs Orig)

%% -------------------------- حلقة معالجة الصور --------------------------
for i = 1:num_images
    % ----- القراءة والتحضير -----
    img_path = fullfile(folder_path, image_files(i).name);
    I = imread(img_path);
    if ~isempty(targetSize), I = imresize(I, targetSize); end
    I = im2gray(I);                % معالجة رمادية (قناة واحدة)
    Iu8 = im2uint8(I);
    [M,N] = size(Iu8);

    % ----- بذور فوضوية آمنة من هاش ثابت للصورة + السر -----
    seedBytes = sha256_bytes([userSecret ':' num2str(i) ':' char(Iu8(:).')]);
    [x0,y0,z0,logistic_x0] = seeds_from_hash(seedBytes);
    logistic_x0 = sanitize01(logistic_x0);
    x0 = sanitizeReal(x0); y0 = sanitizeReal(y0); z0 = sanitizeReal(z0);

    % ----- Permutation (Lorenz) -----
    [rowPerm, colPerm] = lorenz_permutation_indices(M, N, [x0 y0 z0], lorenzParams);
    Iperm = permute2d(Iu8, rowPerm, colPerm);

    % ----- Diffusion (Logistic + تمريرين لرفع UACI/NPCR) -----
    K = logistic_keystream_uint8(M, N, logisticParams.r, logistic_x0, logisticParams.burn);
    Krot = rot90(K, 2);                                % مفتاح مقلوب لتمرير رجعي

    % IVs مشتقة بثبات من الهاش لضمان قابلية الاسترجاع
    iv1 = uint8(bitxor(seedBytes(1),  seedBytes(16)));
    iv2 = uint8(bitxor(seedBytes(8),  seedBytes(24)));

    % تشفير بتمريرين (CBC أمامي + CBC رجعي)
    C = diffuse_two_pass_uint8(Iperm, K, Krot, iv1, iv2);   % مشفّرة (uint8)

    % ----- فك التشفير للتحقق (عكس التمريرين ثم عكس الـPermutation) -----
    Iperm_rec = inverse_diffuse_two_pass_uint8(C, K, Krot, iv1, iv2);
    Irec      = inverse_permute2d(Iperm_rec, rowPerm, colPerm);

    % ----- القياسات -----
    Entropy_enc(i) = entropy(C);
    PSNR_dec(i)    = psnr(Irec, Iu8);
    SSIM_dec(i)    = ssim(Irec, Iu8);
    PSNR_enc(i)    = psnr(C, Iu8);        % بين الأصل والمشفّر (منخفض = تشفير أقوى)
    SSIM_enc(i)    = ssim(C, Iu8);        % بين الأصل والمشفّر (منخفض = تشفير أقوى)
    MSE_enc(i)     = immse(im2double(C),    im2double(Iu8));  % NEW
    MSE_dec(i)     = immse(im2double(Irec), im2double(Iu8));  % NEW

    % ----- هجوم تفاضلي (قلب بكسل واحد في الوسط تقريبًا) -----
    Iu8_pert = Iu8;
    rr = max(1, floor(M/2)); cc = max(1, floor(N/2));
    Iu8_pert(rr,cc) = uint8(mod(double(Iu8_pert(rr,cc)) + 1, 256));

    Iperm2 = permute2d(Iu8_pert, rowPerm, colPerm);
    C2     = diffuse_two_pass_uint8(Iperm2, K, Krot, iv1, iv2);
    [NPCR_scores(i), UACI_scores(i)] = npcr_uaci(C, C2);

    % ----- عرض اختياري -----
    figure('Name', sprintf('Image %d : %s', i, image_files(i).name), 'NumberTitle','off');
    subplot(1,3,1); imshow(Iu8, []); title('Original');
    subplot(1,3,2); imshow(C,   []); title('Encrypted');
    subplot(1,3,3); imshow(Irec,[]); title('Decrypted');

    % ----- طباعة نتائج الصورة -----
    fprintf('--- Image %d: %s ---\n', i, image_files(i).name);
    fprintf('Entropy(Encrypted)            : %.4f (هدف ≥ 7.8)\n', Entropy_enc(i));
    fprintf('PSNR_Encrypted (Orig vs Enc)  : %.4f dB (منخفض = تشفير جيد)\n', PSNR_enc(i));
    fprintf('SSIM_Encrypted (Orig vs Enc)  : %.4f   (منخفض = تشفير جيد)\n', SSIM_enc(i));
    fprintf('MSE_Encrypted (Orig vs Enc)   : %.6f   (أكبر = تشفير أقوى)\n', MSE_enc(i));
    fprintf('PSNR(Decrypted vs Original)   : %.4f dB (هدف > 30 dB)\n', PSNR_dec(i));
    fprintf('SSIM(Decrypted vs Original)   : %.4f   (هدف ≥ 0.9)\n', SSIM_dec(i));
    fprintf('MSE(Decrypted vs Original)    : %.6f   (أقرب للصفر = أفضل)\n', MSE_dec(i));
    fprintf('NPCR (C vs C2)                : %.4f %% (هدف > 99%%)\n', NPCR_scores(i));
    fprintf('UACI (C vs C2)                : %.4f %% (هدف ≈ 33%%)\n\n', UACI_scores(i));
end

%% -------------------------- ملخّص نهائي --------------------------
fprintf('==================== SUMMARY ====================\n');
fprintf('Images processed            : %d\n', num_images);
fprintf('Entropy (enc)               : mean=%.4f  min=%.4f  max=%.4f\n', mean(Entropy_enc), min(Entropy_enc), max(Entropy_enc));
fprintf('PSNR_Encrypted (dB)         : mean=%.4f  min=%.4f  max=%.4f\n', mean(PSNR_enc),    min(PSNR_enc),    max(PSNR_enc));
fprintf('SSIM_Encrypted              : mean=%.4f  min=%.4f  max=%.4f\n', mean(SSIM_enc),    min(SSIM_enc),    max(SSIM_enc));
fprintf('MSE_Encrypted               : mean=%.6f  min=%.6f  max=%.6f\n', mean(MSE_enc),     min(MSE_enc),     max(MSE_enc));
fprintf('PSNR (dec) dB               : mean=%.4f  min=%.4f  max=%.4f\n', mean(PSNR_dec),     min(PSNR_dec),    max(PSNR_dec));
fprintf('SSIM (dec)                  : mean=%.4f  min=%.4f  max=%.4f\n', mean(SSIM_dec),     min(SSIM_dec),    max(SSIM_dec));
fprintf('MSE(Decrypted)              : mean=%.6f  min=%.6f  max=%.6f\n', mean(MSE_dec),     min(MSE_dec),     max(MSE_dec));
fprintf('NPCR  %%                     : mean=%.4f  min=%.4f  max=%.4f\n', mean(NPCR_scores),  min(NPCR_scores), max(NPCR_scores));
fprintf('UACI  %%                     : mean=%.4f  min=%.4f  max=%.4f\n', mean(UACI_scores),  min(UACI_scores), max(UACI_scores));
fprintf('=================================================\n');

%% ========================= الدوال المساعدة =========================
function x = sanitize01(x)
    if isstring(x) || ischar(x), x = double(uint8(x(:))).'; end
    x = double(x(:)); if isempty(x), x = 0.123456; else, x = x(1); end
    if ~isfinite(x), x = 0.654321; end
    x = x - floor(x);
    if x <= 0, x = eps; end
    if x >= 1, x = 1 - eps; end
end

function x = sanitizeReal(x)
    if isstring(x) || ischar(x), x = double(uint8(x(:))).'; end
    x = double(x(:)); if isempty(x), x = 0; else, x = x(1); end
    if ~isfinite(x), x = 0; end
end

function bytes = sha256_bytes(str)
    md = java.security.MessageDigest.getInstance('SHA-256');
    md.update(uint8(str));
    digest = typecast(md.digest(),'uint8');
    bytes = uint8(digest);
end

function [x0,y0,z0,logx0] = seeds_from_hash(hashBytes)
    assert(numel(hashBytes) >= 16, 'Hash too short');
    w = double(hashBytes(:)).';
    s1 = sum(w(1:4)  .* [1 256 65536 16777216]);
    s2 = sum(w(5:8)  .* [1 256 65536 16777216]);
    s3 = sum(w(9:12) .* [1 256 65536 16777216]);
    s4 = sum(w(13:16).* [1 256 65536 16777216]);

    u1 = (s1+1)/ (2^32); u2 = (s2+1)/ (2^32);
    u3 = (s3+1)/ (2^32); u4 = (s4+1)/ (2^32);

    x0 = -15 + 30*u1;
    y0 = -15 + 30*u2;
    z0 =  10 + 30*u3;
    logx0 = max(min(u4, 1 - eps), eps);
end

function [rowPerm, colPerm] = lorenz_permutation_indices(M, N, initXYZ, p)
    x = initXYZ(1); y = initXYZ(2); z = initXYZ(3);
    L = max(M,N) + p.burn + 10;
    seq = zeros(L,3);
    for k = 2:L
        dx = p.sigma*(y - x);
        dy = x*(p.rho - z) - y;
        dz = x*y - p.beta*z;
        x = x + p.h*dx;
        y = y + p.h*dy;
        z = z + p.h*dz;
        seq(k,:) = [x,y,z];
    end
    S = seq(p.burn+1:end,:);
    S = S - floor(S);
    Sr = S(1:M,1);
    Sc = S(1:N,2);
    [~, rowPerm] = sort(Sr, 'ascend');
    [~, colPerm] = sort(Sc, 'ascend');
end

function K = logistic_keystream_uint8(M, N, r, x0, burn)
    r  = double(r); if ~isscalar(r),  r  = r(1);  end
    x0 = sanitize01(x0);
    L = M*N + burn; if L < 2, L = 2; end
    x = zeros(L,1);
    x(1) = x0;
    for k = 2:L
        x(k) = r .* x(k-1) .* (1 - x(k-1));
        if x(k) <= 0, x(k) = eps; elseif x(k) >= 1, x(k) = 1 - eps; end
    end
    stream = x(burn+1:end);
    K = uint8(floor(stream * 256));
    K = reshape(K, M, N);
end

function J = permute2d(I, rowPerm, colPerm)
    J = I(rowPerm, :);
    J = J(:, colPerm);
end

function I = inverse_permute2d(J, rowPerm, colPerm)
    Itemp = zeros(size(J), 'like', J);
    I = zeros(size(J), 'like', J);
    Itemp(rowPerm,:) = J;   % إعادة الصفوف
    I(:, colPerm)    = Itemp; % ثم الأعمدة
end

% ===================== انتشار ثنائي التمرير (لرفع UACI/NPCR) =====================
function C = diffuse_two_pass_uint8(P, Kf, Kr, iv1, iv2)
    % تمرير أمامي (CBC)
    [M,N] = size(P);
    C1 = zeros(M,N,'uint8');
    prevC = iv1; prevP = bitxor(iv1, uint8(171));
    for r = 1:M
        for c = 1:N
            x = bitxor(P(r,c), Kf(r,c));
            s = uint16(x) + uint16(prevC) + uint16(prevP);   % mod 256
            C1(r,c) = uint8(mod(s, 256));
            prevP = P(r,c); prevC = C1(r,c);
        end
    end

    % تمرير رجعي (CBC) بمفتاح مختلف (Kr)
    C2 = zeros(M,N,'uint8');
    prevC = iv2; prevP = bitxor(iv2, uint8(219));
    for r = M:-1:1
        for c = N:-1:1
            x = bitxor(C1(r,c), Kr(r,c));
            s = uint16(x) + uint16(prevC) + uint16(prevP);
            C2(r,c) = uint8(mod(s, 256));
            prevP = C1(r,c); prevC = C2(r,c);
        end
    end
    C = C2;
end

function P = inverse_diffuse_two_pass_uint8(C2, Kf, Kr, iv1, iv2)
    % عكس التمرير الرجعي أولاً
    [M,N] = size(C2);
    C1 = zeros(M,N,'uint8');
    prevC = iv2; prevP = bitxor(iv2, uint8(219));
    for r = M:-1:1
        for c = N:-1:1
            % C2 = (C1 ⊕ Kr) + prevC + prevP (mod 256)
            t = int16(C2(r,c)) - int16(prevC) - int16(prevP);
            x = uint8(mod(t, 256));
            C1(r,c) = bitxor(x, Kr(r,c));
            prevP = C1(r,c); prevC = C2(r,c);
        end
    end

    % ثم عكس التمرير الأمامي
    P = zeros(M,N,'uint8');
    prevC = iv1; prevP = bitxor(iv1, uint8(171));
    for r = 1:M
        for c = 1:N
            % C1 = (P ⊕ Kf) + prevC + prevP (mod 256)
            t = int16(C1(r,c)) - int16(prevC) - int16(prevP);
            x = uint8(mod(t, 256));
            P(r,c) = bitxor(x, Kf(r,c));
            prevP = P(r,c); prevC = C1(r,c);
        end
    end
end

% ============================ NPCR & UACI ============================
function [NPCRp, UACIp] = npcr_uaci(C1, C2)
    assert(isequal(size(C1), size(C2)), 'أبعاد غير متساوية');
    [M,N] = size(C1);
    D = C1 ~= C2;
    NPCRp = 100 * sum(D(:)) / (M*N);
    UACIp = 100 * mean(abs(double(C1(:)) - double(C2(:))) / 255);
end
