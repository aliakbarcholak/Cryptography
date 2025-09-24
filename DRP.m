clc; clear; close all;

%% ========= إعدادات عامة =========
folder_path = 'F:\رسالة ماجستير\البيانات المطلوبة\DataBase\png';      % <-- عدّل المسار
passphrase  = 'Strong@Pass-Example!'; % <-- عدّل عبارة المرور إن رغبت

image_files = dir(fullfile(folder_path, '*.PNG'));
if isempty(image_files)
    error('No images found. Check folder path/extension.');
end

% مصفوفات تخزين القياسات
num_images = numel(image_files);
MSE_encrypted  = zeros(1, num_images);
MSE_decrypted  = zeros(1, num_images);
PSNR_enc       = zeros(1, num_images);
PSNR_dec       = zeros(1, num_images);
SSIM_enc       = zeros(1, num_images);
SSIM_dec       = zeros(1, num_images);
ENT_enc        = zeros(1, num_images);
NPCR_vals      = zeros(1, num_images);
UACI_vals      = zeros(1, num_images);

for i = 1:num_images
    %% ======== قراءة وتحضير الصورة ========
    img_path = fullfile(folder_path, image_files(i).name);
    I0 = imread(img_path);
    I0 = imresize(I0, [1024 1024]); % توحيد الحجم
    I0 = im2gray(I0);
    Iu8 = uint8(I0);
    [R,C] = size(Iu8);

    %% ======== اشتقاق مفاتيح/سلاسل من العبارة والاسم ========
    key_material = [passphrase '|' image_files(i).name];
    key = sha256_bytes(key_material);          % 32 بايت
    need_len = R*C*8;
    [seq1, seq2] = chaotic_coupled_logsin(key, need_len);

    % فهارس Permutation
    p_seq = seq1(1:R);  [~, row_idx] = sort(p_seq);  inv_row(row_idx) = 1:R;
    q_seq = seq2(1:C);  [~, col_idx] = sort(q_seq);  inv_col(col_idx) = 1:C;

    % Keystream بايتات
    ks = seq1(R+1 : R + C*R);
    ks = uint8( floor( mod(ks*1e14, 256) ) );
    ks = reshape(ks, R, C);

    % S-Box ديناميكية مشتقة من الفوضى (Permutation لـ 0..255) وعكسها
    [S, Sinv] = sbox_from_seq(seq2(R+1:R+256));

    % IVs من المفتاح
    iv1 = uint8(bitxor(key(1),  key(16)));  % مبدئيات التغذية الراجعة
    iv2 = uint8(bitxor(key(8),  key(24)));

    %% ======== التشفير في المجال المكاني (Permutation → Subst → 2-Pass Diffusion) ========
    P = Iu8(row_idx, :);           % Permute rows
    P = P(:, col_idx);             % Permute cols

    % استبدال (Substitution) عبر S-Box
    Psub = S(double(P)+1);         % قيم 1..256
    Psub = uint8(Psub-1);

    % انتشار قوي: تمرير أمامي ثم رجعي مع XOR + جمع نمطي + تغذية راجعة
    C1  = diffuse_pass(Psub, ks, iv1, 'forward');
    Csp = diffuse_pass(C1,  rot90(ks,2), iv2, 'backward'); % أقوى حساسية

    %% ======== DRPE لطمس طيفي إضافي (قابل للعكس) ========
    P1 = exp(1i * 2*pi * reshape(mod(seq1(end-R*C+1:end),1), R, C));
    P2 = exp(1i * 2*pi * reshape(mod(seq2(end-R*C+1:end),1), R, C));
    F1 = fft2(double(Csp)/255);
    E  = fft2(F1 .* P1) .* P2;                  % مُعقّد
    % قياس Entropy على |E|
    E_mag_u8 = uint8( round( 255 * mat2gray(abs(E)) ) );
    ENT_enc(i) = image_entropy(E_mag_u8);

    %% ======== فك التشفير للتحقق من SSIM/PSNR ========
    F1_rec = ifft2(E ./ P2);
    Dsp_d  = ifft2(F1_rec ./ P1);
    Dsp    = uint8( round( min(max(Dsp_d,0),1) * 255 ) );

    % عكس الانتشار ثنائي المرور
    C1_rec  = undiffuse_pass(Dsp,  rot90(ks,2), iv2, 'backward');
    Psubrec = undiffuse_pass(C1_rec, ks, iv1, 'forward');

    % عكس الاستبدال والـ Permutation
    Prec_u8  = uint8(Sinv(double(Psubrec)+1)-1);
    I_rec_u8 = Prec_u8(inv_row, :);
    I_rec_u8 = I_rec_u8(:, inv_col);

    %% ======== القياسات ========
    % SSIM/PSNR بين الأصلية والمشفّرة (منخفضة) – نقيس على Csp (المجال المكاني) لملاءمة NPCR/UACI
    PSNR_enc(i) = psnr(Csp, Iu8);
    SSIM_enc(i) = ssim(Csp, Iu8);

    % SSIM/PSNR بين الأصلية والمفكوكة (مرتفعة)
    PSNR_dec(i) = psnr(I_rec_u8, Iu8);
    SSIM_dec(i) = ssim(I_rec_u8, Iu8);

    % MSE
    MSE_encrypted(i) = immse(im2double(Iu8), im2double(Csp));
    MSE_decrypted(i) = immse(im2double(Iu8), im2double(I_rec_u8));

    % NPCR/UACI على النص المُعمّى المكاني Csp
    I0_f = Iu8; I0_f(1,1) = bitxor(I0_f(1,1), uint8(1)); % تغيير بكسل واحد
    Pf    = I0_f(row_idx, :); Pf = Pf(:, col_idx);
    Psubf = S(double(Pf)+1); Psubf = uint8(Psubf-1);
    C1f   = diffuse_pass(Psubf, ks, iv1, 'forward');
    Csp_f = diffuse_pass(C1f,  rot90(ks,2), iv2, 'backward');

    [NPCR_vals(i), UACI_vals(i)] = npcr_uaci(Csp, Csp_f);

    %% ======== طباعة وعرض مختصر ========
    fprintf('Image %d: %s\n', i, image_files(i).name);
    fprintf('  NPCR: %.2f%% | UACI: %.2f%% | SSIM(enc): %.4f\n', NPCR_vals(i), UACI_vals(i), SSIM_enc(i));
    fprintf('  SSIM(dec): %.4f | PSNR(dec): %.2f dB | Entropy(|E|): %.4f\n\n', SSIM_dec(i), PSNR_dec(i), ENT_enc(i));

    figure('Name', sprintf('Image %d', i), 'NumberTitle','off');
    subplot(1,3,1); imshow(Iu8);       title('Original');
    subplot(1,3,2); imshow(Csp);       title('Encrypted (spatial)');
    subplot(1,3,3); imshow(I_rec_u8);  title('Decrypted');
    drawnow;
end

%% ======== ملخص النتائج ========
disp('--- NPCR (%) ---');                 disp(NPCR_vals);
disp('--- UACI (%) ---');                 disp(UACI_vals);
disp('--- SSIM(enc) ---');                disp(SSIM_enc);
disp('--- SSIM(dec) ---');                disp(SSIM_dec);
disp('--- PSNR(dec) ---');                disp(PSNR_dec);
disp('--- MSE Encrypted ---');            disp(MSE_encrypted);
disp('--- MSE Decrypted ---');            disp(MSE_decrypted);
disp('--- Entropy of |E| ---');           disp(ENT_enc);

%% ================== الدوال المساعدة ==================

function H = image_entropy(Iu8)
    counts = imhist(Iu8, 256);
    p = counts / sum(counts);
    p = p(p>0);
    H = -sum(p .* log2(p));
end

function [NPCR, UACI] = npcr_uaci(C1, C2)
    C1 = uint8(C1); C2 = uint8(C2);
    [R,CC] = size(C1);
    D = C1 ~= C2;
    NPCR = 100 * sum(D(:)) / (R*CC);
    UACI = 100 * mean( abs(double(C1(:)) - double(C2(:))) / 255 );
end

function key = sha256_bytes(str)
    md = java.security.MessageDigest.getInstance('SHA-256');
    md.update(uint8(str));
    key = typecast(md.digest, 'uint8');
end

function [x, y] = chaotic_coupled_logsin(key_bytes, L)
% مُولِّد فوضوي مُقترن Logistic–Sine مستهلَك من مفتاح 256-بت
% ينتج سلسلتين x,y بطول L بقيم ضمن (0,1)
    u32 = typecast(key_bytes, 'uint32');
    s1 = double(bitxor(u32(1), u32(3))) / double(intmax('uint32')) + 0.11;
    s2 = double(bitxor(u32(2), u32(4))) / double(intmax('uint32')) + 0.17;
    a  = 3.8 + 0.19*frac(s1*7.1);
    b  = 3.8 + 0.19*frac(s2*9.3);
    x0 = frac( sum(double(u32(5:6)))/2^32 + 0.123456 );
    y0 = frac( sum(double(u32(7:8)))/2^32 + 0.654321 );
    burn = 1000;
    [x0,y0] = iter(x0, y0, a, b, burn);
    x = zeros(L,1); y = zeros(L,1);
    xn = x0; yn = y0;
    for k = 1:L
        [xn, yn] = step_map(xn, yn, a, b);
        xn = frac(xn + 0.5*frac(s1*xn + s2*yn));
        yn = frac(yn + 0.5*frac(s2*yn + s1*xn));
        x(k) = frac(xn); y(k) = frac(yn);
    end

    function [xo, yo] = iter(xi, yi, aa, bb, n)
        xo = xi; yo = yi;
        for t = 1:n
            [xo, yo] = step_map(xo, yo, aa, bb);
            xo = frac(xo + 0.3*frac(s1*yo));
            yo = frac(yo + 0.3*frac(s2*xo));
        end
    end
    function [xn1, yn1] = step_map(xn, yn, aa, bb)
        xn1 = sin(pi*( aa*xn*(1-xn) + 4*yn*(1-yn) ));
        yn1 = sin(pi*( bb*yn*(1-yn) + 4*xn1*(1-xn1) ));
        xn1 = frac(0.5*(xn1+1));
        yn1 = frac(0.5*(yn1+1));
    end
    function z = frac(v); z = v - floor(v); end
end

function out = diffuse_pass(P, ks, iv, dir)
% نشر قوي بتمرير واحد (أمامي/رجعي) مع XOR + جمع نمطي + تغذية راجعة
% out = C ؛ P و ks uint8 بنفس الأبعاد
    [R,CC] = size(P);
    out = zeros(R,CC,'uint8');
    if strcmpi(dir,'forward')
        prevC = iv; prevP = uint8(bitxor(iv, 171));
        for r = 1:R
            for c = 1:CC
                t = bitxor(P(r,c), ks(r,c));             % XOR
                s = uint16(t) + uint16(prevC) + uint16(prevP); % جمع نمطي
                out(r,c) = uint8( mod(s, 256) );
                prevC = out(r,c);
                prevP = P(r,c);
            end
        end
    else % backward
        prevC = iv; prevP = uint8(bitxor(iv, 219));
        for r = R:-1:1
            for c = CC:-1:1
                t = bitxor(P(r,c), ks(r,c));
                s = uint16(t) + uint16(prevC) + uint16(prevP);
                out(r,c) = uint8( mod(s, 256) );
                prevC = out(r,c);
                prevP = P(r,c);
            end
        end
    end
end

function P = undiffuse_pass(Cin, ks, iv, dir)
% عكس انتشار diffuse_pass بدقّة (نسخة مُصحّحة)
% Cin, ks: uint8 وبنفس الأبعاد
    [R,CC] = size(Cin);
    P = zeros(R,CC,'uint8');

    if strcmpi(dir,'forward')
        prevC = iv;
        prevP = uint8(bitxor(iv, 171));
        for r = 1:R
            for c = 1:CC
                s = int16(Cin(r,c)) - int16(prevC) - int16(prevP);
                t = uint8( mod(s, 256) );
                P(r,c) = bitxor(t, ks(r,c));
                prevC = Cin(r,c);
                prevP = P(r,c);
            end
        end
    else % 'backward'
        prevC = iv;
        prevP = uint8(bitxor(iv, 219));
        for r = R:-1:1
            for c = CC:-1:1
                s = int16(Cin(r,c)) - int16(prevC) - int16(prevP);
                t = uint8( mod(s, 256) );
                P(r,c) = bitxor(t, ks(r,c));
                prevC = Cin(r,c);
                prevP = P(r,c);
            end
        end
    end
end

function [S, Sinv] = sbox_from_seq(seq256)
% توليد S-Box ديناميكية قابلة للعكس من تسلسل فوضوي
    [~, ord] = sort(mod(seq256(:),1), 'ascend');
    S = zeros(256,1,'uint16');
    S(ord) = uint16(0:255)+1;     % قيم 1..256
    % العكس
    Sinv = zeros(256,1,'uint16');
    Sinv(S) = uint16(0:255)+1;
end
