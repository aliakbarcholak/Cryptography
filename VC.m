clc; clear; close all;

%% ========= إعدادات عامة =========
folder_path = 'F:\DataBase\PNG';      % <-- عدّل المسار
passphrase  = 'Strong@Pass-Example!'; % <-- يمكن تغييره

image_files = dir(fullfile(folder_path, '*.PNG'));  % يمكنك تغيير الامتداد
if isempty(image_files)
    error('No images found in the folder. Make sure the folder path is correct and contains images.');
end

% مصفوفات القياسات لكل صورة
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
    %% ======== قراءة وتحضير الصورة (RGB أو رمادية) ========
    img_path = fullfile(folder_path, image_files(i).name);
    Iraw = imread(img_path);
    Iraw = imresize(Iraw, [1024 1024]);    % حجم موحّد
    Iu8  = im2uint8(Iraw);                 % إلى 0..255
    [R, CC, CH] = size(Iu8);               % CH=1 للرمادية، 3 للملوّنة

    % حاويات للمخرجات لكل قناة
    Csp      = zeros(R,CC,CH,'uint8');   % النص المعمّى المكاني (للتقييم وتقسيم الأسهم)
    Irec_u8  = zeros(R,CC,CH,'uint8');   % صورة مفكوكة
    ENT_ch   = zeros(CH,1);              % Entropy القنوية على |E|

    %% ======== تعمية كل قناة مستقلة (مفاتيح مشتقة من الاسم والقناة) ========
    for ch = 1:CH
        Ich = Iu8(:,:,ch);

        % مفاتيح/سلاسل من عبارة المرور + اسم الملف + رقم القناة
        key_material = [passphrase '|' image_files(i).name '|ch=' num2str(ch)];
        key = sha256_bytes(key_material);                 % 32 بايت
        need_len = R*CC*8;
        [seq1, seq2] = chaotic_coupled_logsin(key, need_len);

        % فهارس Permutation
        p_seq = seq1(1:R);  [~, row_idx] = sort(p_seq, 'ascend');  inv_row(row_idx) = 1:R;
        q_seq = seq2(1:CC); [~, col_idx] = sort(q_seq, 'ascend');  inv_col(col_idx) = 1:CC;

        % Keystream (R×C)
        ks = seq1(R+1 : R + R*CC);
        ks = uint8( floor( mod(ks*1e14, 256) ) );
        ks = reshape(ks, R, CC);

        % S-Box ديناميكية وعكسها
        [S, Sinv] = sbox_from_seq(seq2(R+1:R+256));

        % IVs للتغذية الراجعة
        iv1 = uint8(bitxor(key(1),  key(16)));
        iv2 = uint8(bitxor(key(8),  key(24)));

        % ======== التشفير المكاني: Permute → Substitution → 2-Pass Diffusion ========
        P = Ich(row_idx, :);
        P = P(:, col_idx);

        Psub = S(double(P)+1);             % 1..256
        Psub = uint8(Psub - 1);

        C1          = diffuse_pass(Psub, ks, iv1, 'forward');
        Csp(:,:,ch) = diffuse_pass(C1,  rot90(ks,2), iv2, 'backward');

        % ======== DRPE (اختياري للطمس وقياس Entropy) ========
        % لا تُستخدم DRPE في بناء الأسهم، لكنها تقدّم طمساً طيفياً ومؤشر Entropy قوي.
        P1 = exp(1i * 2*pi * reshape(mod(seq1(end-R*CC+1:end),1), R, CC));
        P2 = exp(1i * 2*pi * reshape(mod(seq2(end-R*CC+1:end),1), R, CC));
        F1 = fft2(double(Csp(:,:,ch))/255);
        E  = fft2(F1 .* P1) .* P2;                         % مُعقّد
        E_mag_u8 = uint8( round( 255 * mat2gray(abs(E)) ) );
        ENT_ch(ch) = image_entropy(E_mag_u8);

        % ======== فك التشفير للتحقق ========
        F1_rec = ifft2(E ./ P2);
        Dsp_d  = ifft2(F1_rec ./ P1);
        Dsp    = uint8( round( min(max(Dsp_d,0),1) * 255 ) );

        C1_rec   = undiffuse_pass(Dsp,  rot90(ks,2), iv2, 'backward');
        Psub_rec = undiffuse_pass(C1_rec, ks, iv1, 'forward');

        Prec_u8        = uint8(Sinv(double(Psub_rec)+1)-1);
        % عكس الترتيب: الصفوف ثم الأعمدة (بدون فهرسة متسلسلة غير مدعومة)
        tmp = Prec_u8(inv_row, :);  % استرجاع الصفوف
        tmp = tmp(:, inv_col);      % ثم الأعمدة
        Irec_u8(:,:,ch) = tmp;
    end

    % ======== بناء الأسهم من Csp (Visual Sharing) ========
    share1 = zeros(R,CC,CH,'uint8');
    share2 = zeros(R,CC,CH,'uint8');
    for ch = 1:CH
        % Keystream للأسهم من seq2 (جزء مختلف) لكل قناة
        key_material = [passphrase '|' image_files(i).name '|share|ch=' num2str(ch)];
        key_sh = sha256_bytes(key_material);
        [sqA, ~] = chaotic_coupled_logsin(key_sh, R*CC + 10);
        ks_share = uint8( floor( mod(sqA(1:R*CC)*1e14, 256) ) );
        ks_share = reshape(ks_share, R, CC);

        share1(:,:,ch) = ks_share;                         % سهم 1
        share2(:,:,ch) = bitxor(Csp(:,:,ch), ks_share);    % سهم 2 بحيث share1 XOR share2 = Csp
    end
    Csp_from_shares = bitxor(share1, share2);              % تحقق: يجب أن يساوي Csp

    % ======== قياسات الجودة/الأمان ========
    PSNR_enc(i) = psnr(Csp, Iu8);
    SSIM_enc(i) = ssim(Csp, Iu8);

    PSNR_dec(i) = psnr(Irec_u8, Iu8);
    SSIM_dec(i) = ssim(Irec_u8, Iu8);

    MSE_encrypted(i) = immse(im2double(Iu8), im2double(Csp));
    MSE_decrypted(i) = immse(im2double(Iu8), im2double(Irec_u8));

    ENT_enc(i) = mean(ENT_ch);

    % NPCR/UACI: قلب بكسل واحد ثم إعادة السلسلة حتى Csp_f
    Iflip = Iu8;
    Iflip(1,1,1) = bitxor(Iflip(1,1,1), uint8(1));   % قلب LSB لقناة أولى
    Csp_f = zeros(R,CC,CH,'uint8');
    for ch = 1:CH
        Ich = Iflip(:,:,ch);
        key_material = [passphrase '|' image_files(i).name '|ch=' num2str(ch)];
        key = sha256_bytes(key_material);
        need_len = R*CC*8;
        [seq1, seq2] = chaotic_coupled_logsin(key, need_len);

        p_seq = seq1(1:R);  [~, row_idx] = sort(p_seq, 'ascend');
        q_seq = seq2(1:CC); [~, col_idx] = sort(q_seq, 'ascend');

        ks = seq1(R+1 : R + R*CC);
        ks = uint8( floor( mod(ks*1e14, 256) ) );
        ks = reshape(ks, R, CC);

        [S, ~] = sbox_from_seq(seq2(R+1:R+256));
        iv1 = uint8(bitxor(key(1),  key(16)));
        iv2 = uint8(bitxor(key(8),  key(24)));

        P = Ich(row_idx, :); P = P(:, col_idx);
        Psub = S(double(P)+1); Psub = uint8(Psub - 1);
        C1 = diffuse_pass(Psub, ks, iv1, 'forward');
        Csp_f(:,:,ch) = diffuse_pass(C1, rot90(ks,2), iv2, 'backward');
    end
    [NPCR_vals(i), UACI_vals(i)] = npcr_uaci(Csp, Csp_f);

    %% ======== عرض مختصر ========
    fprintf('Image %d: %s\n', i, image_files(i).name);
    fprintf('  NPCR: %.2f%% | UACI: %.2f%% | SSIM(enc): %.4f\n', NPCR_vals(i), UACI_vals(i), SSIM_enc(i));
    fprintf('  SSIM(dec): %.4f | PSNR(dec): %.2f dB | Entropy(|E|): %.4f\n\n', SSIM_dec(i), PSNR_dec(i), ENT_enc(i));

    figure('Name', sprintf('Image %d', i), 'NumberTitle','off');
    subplot(2,3,1); imshow(Iu8);                title('Original');
    subplot(2,3,2); imshow(Csp);                title('Encrypted (spatial)');
    subplot(2,3,3); imshow(Irec_u8);            title('Decrypted');
    subplot(2,3,4); imshow(share1);             title('Share 1');
    subplot(2,3,5); imshow(share2);             title('Share 2');
    subplot(2,3,6); imshow(Csp_from_shares);    title('XOR(share1,share2)');
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
% Entropy لصورة 8-بت (قناة واحدة أو 3 قنوات)
    if ndims(Iu8) == 3 && size(Iu8,3) == 3
        H = mean([image_entropy(Iu8(:,:,1)), image_entropy(Iu8(:,:,2)), image_entropy(Iu8(:,:,3))]);
        return;
    end
    counts = imhist(Iu8, 256);
    p = counts / sum(counts);
    p = p(p>0);
    H = -sum(p .* log2(p));
end

function [NPCR, UACI] = npcr_uaci(C1, C2)
% NPCR & UACI بين صورتين بنفس الأبعاد (2D أو 3D)
    C1 = uint8(C1); C2 = uint8(C2);
    assert(isequal(size(C1), size(C2)), 'npcr_uaci: size mismatch');
    N = numel(C1);
    D = C1 ~= C2;
    NPCR = 100 * sum(D(:)) / N;
    UACI = 100 * mean( abs(double(C1(:)) - double(C2(:))) / 255 );
end

function key = sha256_bytes(str)
% 32 بايت من SHA-256 لعبارة نصية
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
    burn = 1000; [x0,y0] = iter(x0, y0, a, b, burn);
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
        xn1 = frac(0.5*(xn1+1));  yn1 = frac(0.5*(yn1+1));
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
                t = bitxor(P(r,c), ks(r,c));                 % XOR
                s = uint16(t) + uint16(prevC) + uint16(prevP); % جمع نمطي (mod 256)
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
        prevC = iv; prevP = uint8(bitxor(iv, 171));
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
        prevC = iv; prevP = uint8(bitxor(iv, 219));
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
