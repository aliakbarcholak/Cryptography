clc; clear; close all;

%% ========= إعدادات =========
folder_path   = 'F:\DataBase\PNG';  % <-- عدّل المسار
passphrase    = 'Strong@Pass-Example!';   % <-- عدّل العبارة السرّية
resize_to     = [1024 1024];              % حجم موحّد
block_size    = 4;                        % يجب أن يقسّم الأبعاد
num_clusters  = 256;                      % <=256 (نطاق بايت)
rounds        = 3;                        % 2–4 مناسب (4 يرفع NPCR/UACI/Entropy)
use_drpe      = true;                     % طبقة DRPE لقياس Entropy

% جمع الصور
image_files = dir(fullfile(folder_path, '*.PNG'));
if isempty(image_files), error('No images found in %s', folder_path); end

%% ========= حاويات القياسات =========
n = numel(image_files);
MSE_encrypted = zeros(1,n);
MSE_decrypted = zeros(1,n);
PSNR_enc      = zeros(1,n);
PSNR_dec      = zeros(1,n);
SSIM_enc      = zeros(1,n);
SSIM_dec      = zeros(1,n);
NPCR_vals     = zeros(1,n);
UACI_vals     = zeros(1,n);
ENT_vals      = zeros(1,n);

for i = 1:n
    %% ======== قراءة وتحضير الصورة (رمادي) ========
    p = fullfile(folder_path, image_files(i).name);
    Iraw = imread(p);
    if ~isempty(resize_to), Iraw = imresize(Iraw, resize_to); end
    Igray = im2double(im2gray(Iraw));           % [0,1] double
    [R,C] = size(Igray);
    bs = block_size; K = num_clusters;
    assert(mod(R,bs)==0 && mod(C,bs)==0, 'Image must be divisible by block_size.');
    assert(K<=256, 'num_clusters must be <= 256.');

    %% ======== قاموس حتمي (KDF) وبناء فهارس VQ ========
    key_material = [passphrase '|' image_files(i).name];
    key_bytes    = sha256_bytes(key_material);   % 32 بايت

    % نحفظ حالة المولّد العشوائي ونثبت البداية بمفتاح مشتق حتى تكون النتيجة قابلة للتكرار
    rng_state = rng;
    rng(double(typecast(key_bytes(1:4),'uint32')),'twister'); %#ok<RNGR>
    [codebook, idx_mat] = vq_build_codebook(Igray, bs, K);  % idx_mat: 1..K
    rng(rng_state); % نعيد الحالة

    Rb = R/bs; Cb = C/bs;
    idx0 = uint8(idx_mat - 1);                % 0..K-1 على شبكة الكتل

    %% ======== مفاتيح فوضوية على شبكة الكتل ========
    need_len_blk = Rb*Cb*8 + Rb + max(256,K) + 64*rounds;
    [seq1_blk, seq2_blk] = chaotic_coupled_logsin(key_bytes, need_len_blk);

    % Permutation indices
    [~, row_idx] = sort(seq1_blk(1:Rb),'ascend');
    [~, col_idx] = sort(seq2_blk(1:Cb),'ascend');
    inv_row = zeros(1,Rb); inv_row(row_idx)=1:Rb;
    inv_col = zeros(1,Cb); inv_col(col_idx)=1:Cb;

    % Keystreams rounds + IVs
    [ksF, ksB, ivF, ivB] = derive_round_keys(key_bytes, Rb, Cb, rounds);

    % S-Box ديناميكية بطول K وعكسها
    [S, Sinv] = sbox_from_seq_len(seq2_blk(Rb+1:Rb+max(256,K)), K);

    %% ======== تشفير فهارس الكتل ========
    P    = idx0(row_idx, :);                 % تبديل صفوف
    P    = P(:, col_idx);                    % تبديل أعمدة
    Psub = uint8(S(double(P)+1)-1);          % استبدال (0..K-1)

    X = Psub;
    for r = 1:rounds
        X = diffuse_pass(X, ksF{r}, ivF{r}, 'forward');
        X = diffuse_pass(X, ksB{r}, ivB{r}, 'backward');
    end
    Cidx = X;                                % نص مُعمّى على شبكة الكتل (Rb×Cb)، uint8

    %% ======== صورة مُعمّاة قابلة للعرض (فك VQ فقط) ========
    Csp_px = vq_decode(codebook, double(Cidx)+1, bs, [R C]);   % double [0,1]

    %% ======== DRPE مستقل لقياس Entropy(|E|) ========
    if use_drpe
        key_drpe = sha256_bytes([key_material '|DRPE']);
        [s1d, s2d] = chaotic_coupled_logsin(key_drpe, R*C);
        P1 = exp(1i * 2*pi * reshape(mod(s1d,1), R, C));
        P2 = exp(1i * 2*pi * reshape(mod(s2d,1), R, C));
        E  = fft2(fft2(Csp_px).*P1).*P2;
        ENT_vals(i) = image_entropy(uint8(round(255*mat2gray(abs(E)))));
    else
        ENT_vals(i) = NaN;
    end

    %% ======== فك التشفير (عكس كامل الشبكة) ========
    X = Cidx;
    for r = rounds:-1:1
        X = undiffuse_pass(X, ksB{r}, ivB{r}, 'backward');
        X = undiffuse_pass(X, ksF{r}, ivF{r}, 'forward');
    end
    Psub_rec = X;
    Prec     = uint8(Sinv(double(Psub_rec)+1)-1);
    tmp      = Prec(inv_row, :);
    idx_rec0 = tmp(:, inv_col);              % 0..K-1
    idx_rec  = double(idx_rec0)+1;           % 1..K
    Irec     = vq_decode(codebook, idx_rec, bs, [R C]);  % double [0,1]

    %% ======== NPCR/UACI: قلب بكسل واحد بنفس القاموس والمفاتيح ========
    Iu8 = uint8(round(Igray*255)); Iu8f = Iu8; Iu8f(1,1)=bitxor(Iu8f(1,1),uint8(1));
    Iflip    = double(Iu8f)/255;
    idx_flip = vq_encode_with_codebook(Iflip, codebook, bs);  % 1..K
    Pf = uint8(idx_flip-1); Pf = Pf(row_idx,:); Pf = Pf(:,col_idx);
    Pf = uint8(S(double(Pf)+1)-1);
    Cf = Pf;
    for r = 1:rounds
        Cf = diffuse_pass(Cf, ksF{r}, ivF{r}, 'forward');
        Cf = diffuse_pass(Cf, ksB{r}, ivB{r}, 'backward');
    end
    % تكبير لخريطة الفهارس للحجم الكامل لتحقيق مقارنة بكسلية عادلة
    Cidx_full   = repelem(Cidx, bs, bs);
    Cidx_full_f = repelem(Cf,    bs, bs);
    [NPCR_vals(i), UACI_vals(i)] = npcr_uaci(Cidx_full, Cidx_full_f);

    %% ======== قياسات الجودة ========
    PSNR_enc(i) = psnr(Csp_px, Igray);
    SSIM_enc(i) = ssim(Csp_px, Igray);
    PSNR_dec(i) = psnr(Irec,   Igray);
    SSIM_dec(i) = ssim(Irec,   Igray);
    MSE_encrypted(i) = immse(Igray, Csp_px);
    MSE_decrypted(i) = immse(Igray, Irec);

    %% ======== طباعة وعرض ========
    fprintf('Image %d/%d: %s\n', i, n, image_files(i).name);
    fprintf('  Entropy(|E|): %.4f | NPCR: %.2f%% | UACI: %.2f%%\n', ENT_vals(i), NPCR_vals(i), UACI_vals(i));
    fprintf('  PSNR(enc): %.2f dB | SSIM(enc): %.4f | MSE(enc): %.6f\n', PSNR_enc(i), SSIM_enc(i), MSE_encrypted(i));
    fprintf('  PSNR(dec): %.2f dB | SSIM(dec): %.4f | MSE(dec): %.6f\n\n', PSNR_dec(i), SSIM_dec(i), MSE_decrypted(i));

    figure('Name', sprintf('VQ-Secure: %s', image_files(i).name), 'NumberTitle','off');
    subplot(1,3,1); imshow(Igray);  title('Original');
    subplot(1,3,2); imshow(Csp_px); title('Encrypted (spatial)');
    subplot(1,3,3); imshow(Irec);   title('Decrypted (VQ)');
    drawnow;
end

%% ======== ملخص النتائج ========
disp('--- Entropy(|E|) ---'); disp(ENT_vals);
disp('--- NPCR (%) ---');     disp(NPCR_vals);
disp('--- UACI (%) ---');     disp(UACI_vals);
disp('--- PSNR(enc) ---');    disp(PSNR_enc);
disp('--- PSNR(dec) ---');    disp(PSNR_dec);
disp('--- SSIM(enc) ---');    disp(SSIM_enc);
disp('--- SSIM(dec) ---');    disp(SSIM_dec);
disp('--- MSE Enc ---');      disp(MSE_encrypted);
disp('--- MSE Dec ---');      disp(MSE_decrypted);

%% ================== الدوال المساعدة ==================
function [codebook, idx_mat] = vq_build_codebook(Igray, bs, K)
% Build codebook with deterministic KMeans and encode to indices (1..K).
    [R,C] = size(Igray);
    blocks = im2col(Igray, [bs bs], 'distinct')';     % N×D
    opts = statset('MaxIter', 100, 'UseParallel', false, 'Display','off');
    [idx, cb] = kmeans(blocks, K, 'Options', opts, 'Replicates', 1, ...
                       'Start','plus', 'Distance','sqeuclidean');
    codebook = cb;                                     % K×D
    idx_mat  = reshape(idx, R/bs, C/bs);               % 1..K
end

function idx_mat = vq_encode_with_codebook(Igray, codebook, bs)
% Encode to indices (1..K) using given codebook.
    [R,C] = size(Igray);
    blocks = im2col(Igray, [bs bs], 'distinct');       % D×N
    D = size(codebook,2);
    assert(D==bs*bs, 'Codebook dimension mismatch.');
    x2 = sum(blocks.^2, 1);
    c2 = sum(codebook.^2, 2);
    M  = (-2)*(codebook*blocks) + c2 + x2;             % K×N
    [~, idx] = min(M, [], 1);
    idx_mat = reshape(idx, R/bs, C/bs);                % 1..K
end

function I = vq_decode(codebook, idx_mat, bs, orig_size)
% Decode indices to double image [0,1].
    blocks = codebook(idx_mat(:), :)';
    I = col2im(blocks, [bs bs], orig_size, 'distinct');
    I = min(max(I,0),1);
end

function H = image_entropy(Iu8)
% Entropy على صورة 8-بت
    if ~isa(Iu8,'uint8'), Iu8 = uint8(Iu8); end
    counts = imhist(Iu8, 256);
    p = counts / sum(counts);
    p = p(p>0);
    H = -sum(p .* log2(p));
end

function [NPCR, UACI] = npcr_uaci(A, B)
% NPCR & UACI لصورتين بنفس الأبعاد (uint8)
    A = uint8(A); B = uint8(B);
    assert(isequal(size(A), size(B)), 'npcr_uaci: size mismatch');
    N = numel(A);
    D = A ~= B;
    NPCR = 100 * sum(D(:)) / N;
    UACI = 100 * mean( abs(double(A(:)) - double(B(:))) / 255 );
end

function key = sha256_bytes(data)
% 32-byte SHA-256 (يقبل نصوصًا أو بايتات)
    md = java.security.MessageDigest.getInstance('SHA-256');
    if isstring(data) || ischar(data)
        bytes = uint8(char(data));
    else
        bytes = uint8(data);
    end
    md.update(bytes(:));                  % عمود بايتات
    key = typecast(md.digest, 'uint8');   % 1×32
end

function [x, y] = chaotic_coupled_logsin(key_bytes, L)
% Coupled Logistic–Sine chaotic generator seeded from 256-bit key.
    u32 = typecast(uint8(key_bytes), 'uint32');
    s1 = double(bitxor(u32(1), u32(3)))/double(intmax('uint32')) + 0.11;
    s2 = double(bitxor(u32(2), u32(4)))/double(intmax('uint32')) + 0.17;
    a  = 3.8 + 0.19*frac(s1*7.1);
    b  = 3.8 + 0.19*frac(s2*9.3);
    x0 = frac(sum(double(u32(5:6)))/2^32 + 0.123456);
    y0 = frac(sum(double(u32(7:8)))/2^32 + 0.654321);
    [x0,y0] = iter(x0,y0,a,b,1000);
    x=zeros(L,1); y=zeros(L,1); xn=x0; yn=y0;
    for k=1:L
        [xn,yn] = step_map(xn,yn,a,b);
        xn = frac(xn + 0.5*frac(s1*xn + s2*yn));
        yn = frac(yn + 0.5*frac(s2*yn + s1*xn));
        x(k)=frac(xn); y(k)=frac(yn);
    end
    function [xo,yo]=iter(xi,yi,aa,bb,n)
        xo=xi; yo=yi;
        for t=1:n
            [xo,yo]=step_map(xo,yo,aa,bb);
            xo = frac(xo + 0.3*frac(s1*yo));
            yo = frac(yo + 0.3*frac(s2*xo));
        end
    end
    function [xn1,yn1]=step_map(xn,yn,aa,bb)
        xn1 = sin(pi*( aa*xn*(1-xn) + 4*yn*(1-yn) ));
        yn1 = sin(pi*( bb*yn*(1-yn) + 4*xn1*(1-xn1) ));
        xn1 = frac(0.5*(xn1+1));
        yn1 = frac(0.5*(yn1+1));
    end
    function z=frac(v), z=v-floor(v); end
end

function out = diffuse_pass(P, ks, iv, dir)
% انتشار بتمرير واحد (XOR + جمع نمطي + تغذية راجعة)
    [R,CC] = size(P); out = zeros(R,CC,'uint8');
    if strcmpi(dir,'forward')
        prevC=iv; prevP=uint8(bitxor(iv,171));
        for r=1:R, for c=1:CC
            t = bitxor(P(r,c), ks(r,c));
            s = uint16(t)+uint16(prevC)+uint16(prevP);
            out(r,c) = uint8(mod(s,256));
            prevC = out(r,c); prevP = P(r,c);
        end, end
    else
        prevC=iv; prevP=uint8(bitxor(iv,219));
        for r=R:-1:1, for c=CC:-1:1
            t = bitxor(P(r,c), ks(r,c));
            s = uint16(t)+uint16(prevC)+uint16(prevP);
            out(r,c) = uint8(mod(s,256));
            prevC = out(r,c); prevP = P(r,c);
        end, end
    end
end

function P = undiffuse_pass(Cin, ks, iv, dir)
% العكس الدقيق لـ diffuse_pass
    [R,CC] = size(Cin); P = zeros(R,CC,'uint8');
    if strcmpi(dir,'forward')
        prevC=iv; prevP=uint8(bitxor(iv,171));
        for r=1:R, for c=1:CC
            s = int16(Cin(r,c)) - int16(prevC) - int16(prevP);
            t = uint8(mod(s,256));
            P(r,c) = bitxor(t, ks(r,c));
            prevC = Cin(r,c); prevP = P(r,c);
        end, end
    else
        prevC=iv; prevP=uint8(bitxor(iv,219));
        for r=R:-1:1, for c=CC:-1:1
            s = int16(Cin(r,c)) - int16(prevC) - int16(prevP);
            t = uint8(mod(s,256));
            P(r,c) = bitxor(t, ks(r,c));
            prevC = Cin(r,c); prevP = P(r,c);
        end, end
    end
end

function [S, Sinv] = sbox_from_seq_len(seq_in, M)
% S-Box بطول M (قيم 0..M-1 تُفهرس 1..M) وعكسها
    seq_in = seq_in(:);
    if numel(seq_in) < M
        seq_in = repmat(seq_in, ceil(M/numel(seq_in)), 1);
    end
    [~, ord] = sort(mod(seq_in(1:M),1), 'ascend');
    S = zeros(M,1,'uint16');  S(ord) = uint16(0:M-1)+1;
    Sinv = zeros(M,1,'uint16'); Sinv(S) = uint16(0:M-1)+1;
end

function [ksF, ksB, ivF, ivB] = derive_round_keys(key_bytes, Rb, Cb, rounds)
% توليد مفاتيح/IVs مستقلة لكل جولة على شبكة الكتل (Rb×Cb)
    ksF = cell(rounds,1); ksB = cell(rounds,1);
    ivF = cell(rounds,1); ivB = cell(rounds,1);
    N   = Rb*Cb;
    for r = 1:rounds
        % اجعل كل شيء صفاً من uint8 قبل التجميع (إصلاح أي لبس أنواع)
        kf = sha256_bytes([key_bytes(:).' uint8('|F') uint8(sprintf('%d', r))]);
        kb = sha256_bytes([key_bytes(:).' uint8('|B') uint8(sprintf('%d', r))]);
        [sF, ~] = chaotic_coupled_logsin(kf, N + 16);
        [sB, ~] = chaotic_coupled_logsin(kb, N + 16);
        ksF{r} = uint8(reshape(uint8(floor(mod(sF(1:N)*1e14,256))), Rb, Cb));
        ksB{r} = uint8(reshape(uint8(floor(mod(sB(1:N)*1e14,256))), Rb, Cb));
        ivF{r} = uint8(bitxor(kf(1), kf(end)));
        ivB{r} = uint8(bitxor(kb(1), kb(end)));
    end
end
