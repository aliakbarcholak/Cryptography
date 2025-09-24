function results = MIE_Secure_Clean()
% MIE_Secure_Clean
% ------------------------------------------------------------
% Secure, modular implementation that upgrades a mirror-like image
% encryption (MIE) into MIE+ with:
% - Keyed permutation (row/col) + optional mirror flip
% - Dynamic per-image S-Box (+ inverse)
% - Multi-round diffusion (forward/backward) with independent keys/IVs
% - MAC (keyed SHA-256) for integrity checking
% - Full metrics: Entropy (DRPE-probe), NPCR, UACI, SSIM, PSNR, MSE
%
% ملاحظات:
% - mirror-like (flip) أصبح خطوة اختيارية داخل الـ Permutation (cfg.useMirror)
% - النتائج تُحفظ تلقائياً في Workspace باسم cfg.workspaceVarName (افتراضياً: 'MIE_results')
%
% Save as: MIE_Secure_Clean.m
% Call:    results = MIE_Secure_Clean;
% ------------------------------------------------------------

%% ===================== Configuration =====================
cfg.folderPath       = 'F:\DataBase\PNG';      % <-- عدّل المسار
cfg.passphrase       = 'Strong@Pass-Example!'; % <-- عدّل عبارة المرور
cfg.resizeTo         = [1024 1024];              % ثابت لتوحيد القياس (أو [] لتعطيل)
cfg.rounds           = 4;                      % 2..4 يرفع NPCR/UACI/Entropy
cfg.useDRPE          = true;                   % قياس Entropy(|E|) على |E|
cfg.showFigs         = true;                   % عرض الصور
cfg.writeCSV         = false;                  % حفظ CSV
cfg.csvName          = 'MIE_metrics.csv';
cfg.workspaceVarName = 'MIE_results';
cfg.useMirror        = true;                   % خطوة mirror داخل الـ permutation
cfg.patterns         = {'*.png','*.jpg','*.jpeg','*.tif','*.tiff','*.bmp', ...
                        '*.PNG','*.JPG','*.JPEG','*.TIF','*.TIFF','*.BMP'};

%% =================== Collect Image Files ==================
files = findImages(cfg.folderPath, cfg.patterns);
if isempty(files)
    error('No images found in "%s".', cfg.folderPath);
end

%% =============== Prepare Results Table (prealloc) =========
n = numel(files);
results = table('Size',[n 11], 'VariableTypes', ...
    {'string','double','double','double','double','double','double','double','double','double','string'}, ...
    'VariableNames', {'File','Entropy','NPCR','UACI','SSIM_enc','PSNR_Encryption','SSIM_dec','PSNR_dec','MSE_enc','MSE_dec','MAC'});

%% ====================== Main Loop =========================
for i = 1:n
    filePath = fullfile(cfg.folderPath, files(i).name);
    Iu8 = loadToUint8Gray(filePath, cfg.resizeTo);
    [H, W] = size(Iu8);

    % سياق المفاتيح لكل صورة
    ctx.passphrase = cfg.passphrase;
    ctx.fileId     = files(i).name;
    ctx.rounds     = cfg.rounds;
    ctx.useMirror  = cfg.useMirror;

    % ---- Encrypt (MIE+) ----
    [Ct, macHex] = miePlusEncrypt(Iu8, ctx);

    % ---- Decrypt + Verify MAC ----
    [Rec, macOk] = miePlusDecrypt(Ct, macHex, ctx); %#ok<NASGU>

    % ---- Entropy probe via DRPE ----
    if cfg.useDRPE
        keyDRPE = sha256Bytes([ctx.passphrase '|' ctx.fileId '|DRPE']);
        [s1, s2] = chaoticCoupledLogSin(keyDRPE, H*W);
        P1 = exp(1i*2*pi*reshape(mod(s1,1), H, W));
        P2 = exp(1i*2*pi*reshape(mod(s2,1), H, W));
        E  = fft2( fft2(im2double(Ct)).*P1 ).*P2;
        entropyVal = imageEntropy(uint8(round(255*mat2gray(abs(E)))));
    else
        entropyVal = NaN;
    end

    % ---- NPCR & UACI (flip single pixel in input) ----
    If = Iu8; If(1,1) = bitxor(If(1,1), uint8(1));
    Cf = miePlusEncrypt(If, ctx);  % MAC غير مطلوب هنا
    [npcrVal, uaciVal] = npcrUaci(Ct, Cf);

    % ---- Quality Metrics ----
    I_d  = im2double(Iu8);
    R_d  = im2double(Rec);
    Ct_d = im2double(Ct);

    psnrEnc = psnr(Ct_d, I_d);
    ssimEnc = ssim(Ct_d, I_d);
    ssimDec = ssim(R_d,  I_d);
    psnrDec = psnr(R_d,  I_d);
    mseEnc  = immse(I_d, Ct_d);
    mseDec  = immse(I_d, R_d);

    % ---- Print & Plot ----
    fprintf('Image %d/%d: %s\n', i, n, files(i).name);
    fprintf('  Entropy(|E|): %.4f | NPCR: %.2f%% | UACI: %.2f%%\n', entropyVal, npcrVal, uaciVal);
    fprintf(['  PSNR(enc): %.2f dB | SSIM(enc): %.4f | ' ...
             'PSNR(dec): %.2f dB | SSIM(dec): %.4f | ' ...
             'MSE(enc): %.6f | MSE(dec): %.6f | MAC: %s\n\n'], ...
             psnrEnc, ssimEnc, psnrDec, ssimDec, mseEnc, mseDec, tern(macOk,'OK','FAIL'));

    if cfg.showFigs
        showTriptych(Iu8, Ct, Rec, sprintf('MIE+: %s', files(i).name), tern(macOk,'MAC OK','MAC FAIL'));
    end

    % ---- Save Row ----
    results.File(i)            = string(files(i).name);
    results.Entropy(i)         = entropyVal;
    results.NPCR(i)            = npcrVal;
    results.UACI(i)            = uaciVal;
    results.SSIM_enc(i)        = ssimEnc;
    results.PSNR_Encryption(i) = psnrEnc;
    results.SSIM_dec(i)        = ssimDec;
    results.PSNR_dec(i)        = psnrDec;
    results.MSE_enc(i)         = mseEnc;
    results.MSE_dec(i)         = mseDec;
    results.MAC(i)             = string(tern(macOk,'OK','FAIL'));
end

%% ===================== Dump Summary =======================
disp(results);
if cfg.writeCSV
    writetable(results, fullfile(cfg.folderPath, cfg.csvName));
    fprintf('Saved CSV -> %s\n', fullfile(cfg.folderPath, cfg.csvName));
end

% --- Save to workspace automatically ---
try
    if ~isvarname(cfg.workspaceVarName)
        warning('Invalid workspaceVarName "%s". Using default "MIE_results".', cfg.workspaceVarName);
        cfg.workspaceVarName = 'MIE_results';
    end
    assignin('base', cfg.workspaceVarName, results);
    fprintf('Saved results to workspace variable: %s\n', cfg.workspaceVarName);
catch ME
    warning('Failed to assign results to base workspace: %s', ME.message);
end

end  % ===== End of main function =====

%% ==================== Helper: I/O =========================
function files = findImages(folderPath, patterns)
    files = [];
    for k = 1:numel(patterns)
        files = [files; dir(fullfile(folderPath, patterns{k}))]; %#ok<AGROW>
    end
end

function Iu8 = loadToUint8Gray(path, resizeTo)
    I = imread(path);
    if ~isempty(resizeTo)
        I = imresize(I, resizeTo);
    end
    if ndims(I)==3
        I = rgb2gray(I);
    end
    if isa(I,'uint16')
        Iu8 = uint8( round(double(I)/double(intmax('uint16'))*255) );
    else
        Iu8 = im2uint8(I);
    end
end

function showTriptych(A, B, C, figTitle, sub3Title)
    figure('Name', figTitle, 'NumberTitle', 'off');
    subplot(1,3,1); imshow(A); title('Original');
    subplot(1,3,2); imshow(B); title('Encrypted');
    subplot(1,3,3); imshow(C); title(['Decrypted (' sub3Title ')']);
    drawnow;
end

%% ==================== Helper: Ternary =====================
function out = tern(cond, a, b)
    if cond, out = a; else, out = b; end
end

%% ==================== MIE+ (Top-Level) ====================
function [Ct, macHex] = miePlusEncrypt(Iu8, ctx)
% Enc: Permute(+mirror optional) → S-Box → Multi-Round Diffusion → MAC
    [H, W]  = size(Iu8);
    rounds  = ctx.rounds;
    keyBytes = sha256Bytes([ctx.passphrase '|' ctx.fileId]);

    [ksF, ksB, ivF, ivB] = deriveRoundKeys(keyBytes, H, W, rounds);
    needLen = H + W + 512;
    [seq1, seq2] = chaoticCoupledLogSin(keyBytes, needLen);

    % مرآة اختيارية قبل/بعد الفرز ليبقى "طابع" MIE
    P = Iu8;
    if ctx.useMirror
        P = fliplr(P); % الخطوة المرآوية
    end
    % Permute (row/col) بمؤشرات فوضوية
    [~, rowIdx] = sort(seq1(1:H), 'ascend');
    [~, colIdx] = sort(seq2(1:W), 'ascend');
    P = P(rowIdx, :);
    P = P(:, colIdx);

    % S-Box ديناميكية
    [S, ~] = sboxFromSeq(seq2(H+1:H+256));
    Psub = uint8(S(double(P)+1) - 1);

    % Multi-round diffusion
    X = Psub;
    for r = 1:rounds
        X = diffusePass(X, ksF{r}, ivF{r}, 'forward');
        X = diffusePass(X, ksB{r}, ivB{r}, 'backward');
    end
    Ct = X;

    % MAC على Ct
    macHex = toHex( sha256Bytes([keyBytes(:).' uint8('|MAC|') Ct(:).']) );
end

function [Irec, ok] = miePlusDecrypt(Ct, macHex, ctx)
% Dec: Verify MAC → Undiffuse → Inverse S-Box → Inverse Perm (+mirror opt.)
    [H, W]  = size(Ct);
    rounds  = ctx.rounds;
    keyBytes = sha256Bytes([ctx.passphrase '|' ctx.fileId]);

    macCheck = toHex( sha256Bytes([keyBytes(:).' uint8('|MAC|') Ct(:).']) );
    ok = strcmp(macHex, macCheck);

    [ksF, ksB, ivF, ivB] = deriveRoundKeys(keyBytes, H, W, rounds);
    needLen = H + W + 512;
    [seq1, seq2] = chaoticCoupledLogSin(keyBytes, needLen);

    % مؤشرات عكسية آمنة
    [~, rowIdx] = sort(seq1(1:H), 'ascend'); invRow = zeros(1,H); invRow(rowIdx) = 1:H;
    [~, colIdx] = sort(seq2(1:W), 'ascend'); invCol = zeros(1,W); invCol(colIdx) = 1:W;

    [~, Sinv] = sboxFromSeq(seq2(H+1:H+256));

    % إزالة النشر متعدد الجولات
    X = Ct;
    for r = rounds:-1:1
        X = undiffusePass(X, ksB{r}, ivB{r}, 'backward');
        X = undiffusePass(X, ksF{r}, ivF{r}, 'forward');
    end
    PsubRec = X;
    Prec    = uint8(Sinv(double(PsubRec)+1) - 1);

    % عكس الترتيب الصف/العمود
    tmp  = Prec(invRow, :);
    Irec = tmp(:, invCol);

    % استرجاع المرآة إن كانت مفعّلة
    if ctx.useMirror
        Irec = fliplr(Irec);
    end
end

%% ==================== Crypto Primitives ===================
function [ksF, ksB, ivF, ivB] = deriveRoundKeys(keyBytes, H, W, rounds)
    ksF = cell(rounds,1);
    ksB = cell(rounds,1);
    ivF = cell(rounds,1);
    ivB = cell(rounds,1);
    N = H*W;
    for r = 1:rounds
        kf = sha256Bytes([keyBytes(:).' uint8('|F') uint8(sprintf('%d',r))]);
        kb = sha256Bytes([keyBytes(:).' uint8('|B') uint8(sprintf('%d',r))]);
        [sF, ~] = chaoticCoupledLogSin(kf, N+16);
        [sB, ~] = chaoticCoupledLogSin(kb, N+16);
        ksF{r} = reshape( uint8(floor(mod(sF(1:N)*1e14, 256))), H, W );
        ksB{r} = reshape( uint8(floor(mod(sB(1:N)*1e14, 256))), H, W );
        ivF{r} = uint8(bitxor(kf(1), kf(end)));
        ivB{r} = uint8(bitxor(kb(1), kb(end)));
    end
end

function [S, Sinv] = sboxFromSeq(seq256)
    seq256 = seq256(:);
    if numel(seq256) < 256
        seq256 = repmat(seq256, ceil(256/numel(seq256)), 1);
    end
    [~, ord] = sort(mod(seq256(1:256),1), 'ascend');
    S    = zeros(256,1,'uint16'); S(ord) = uint16(0:255)+1;
    Sinv = zeros(256,1,'uint16'); Sinv(S) = uint16(0:255)+1;
end

function out = diffusePass(P, ks, iv, direction)
    [H, W] = size(P);
    out = zeros(H, W, 'uint8');
    if strcmpi(direction, 'forward')
        prevC = iv;
        prevP = uint8(bitxor(iv, 171));
        for r = 1:H
            for c = 1:W
                t = bitxor(P(r,c), ks(r,c));
                s = uint16(t) + uint16(prevC) + uint16(prevP);
                out(r,c) = uint8(mod(s, 256));
                prevC = out(r,c);
                prevP = P(r,c);
            end
        end
    else
        prevC = iv;
        prevP = uint8(bitxor(iv, 219));
        for r = H:-1:1
            for c = W:-1:1
                t = bitxor(P(r,c), ks(r,c));
                s = uint16(t) + uint16(prevC) + uint16(prevP);
                out(r,c) = uint8(mod(s, 256));
                prevC = out(r,c);
                prevP = P(r,c);
            end
        end
    end
end

function P = undiffusePass(Cin, ks, iv, direction)
    [H, W] = size(Cin);
    P = zeros(H, W, 'uint8');
    if strcmpi(direction, 'forward')
        prevC = iv;
        prevP = uint8(bitxor(iv, 171));
        for r = 1:H
            for c = 1:W
                s = int16(Cin(r,c)) - int16(prevC) - int16(prevP);
                t = uint8(mod(s, 256));
                P(r,c) = bitxor(t, ks(r,c));
                prevC = Cin(r,c);
                prevP = P(r,c);
            end
        end
    else
        prevC = iv;
        prevP = uint8(bitxor(iv, 219));
        for r = H:-1:1
            for c = W:-1:1
                s = int16(Cin(r,c)) - int16(prevC) - int16(prevP);
                t = uint8(mod(s, 256));
                P(r,c) = bitxor(t, ks(r,c));
                prevC = Cin(r,c);
                prevP = P(r,c);
            end
        end
    end
end

%% ==================== Metrics & Utils =====================
function Hh = imageEntropy(Iu8)
    if ~isa(Iu8,'uint8'), Iu8 = uint8(Iu8); end
    counts = imhist(Iu8, 256);
    p = counts / sum(counts);
    p = p(p > 0);
    Hh = -sum(p .* log2(p));
end

function [NPCR, UACI] = npcrUaci(A, B)
    A = uint8(A); B = uint8(B);
    assert(isequal(size(A), size(B)), 'npcrUaci: size mismatch');
    N = numel(A);
    NPCR = 100 * sum(A(:) ~= B(:)) / N;
    UACI = 100 * mean( abs(double(A(:)) - double(B(:))) / 255 );
end

function key = sha256Bytes(data)
    md = java.security.MessageDigest.getInstance('SHA-256');
    if isstring(data) || ischar(data)
        bytes = uint8(char(data));
    else
        bytes = uint8(data);
    end
    md.update(bytes(:));
    key = typecast(md.digest, 'uint8');
end

function [x, y] = chaoticCoupledLogSin(keyBytes, L)
% Coupled Logistic–Sine (0..1) seeded from 256-bit key.
    if ~isscalar(L), L = prod(double(L)); end
    L = max(1, round(double(L)));

    u32 = typecast(uint8(keyBytes), 'uint32');
    s1 = double(bitxor(u32(1),u32(3))) / double(intmax('uint32')) + 0.11;
    s2 = double(bitxor(u32(2),u32(4))) / double(intmax('uint32')) + 0.17;

    a  = 3.8 + 0.19*frac(s1*7.1);
    b  = 3.8 + 0.19*frac(s2*9.3);

    x0 = frac(sum(double(u32(5:6)))/2^32 + 0.123456);
    y0 = frac(sum(double(u32(7:8)))/2^32 + 0.654321);
    [x0,y0] = iter(x0,y0,a,b,1000);

    x  = zeros(L,1);
    y  = zeros(L,1);
    xn = x0;
    yn = y0;

    for k = 1:L
        [xn,yn] = stepMap(xn,yn,a,b);
        xn = frac(xn + 0.5*frac(s1*xn + s2*yn));
        yn = frac(yn + 0.5*frac(s2*yn + s1*xn));
        x(k) = frac(xn);
        y(k) = frac(yn);
    end

    function [xo,yo] = iter(xi,yi,aa,bb,n)
        xo = xi; yo = yi;
        for t = 1:n
            [xo,yo] = stepMap(xo,yo,aa,bb);
            xo = frac(xo + 0.3*frac(s1*yo));
            yo = frac(yo + 0.3*frac(s2*xo));
        end
    end

    function [xn1,yn1] = stepMap(xn0,yn0,aa,bb)
        xn1 = sin(pi*( aa*xn0*(1-xn0) + 4*yn0*(1-yn0) ));
        yn1 = sin(pi*( bb*yn0*(1-yn0) + 4*xn1*(1-xn1) ));
        xn1 = frac(0.5*(xn1+1));
        yn1 = frac(0.5*(yn1+1));
    end

    function z = frac(v), z = v - floor(v); end
end

function hx = toHex(u8)
    hx = lower(reshape(dec2hex(uint8(u8),2).',1,[]));
end
