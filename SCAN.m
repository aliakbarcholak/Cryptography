function SCAN_Secure_Clean()
% SCAN_Secure_Clean
% ------------------------------------------------------------
% Secure SCAN-based image encryption with integrity (MAC) and
% full metrics: Entropy (via DRPE), NPCR, UACI, SSIM, PSNR
% (encryption & decryption), and MSE.
% Pipeline (encryption):
%  0) SCAN pre-mix (serpentine order reverse) – reversible
%  1) Keyed row/column permutation (chaotic, deterministic)
%  2) Dynamic invertible S-Box (256)
%  3) Multi-round diffusion (forward/backward) with per-round keys/IVs
%  4) MAC (keyed SHA-256) over ciphertext
% ------------------------------------------------------------

%% ===================== Configuration =====================
cfg.folderPath     = 'F:\رسالة ماجستير\البيانات المطلوبة\DataBase\png';  % <-- عدّل المسار
cfg.passphrase     = 'Strong@Pass-Example!';  % <-- عبارة المرور
cfg.resizeTo       = [1024 1024];             % ثابت لتوحيد القياس (أو [] لتعطيل)
cfg.rounds         = 3;                       % 2..4 (4 يرفع NPCR/UACI/Entropy)
cfg.useDRPE        = true;                    % لقياس Entropy(|E|) على |E|
cfg.showFigs       = true;                    % عرض صور لكل ملف
cfg.writeCSV       = false;                   % حفظ النتائج CSV بجانب البيانات
cfg.csvName        = 'SCAN_metrics.csv';

% === خيارات Workspace/حفظ ملف MAT ===
cfg.pushToWorkspace = true;                   % يدفع النتائج للـ Base Workspace
cfg.saveMat         = false;                  % يحفظ .mat
cfg.matName         = 'SCAN_results.mat';     % اسم ملف MAT عند الحفظ

%% =================== Collect Image Files ==================
files = findImages(cfg.folderPath, {'*.png'});
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

    % Per-image crypto context
    ctx.passphrase = cfg.passphrase;
    ctx.fileId     = files(i).name;
    ctx.rounds     = cfg.rounds;

    % ---- Encrypt (SCAN-secure) ----
    [Ct, macHex] = scanEncrypt(Iu8, ctx);

    % ---- Decrypt + Verify MAC ----
    [Rec, macOk] = scanDecrypt(Ct, macHex, ctx); %#ok<NASGU>

    % ---- Entropy(|E|) via DRPE ----
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

    % ---- NPCR & UACI (flip 1 input pixel) ----
    If = Iu8; If(1,1) = bitxor(If(1,1), uint8(1));
    Cf = scanEncrypt(If, ctx);  % encryption only (MAC not needed for NPCR/UACI)
    [npcrVal, uaciVal] = npcrUaci(Ct, Cf);

    % ---- Quality Metrics ----
    I_d  = im2double(Iu8);
    R_d  = im2double(Rec);
    Ct_d = im2double(Ct);

    psnrEnc = psnr(Ct_d, I_d);     % PSNR بين الأصل والمُعمّى
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
        showTriptych(Iu8, Ct, Rec, sprintf('SCAN-Secure: %s', files(i).name), tern(macOk,'MAC OK','MAC FAIL'));
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

% ========= Push to Workspace / Save .MAT =========
if cfg.pushToWorkspace
    try
        assignin('base','SCAN_results_table',results);
        assignin('base','SCAN_config',cfg);
        assignin('base','SCAN_files',string({files.name}).');
        fprintf('Pushed variables to Workspace: SCAN_results_table, SCAN_config, SCAN_files\n');
    catch ME
        warning('Could not push variables to Workspace: %s', ME.message);
    end
end
if cfg.saveMat
    try
        save(fullfile(cfg.folderPath, cfg.matName), 'results', 'cfg', 'files','-v7');
        fprintf('Saved MAT -> %s\n', fullfile(cfg.folderPath, cfg.matName));
    catch ME
        warning('Could not save MAT file: %s', ME.message);
    end
end

end
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

%% ==================== Crypto: Top-Level ===================
function [Ct, macHex] = scanEncrypt(Iu8, ctx)
% Enc: SCAN pre-mix → Permutation → S-Box → Multi-Round Diffusion → MAC
    [H, W] = size(Iu8);
    rounds = ctx.rounds;

    keyBytes = sha256Bytes([ctx.passphrase '|' ctx.fileId]);

    [ksF, ksB, ivF, ivB] = deriveRoundKeys(keyBytes, H, W, rounds);
    needLen = H + W + 512;
    [seq1, seq2] = chaoticCoupledLogSin(keyBytes, needLen);

    [~, rowIdx] = sort(seq1(1:H), 'ascend');
    [~, colIdx] = sort(seq2(1:W), 'ascend');

    [S, ~] = sboxFromSeq(seq2(H+1:H+256));

    % 0) SCAN pre-mix (serpentine order reverse) – reversible
    X = scanPremixSerpentine(Iu8);

    % 1) Permute (rows/cols)
    X = X(rowIdx, :);
    X = X(:, colIdx);

    % 2) Substitute (S-Box)
    X = uint8(S(double(X)+1) - 1);

    % 3) Multi-round diffusion
    for r = 1:rounds
        X = diffusePass(X, ksF{r}, ivF{r}, 'forward');
        X = diffusePass(X, ksB{r}, ivB{r}, 'backward');
    end
    Ct = X;

    % 4) MAC over Ct
    macHex = toHex( sha256Bytes([keyBytes(:).' uint8('|MAC|') Ct(:).']) );
end

function [Irec, ok] = scanDecrypt(Ct, macHex, ctx)
% Dec: Verify MAC → Undiffuse → Inverse S-Box → Inverse Perm → Inverse SCAN
    [H, W] = size(Ct);
    rounds  = ctx.rounds;
    keyBytes = sha256Bytes([ctx.passphrase '|' ctx.fileId]);

    macCheck = toHex( sha256Bytes([keyBytes(:).' uint8('|MAC|') Ct(:).']) );
    ok = strcmp(macHex, macCheck);

    [ksF, ksB, ivF, ivB] = deriveRoundKeys(keyBytes, H, W, rounds);
    needLen = H + W + 512;
    [seq1, seq2] = chaoticCoupledLogSin(keyBytes, needLen);

    [~, rowIdx] = sort(seq1(1:H), 'ascend'); invRow(rowIdx) = 1:H;
    [~, colIdx] = sort(seq2(1:W), 'ascend'); invCol(colIdx) = 1:W;
    [~, Sinv]   = sboxFromSeq(seq2(H+1:H+256));

    % عكس الانتشار
    X = Ct;
    for r = rounds:-1:1
        X = undiffusePass(X, ksB{r}, ivB{r}, 'backward');
        X = undiffusePass(X, ksF{r}, ivF{r}, 'forward');
    end

    % عكس الاستبدال
    X = uint8(Sinv(double(X)+1) - 1);

    % عكس الترتيب
    X = X(invRow, :);
    X = X(:, invCol);

    % عكس SCAN pre-mix (نفس العملية)
    Irec = scanPremixSerpentine(X);
end

%% ==================== Crypto: Primitives ==================
function X = scanPremixSerpentine(I)
% Serpentine (boustrophedon) 1D order with reverse assignment (reversible).
% Equivalent to: linearize by serpentine, reverse, and map back.
    [H,W] = size(I);
    idx = serpentineOrder(H, W);     % 1×N linear indices
    X = I;
    X(idx) = X(idx(end:-1:1));       % عكس على ترتيب المسح
end

function idx = serpentineOrder(H, W)
% Returns linear indices of a row-wise serpentine scan.
    idx = zeros(1, H*W, 'uint32');
    k = 1;
    for r = 1:H
        if mod(r,2)==1
            cols = 1:W;
        else
            cols = W:-1:1;
        end
        for c = cols
            idx(k) = sub2ind([H W], r, c);
            k = k + 1;
        end
    end
end

function [ksF, ksB, ivF, ivB] = deriveRoundKeys(keyBytes, H, W, rounds)
% Generate per-round forward/backward keystreams (H×W) and 1-byte IVs.
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
% Build 256-entry S-Box (0..255 -> 1..256) and its inverse from a chaotic seq.
    seq256 = seq256(:);
    if numel(seq256) < 256
        seq256 = repmat(seq256, ceil(256/numel(seq256)), 1);
    end
    [~, ord] = sort(mod(seq256(1:256),1), 'ascend');
    S    = zeros(256,1,'uint16'); S(ord) = uint16(0:255)+1;
    Sinv = zeros(256,1,'uint16'); Sinv(S) = uint16(0:255)+1;
end

function out = diffusePass(P, ks, iv, direction)
% Single pass diffusion with XOR, modular add, and feedback (CFB-like).
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
% Exact inverse of diffusePass.
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
% Returns 1×32 uint8 SHA-256 digest for char/string/uint8 input.
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
% Coupled Logistic–Sine (0..1) seeded from 256-bit key (robust to non-scalar L).
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

    x  = zeros(L,1); y = zeros(L,1);
    xn = x0;        yn = y0;

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
% Converts bytes to lowercase hex string.
    hx = lower(reshape(dec2hex(uint8(u8),2).',1,[]));
end
