function savePingData(dataMics, fileName, bitDepth)
% savePingData  归一化并带文件头保存多通道 Ping 数据（可选位宽）
%
%   savePingData(dataMics)                  % 默认保存为 int16，文件 ./ping.dat
%   savePingData(dataMics,'foo.dat')        % 同上，改文件名
%   savePingData(dataMics,'foo.dat',24)     % 指定位宽 24 bit（使用 int32 容器写入）
%
% 文件头格式（小端）
% -------------------------------------------------------------
% Offset | Size | 内容
%   0    |  8   | 数据类型字符串，如 'int16', 'int24', 'int32', 'int8'（不足补 0x00）
%   8    |  4   | 横向元素数 = 列数 M (int32)
%  12    |  4   | 纵向元素数 = 行数 N (int32)
%  16    | ...  | 按列展开的采样数据（与类型字符串一致）
% -------------------------------------------------------------
% * 支持位宽：8, 16, 24, 32, 64（>32 位通常无必要）。
% * bitDepth 非 8/16/32/64 时将存储到下一档容器（例：24→int32），
%   但头信息仍记录原始 bitDepth ('int24') 方便解析端截取有效位。

    if nargin < 2 || isempty(fileName)
        fileName = './ping.dat';
    end
    if nargin < 3 || isempty(bitDepth)
        bitDepth = 16;  % 默认 16 bit
    end

    % ---------- 参数合法性检查 ----------
    validateattributes(dataMics, {'numeric'}, {'nonempty'}, mfilename, 'dataMics');
    validateattributes(bitDepth, {'numeric'}, {'scalar', 'integer', '>=', 1, '<=', 64}, mfilename, 'bitDepth');

    % ---------- 路径检查 ----------
    [pathStr, ~, ~] = fileparts(fileName);
    if ~isempty(pathStr) && ~exist(pathStr, 'dir')
        mkdir(pathStr);
    end

    % ---------- 步骤 1：归一化 ----------
    dataMics = single(dataMics);           % 统一 precision
    maxAbs   = max(abs(dataMics(:)));
    if maxAbs == 0
        warning('Input data all zeros; output file will be all zeros.');
        normData = dataMics;
    else
        normData = dataMics / maxAbs;      % [-1, 1]
    end

    % ---------- 步骤 2：映射到整数 ----------
    % 选择合适的整数容器类型
    if bitDepth <= 8
        container = 'int8';   bytesPerSample = 1;  effectiveBits = 8;
    elseif bitDepth <= 16
        container = 'int16';  bytesPerSample = 2;  effectiveBits = 16;
    elseif bitDepth <= 32
        container = 'int32';  bytesPerSample = 4;  effectiveBits = 32;
    else
        container = 'int64';  bytesPerSample = 8;  effectiveBits = 64;
    end

    scale     = single(2^(bitDepth-1) - 1);  % 用原始 bitDepth 计算动态范围
    intData   = cast(round(normData * scale), container);

    % ---------- 步骤 3：写文件含头 ----------
    fid = fopen(fileName, 'w', 'ieee-le');  % 强制小端
    if fid < 0, error('Cannot open %s for writing.', fileName); end

    % --- 3-A 写文件头 ---
    dtypeStr   = sprintf('int%d', bitDepth); % 如 'int16', 'int24'
    fixedLen   = 8;                          % 头部固定 8 字节保存 dtype 字符串
    if numel(dtypeStr) > fixedLen
        error('dtypeStr "%s" longer than %d bytes.', dtypeStr, fixedLen);
    end
    dtypeBytes = [uint8(dtypeStr), zeros(1, fixedLen - numel(dtypeStr), 'uint8')];
    fwrite(fid, dtypeBytes, 'uint8');

    M = int32(size(dataMics, 2));           % 横向 (列/通道)
    N = int32(size(dataMics, 1));           % 纵向 (行/采样点)
    fwrite(fid, M, 'int32');
    fwrite(fid, N, 'int32');
    % 文件头到此共 16 字节

    % --- 3-B 写数据正文 ---
    count = fwrite(fid, intData, container);
    fclose(fid);

    fprintf(['Saved %d samples × %d channels to %s\n' ...
             'Header: [dtype=%s, horiz=%d, vert=%d]  Total %d bytes (含头).\n'], ...
            N, M, fileName, dtypeStr, M, N, 16 + count * bytesPerSample);
end
