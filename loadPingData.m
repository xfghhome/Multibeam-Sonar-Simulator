%% ============================ 读取函数 ============================
function [data, info] = loadPingData(fileName, outClass)
% loadPingData  读取 savePingData 写出的 Ping 数据文件
%
%   data = loadPingData(fileName)
%   [data, info] = loadPingData(fileName)
%   [...] = loadPingData(fileName, outClass)
%
% outClass 选项
%   'native' (默认)   -> 返回保存时使用的整数容器类型
%   'single' or 'double' -> 转为浮点并自动归一化到 [-1,1]
%
% 注意：若文件为 24 bit（存 int32），native 模式仍返回 int32，
%       有效位为低 24 bit；如需截断请手动位运算。

    if nargin < 2 || isempty(outClass), outClass = 'native'; end

    fid = fopen(fileName,'r','ieee-le');
    if fid<0, error('Cannot open %s for reading.', fileName); end
    cleaner = onCleanup(@() fclose(fid));

    % ---- 读取头 ----
    dtypeStr = char(deblank(char(fread(fid,8,'*uint8'))')); % 去尾零
    M = fread(fid,1,'int32');
    N = fread(fid,1,'int32');
    if isempty(dtypeStr); error('Invalid dtype string in header.'); end

    bitDepth = sscanf(dtypeStr,'int%d');
    if isempty(bitDepth)
        error('Unsupported dtype "%s" in file.', dtypeStr);
    end

    [container, bytesPerSample] = containerType(bitDepth);

    % ---- 读取数据 ----
    rawData = fread(fid, [N, M], ['*' container]);

    % ---- 类型转换 ----
    switch lower(outClass)
        case 'native'
            data = rawData;
        case {'single','double'}
            scale = double(2^(bitDepth-1)-1);
            data = cast(rawData, outClass) / scale;
        otherwise
            error('Unsupported outClass "%s".', outClass);
    end

    % ---- 信息结构体 ----
    info = struct(...
        'dtype',        dtypeStr, ...
        'bitDepth',     bitDepth, ...
        'container',    container, ...
        'cols',         double(M), ...
        'rows',         double(N), ...
        'headerSize',   16, ...
        'totalBytes',   16 + numel(rawData)*bytesPerSample, ...
        'scale',        (bitDepth<=64) * (2^(bitDepth-1)-1));
end

%% ============================ 辅助函数 ============================
function [ctype, bytes] = containerType(bitDepth)
% 根据有效位数选择整数容器类型
    if bitDepth<=8
        ctype = 'int8';  bytes = 1;
    elseif bitDepth<=16
        ctype = 'int16'; bytes = 2;
    elseif bitDepth<=32
        ctype = 'int32'; bytes = 4;
    else
        ctype = 'int64'; bytes = 8;
    end
end