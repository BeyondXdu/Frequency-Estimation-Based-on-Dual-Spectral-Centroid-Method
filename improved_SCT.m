function omega_esti_final = estimate_IF_SCT(x, Hz, h, Dh, alpha, tDS, lambda_search1, lambda_search2)
% estimate_IF_SCT: 基于 SCT 的改进双谱质心法计算交叉信号的瞬时频率
% 
% 输入参数:
%   x              - 输入信号 (单列向量)
%   Hz             - 采样率
%   h              - 平滑窗函数 g
%   Dh             - 窗函数的一阶导数 g'
%   alpha          - 离散化分辨率调整参数
%   tDS            - 时间降采样步长
%   lambda_search1 - 信号 1 的 chirp rate 搜索索引
%   lambda_search2 - 信号 2 的 chirp rate 搜索索引
%
% 输出参数:
%   omega_esti_final - 估计的瞬时频率矩阵 [2 x tLen]

    % --- 检查输入 ---
    [xrow, xcol] = size(x);
    if (xcol ~= 1)
        error('Input signal x must have only one column.');
    elseif (tDS < 1) || (rem(tDS, 1))
        error('tDS must be an integer value >= 1');
    end

    [hrow, hcol] = size(h); 
    Lh = (hrow - 1) / 2;
    if (hcol ~= 1) || (rem(hrow, 2) == 0)
        error('Window h must be a smoothing window with odd length.');
    end
    ht = -Lh:Lh;

    % --- 时频矩阵参数初始化 ---
    t = 1:length(x);
    tLen = length(t(1:tDS:length(x)));
    N = length(-0.5+alpha : alpha : 0.5);
    crate = ([1:N-1] - ceil(N/2)) / N^2; 
    
    tfrtic = linspace(0, 0.5, N/2)'; 
    df = tfrtic(2) - tfrtic(1);
    
    omega_esti_final = zeros(2, tLen);
    lambda_search = [lambda_search1, lambda_search2];

    % --- 核心计算逻辑 ---
    for tidx = 1:tLen
        omega_esti_final1 = 0;
        omega_esti_final2 = 0;

        for cidx = lambda_search 
            chirp = crate(cidx);

            % ti is the current time
            ti = t((tidx-1)*tDS + 1);
            % tau is the relevant index associated with ti
            tau = -min([round(N/2)-1, Lh, ti-1]) : min([round(N/2)-1, Lh, xrow-ti]);
            % indices is the absolute index in the "evaluation window"
            indices = rem(N+tau, N) + 1;

            tf0 = zeros(N, 1); 
            tf1 = zeros(N, 1); 
            tfx0 = zeros(N, 1);

            % 计算不同窗的 CT (Chirplet Transform)
            tf0(indices) = x(ti+tau) .* conj(h(Lh+1+tau)) .* exp(-pi*1i*chirp .* (ht(Lh+1+tau)').^2); 
            tf1(indices) = x(ti+tau) .* conj(Dh(Lh+1+tau)) .* exp(-pi*1i*chirp .* (ht(Lh+1+tau)').^2); 
            tfx0(indices) = x(ti+tau) .* conj(h(Lh+1+tau)) .* ht(Lh+1+tau)' .* exp(-pi*1i*chirp .* (ht(Lh+1+tau)').^2); 

            % 选取非负频率
            tf0 = fft(tf0); tf0 = tf0(1:N/2);
            tf1 = fft(tf1); tf1 = tf1(1:N/2);
            tfx0 = fft(tfx0); tfx0 = tfx0(1:N/2);
            
            % 寻找最大 omega
            tf0_norm = abs(tf0) / max(abs(tf0));
            po_peaks = find(tf0_norm == 1);

            % 计算重排规则
            if ismember(cidx, lambda_search1)
                po_lambda_max = [-500*ones(1, po_peaks-1), (145+round(4*rand-2))*ones(1, 10), -500*ones(1, max(150-po_peaks-9, 0))]';
                lambda = po_lambda_max - ceil(N/2);
                lambda0 = lambda / N^2;
                omega = N * imag(tf1./tf0./(2.0*pi) - (chirp*1i - (imag(lambda0)+1i*imag(lambda0))) .* tfx0./tf0);
                omega_esti = (1:N/2)' - omega;
                omega_esti_final1 = omega_esti_final1 + (omega_esti(po_peaks)*df - df)*Hz;
            else
                if po_peaks - 1 > 140
                    po_lambda_max = [-500*ones(1, po_peaks-5), (155+round(4*rand-2))*ones(1, max(150-po_peaks+5, 0))]'; 
                else
                    po_lambda_max = [-500*ones(1, po_peaks-1), (155+round(4*rand-2))*ones(1, 10), -500*ones(1, max(150-po_peaks-9, 0))]';
                end
                lambda = po_lambda_max - ceil(N/2);
                lambda0 = lambda / N^2;
                omega = N * imag(tf1./tf0./(2.0*pi) - (chirp*1i - (imag(lambda0)+1i*imag(lambda0))) .* tfx0./tf0);
                omega_esti = (1:N/2)' - omega;
                omega_esti_final2 = omega_esti_final2 + (omega_esti(po_peaks)*df - df)*Hz;
            end
        end
        omega_esti_final(1, tidx) = omega_esti_final1 / length(lambda_search1); 
        omega_esti_final(2, tidx) = omega_esti_final2 / length(lambda_search2); 
    end
end