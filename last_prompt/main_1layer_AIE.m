% main_1layer_AIE
%   캐비티 + 1층 임피던스 strip 메타표면 집속 설계 (AIE + 유한 플랜지)
%
%   구조 (main_1layer 와 같음):  z = 0 metasurface,  기판 -t_sub <= z <= 0,
%                               바닥 z = -t_sub - d,  피드 z = -t_sub - d/2
%   (1) 캐비티 안 (z < z_ap): x 방향 PEC 판 (간격 Lx_3D) 사이 -> 2D MoM (유전체, 바닥, 벽, strip)
%   (2) 개구면 z = z_ap: PEC 로 막고 그 위 접선 E (펄스, 폭 D_ap) 를 미지수로 추가
%         안쪽: 막힌 면의 image 를 모든 MoM 블록에 더함
%         바깥: 개구면 연산자 = 무한 플랜지 연산자 + (유한 - 무한) 보정
%         개구면에서 H_y 연속
%   (3) 바깥 유한 플랜지: 단면 (x-z) 금속 몸체를 ky 성분마다 2.5D MoM 으로 풀어
%       개구면 되먹임 Y(ky) 와 초점 field F(ky) 를 구함. 결과는 ffcache_*.mat 로 재사용
%       (개구면 양 끝 y = +-W/2 의 모서리는 무시)
%   (4) 최적화: strip Nm 개의 커패시턴스를 adjoint gradient + fmincon 으로 최적화
%       목적함수 = PSF 모양 + 세기 중심 + 가장자리 벌점 (w_edge) + 손실 벌점 (w_loss)
clc; clear; close all;
format short e;
format compact;

FS = 18;
Fonts = "Times New Roman";
LW = 2;

%% =========================================================================
% 1. 파라미터 및 목표 PSF
% =========================================================================
raw = readmatrix("SO_PSF\coeff_1D_MW1.20_FOV15_SLL-10.txt");

NA_pos = raw(:,2);
mag_pos = raw(:,3) / max(raw(:,3));
phase_pos = raw(:,4);

nu0      = 10e9;
omega0   = 2 * pi * nu0;
epsilon0 = 8.8541878176e-12;
mu0      = (4 * pi) * 1e-7;
c0       = 1 / sqrt(epsilon0 * mu0);
freq0    = 10e9;
lambda0  = c0 / freq0;
eta0     = sqrt(mu0 / epsilon0);
k0       = 2 * pi / lambda0;

Width = 8.000 * lambda0;
z_f   = Width * 0.7317;
y_ROI = 20 * lambda0;

optimize_R = true;
delta_y    = 0.00 * lambda0;
N_obs      = 20001;

% 메타표면 strip
Nm  = 24;                        % strip 개수
Dy  = Width / Nm;
dyM = lambda0 / 8;

% 기판 / 캐비티 (z = 0: metasurface, 기판 -t_sub <= z <= 0)
t_sub    = 1.52e-3;
eps_sub  = 2.20;
d        = 0.750 * lambda0;
z_ground = -t_sub - d;
z_feed   = -t_sub - d/2;

% AIE 모델 / 실행 설정
Lx_3D         = lambda0 / 3;      % x 방향 PEC 판 간격 = strip 길이 (HFSS dxM)
z_ap          = lambda0 / 40;     % 개구면 높이 = 벽 / x 판 윗끝
D_ap          = lambda0 / 40;     % 개구면 펄스 폭
replay_X_file = "";               % "" = 최적화 실행,  파일 경로 = 그 X 로 결과만 계산
save_folder   = "1layer_aie_results";
w_edge        = 0.0;              % 가장자리 field 벌점 가중치
w_loss        = 0.1;              % 손실 벌점 가중치 (P_ohmic / P_in)
edge_width    = 1.0 * lambda0;    % 벌점 구간: 개구면에서 |y| > W/2 - edge_width

% 유한 플랜지 설정
t_flange      = 15e-3;    % x 방향 두꺼운 벽 두께 = 개구면 앞뒤 플랜지 폭 (HFSS 와 같게)
ff_dl         = lambda0 / 36;     % 단면 금속 둘레 분할 간격
ff_nslot      = 12;               % 슬롯 (폭 Lx) 분할 수
ff_beta_Y     = 4.0;              % 개구면 연산자 보정에 쓰는 ky 최대 [k0 단위]
ff_beta_F     = 1.5;              % 초점면에 쓰는 ky 최대 [k0 단위]
ff_cache      = "";               % "" = 형상으로 이름을 만들어 자동 캐시

H_baffle = z_ap;                  % MoM 벽 높이 = 개구면
H        = d + t_sub + H_baffle;
assert(z_ap > 0, 'z_ap 는 Layer 1 (z = 0) 보다 위여야 합니다.');
Na_ap = round(Width / D_ap);
y_ap  = -Width/2 + ((1:Na_ap) - 0.5) * D_ap;       % 펄스 중심 (행 벡터)
gam_E = exp(0.5772156649);
fprintf('AIE model: Lx = %.3f mm, z_ap = %.3f mm, %d aperture pulses\n', Lx_3D*1e3, z_ap*1e3, Na_ap);

% 목표 PSF: 스펙트럼 (NA) -> IFFT -> 초점면 |E|
NA_num = ( (delta_y + Width/2)/sqrt((delta_y + Width/2)^2 + z_f^2) - (delta_y - Width/2)/sqrt((delta_y - Width/2)^2 + z_f^2) ) / 2;
NA = linspace(-NA_num, NA_num, 11999);
kx = NA * k0;

mag   = [flip(mag_pos(2:end)); mag_pos];
phase = [flip(phase_pos(2:end)); phase_pos];
phase = phase - kx(:) * delta_y * (180 / pi);

N_fine   = 20001;
dy_fine  = lambda0 / 20;
y_coords = linspace(-N_fine*dy_fine/2, N_fine*dy_fine/2, N_fine);
fx       = linspace(-1/(2*dy_fine), 1/(2*dy_fine), N_fine);
kx_prop  = 2 * pi * fx;

mag_interp   = interp1(kx, mag, kx_prop, 'linear', 0);
phase_interp = interp1(kx, phase, kx_prop, 'linear', 0);
A_interp     = mag_interp .* exp(1j * phase_interp * pi / 180);

E_focus_ideal     = fftshift(ifft(ifftshift(A_interp))) * N_fine;
E_focus_ideal_abs = abs(E_focus_ideal);

x_obs = linspace(-y_ROI + delta_y, y_ROI + delta_y, N_obs);
target_PSF_lin = interp1(y_coords, E_focus_ideal_abs, x_obs, 'pchip', 0);
target_PSF_lin = target_PSF_lin / max(target_PSF_lin);
target_PSF_lin = target_PSF_lin(:);

%% =========================================================================
% 2. 유닛셀 커패시턴스 -> (R, X) 테이블
% =========================================================================
x_orig     = [0.051, 0.078, 0.234];
y_res_orig = [2.40, 5.56, 34.3];
y_imp_orig = [-175.35, -139.5, +34.2];
n_fit      = 2;

[y_res_fit, y_imp_fit, x_fit] = resistance(x_orig, y_res_orig, y_imp_orig, n_fit);

%% =========================================================================
% 3. 캐비티 / 기판 / 1층 메타표면 이산화
% =========================================================================
yLW_loc = -Width / 2;
yRW_loc =  Width / 2;

Ng = 81; Nv = 81; Nw = 81;
Ngk = 10; Nvky = 10; Nvkz = 1; Nwk = 10;
Nk  = 10;

DyV  = Width / Nv;
DyG  = Width / Ng;
DzW  = H / Nw;
dyMK = dyM / Nk;
dyGK = DyG / Ngk;
dyVK = DyV / Nvky;
dzVK = t_sub / Nvkz;
dzWK = DzW / Nwk;

idx_vecM   = 0:Nm-1;
idx_vecV   = 0:Nv-1;
idx_vecG   = 0:Ng-1;
idx_vecW   = 0:Nw-1;
idx_vecMK  = 0:Nk-1;
idx_vecGK  = 0:Ngk-1;
idx_vecVKy = 0:Nvky-1;
idx_vecVKz = 0:Nvkz-1;
idx_vecWK  = 0:Nwk-1;

ymM_vec = -Width/2 + (idx_vecM + 1/2) * Dy;
ym1_vec = ymM_vec;

ymG_vec = -Width/2 + (idx_vecG + 1/2) * DyG;
ymV_vec = -Width/2 + (idx_vecV + 1/2) * DyV;
zmW_vec = z_ground + (idx_vecW + 1/2) * DzW;

ymMK_vec = zeros(1, Nm * Nk);
ymGK_vec = zeros(1, Ng * Ngk);
ymVK_vec = zeros(1, Nv * Nvky);
zmWK_vec = zeros(1, Nw * Nwk);

for m = 1:Nm
    ymMK_vec((m-1)*Nk+1 : m*Nk) = ymM_vec(m) - dyM/2 + (idx_vecMK + 1/2) * dyMK;
end
for m = 1:Ng
    ymGK_vec((m-1)*Ngk+1 : m*Ngk) = ymG_vec(m) - DyG/2 + (idx_vecGK + 1/2) * dyGK;
end
for m = 1:Nv
    ymVK_vec((m-1)*Nvky+1 : m*Nvky) = ymV_vec(m) - DyV/2 + (idx_vecVKy + 1/2) * dyVK;
end
for m = 1:Nw
    zmWK_vec((m-1)*Nwk+1 : m*Nwk) = zmW_vec(m) - DzW/2 + (idx_vecWK + 1/2) * dzWK;
end
zmV4K_vec = -t_sub + (idx_vecVKz + 1/2) * dzVK;

%% =========================================================================
% 4. 피드 (x 방향 선전류) 입사파
% =========================================================================
I0 = 1.0;
yi = 0;
const_inc = -1 * I0 * eta0 * k0 / 4;

% 기준 전력: 피드 주위 반지름 d/10 원을 지나는 전력
Exinc = const_inc * besselh(0, 2, k0 * (d/10));
Hyinc = (1i * I0 * k0 / 4) * besselh(1, 2, k0 * (d/10));
Szinc = 0.5 * real(Exinc * conj(Hyinc));
Winc  = Szinc * (2 * pi * d / 10);

Vm1 = const_inc * besselh(0, 2, k0 * sqrt((ymMK_vec - yi).^2 + (0 - z_feed)^2));

Vm4 = zeros(Nv * Nvky * Nvkz, 1);
for nz = 1:Nvkz
    idx4 = (nz-1) * (Nv * Nvky) + (1:Nv * Nvky);
    Vm4(idx4) = const_inc * besselh(0, 2, k0 * sqrt((ymVK_vec - yi).^2 + (zmV4K_vec(nz) - z_feed).^2));
end
VmG  = const_inc * besselh(0, 2, k0 * sqrt((ymGK_vec - yi).^2 + (z_ground - z_feed)^2));
VmLW = const_inc * besselh(0, 2, k0 * sqrt((yLW_loc - yi)^2 + (zmWK_vec - z_feed).^2));
VmRW = const_inc * besselh(0, 2, k0 * sqrt((yRW_loc - yi)^2 + (zmWK_vec - z_feed).^2));

Vi = Vm1(:);                                    % 능동 strip
Vs = [Vm4(:); VmG(:); VmLW(:); VmRW(:)];        % 고정 구조체

%% =========================================================================
% 5. 캐비티 MoM / VSIE 행렬 (GPU)
% =========================================================================
fprintf('--- 캐비티 및 1-Layer 스트립 산란 행렬 구축 시작 (Nvkz = %d) ---\n', Nvkz);

try
    gpuDev = gpuDevice;
catch ME
    error(['GPU를 초기화할 수 없습니다. Parallel Computing Toolbox와 ', ...
           '지원되는 NVIDIA GPU가 필요합니다.\n원래 오류: %s'], ME.message);
end
fprintf('GPU precompute: %s, free VRAM %.2f GB\n', gpuDev.Name, gpuDev.AvailableMemory/1024^3);

tPrecompute = tic;

% (1) 고정 구조체 간 상호작용
Z4_4_cell  = cell(Nvkz, Nvkz);
Z4_G_cell  = cell(Nvkz, 1);
Z4_LW_cell = cell(Nvkz, 1);
Z4_RW_cell = cell(Nvkz, 1);
ZG_4_cell  = cell(1, Nvkz);
ZLW_4_cell = cell(1, Nvkz);
ZRW_4_cell = cell(1, Nvkz);

for mz = 1:Nvkz
    for nz = 1:Nvkz
        Z4_4_cell{mz, nz} = GetScatteringMatrix_v2_GPU(nu0, Nv, Nv, Nvky, Nvky, ymV_vec, ymV_vec, ...
            DyV, dyVK, DyV, dyVK, dzVK, zmV4K_vec(mz), zmV4K_vec(nz), 'Dielectric', 'Dielectric');
    end
    Z4_G_cell{mz, 1}  = GetScatteringMatrix_v2_GPU(nu0, Nv, Ng, Nvky, Ngk, ymV_vec, ymG_vec, ...
        DyV, dyVK, DyG, dyGK, 0, zmV4K_vec(mz), z_ground, 'Dielectric', 'ImpSurf');
    Z4_LW_cell{mz, 1} = GetScatteringMatrix_v2_GPU(nu0, Nv, Nw, Nvky, Nwk, ymV_vec, yLW_loc, ...
        DyV, dyVK, DzW, dzWK, 0, zmV4K_vec(mz), zmW_vec, 'Dielectric', 'GWall');
    Z4_RW_cell{mz, 1} = GetScatteringMatrix_v2_GPU(nu0, Nv, Nw, Nvky, Nwk, ymV_vec, yRW_loc, ...
        DyV, dyVK, DzW, dzWK, 0, zmV4K_vec(mz), zmW_vec, 'Dielectric', 'GWall');
end

for nz = 1:Nvkz
    ZG_4_cell{1, nz}  = GetScatteringMatrix_v2_GPU(nu0, Ng, Nv, Ngk, Nvky, ymG_vec, ymV_vec, ...
        DyG, dyGK, DyV, dyVK, dzVK, z_ground, zmV4K_vec(nz), 'ImpSurf', 'Dielectric');
    ZLW_4_cell{1, nz} = GetScatteringMatrix_v2_GPU(nu0, Nw, Nv, Nwk, Nvky, yLW_loc, ymV_vec, ...
        DzW, dzWK, DyV, dyVK, dzVK, zmW_vec, zmV4K_vec(nz), 'GWall', 'Dielectric');
    ZRW_4_cell{1, nz} = GetScatteringMatrix_v2_GPU(nu0, Nw, Nv, Nwk, Nvky, yRW_loc, ymV_vec, ...
        DzW, dzWK, DyV, dyVK, dzVK, zmW_vec, zmV4K_vec(nz), 'GWall', 'Dielectric');
end

Z4_4  = Cell2MatGPU(Z4_4_cell);
Z4_G  = Cell2MatGPU(Z4_G_cell);
Z4_LW = Cell2MatGPU(Z4_LW_cell);
Z4_RW = Cell2MatGPU(Z4_RW_cell);
ZG_4  = Cell2MatGPU(ZG_4_cell);
ZLW_4 = Cell2MatGPU(ZLW_4_cell);
ZRW_4 = Cell2MatGPU(ZRW_4_cell);

ZG_G   = GetScatteringMatrix_v2_GPU(nu0, Ng, Ng, Ngk, Ngk, ymG_vec, ymG_vec, DyG, dyGK, DyG, dyGK, 0, z_ground, z_ground, 'ImpSurf', 'ImpSurf');
ZG_LW  = GetScatteringMatrix_v2_GPU(nu0, Ng, Nw, Ngk, Nwk, ymG_vec, yLW_loc, DyG, dyGK, DzW, dzWK, 0, z_ground, zmW_vec, 'ImpSurf', 'GWall');
ZG_RW  = GetScatteringMatrix_v2_GPU(nu0, Ng, Nw, Ngk, Nwk, ymG_vec, yRW_loc, DyG, dyGK, DzW, dzWK, 0, z_ground, zmW_vec, 'ImpSurf', 'GWall');

ZLW_G  = GetScatteringMatrix_v2_GPU(nu0, Nw, Ng, Nwk, Ngk, yLW_loc, ymG_vec, DzW, dzWK, DyG, dyGK, 0, zmW_vec, z_ground, 'GWall', 'ImpSurf');
ZLW_LW = GetScatteringMatrix_v2_GPU(nu0, Nw, Nw, Nwk, Nwk, yLW_loc, yLW_loc, DzW, dzWK, DzW, dzWK, 0, zmW_vec, zmW_vec, 'GWall', 'GWall');
ZLW_RW = GetScatteringMatrix_v2_GPU(nu0, Nw, Nw, Nwk, Nwk, yLW_loc, yRW_loc, DzW, dzWK, DzW, dzWK, 0, zmW_vec, zmW_vec, 'GWall', 'GWall');

ZRW_G  = GetScatteringMatrix_v2_GPU(nu0, Nw, Ng, Nwk, Ngk, yRW_loc, ymG_vec, DzW, dzWK, DyG, dyGK, 0, zmW_vec, z_ground, 'GWall', 'ImpSurf');
ZRW_LW = GetScatteringMatrix_v2_GPU(nu0, Nw, Nw, Nwk, Nwk, yRW_loc, yLW_loc, DzW, dzWK, DzW, dzWK, 0, zmW_vec, zmW_vec, 'GWall', 'GWall');
ZRW_RW = GetScatteringMatrix_v2_GPU(nu0, Nw, Nw, Nwk, Nwk, yRW_loc, yRW_loc, DzW, dzWK, DzW, dzWK, 0, zmW_vec, zmW_vec, 'GWall', 'GWall');

% (2) 고정 구조체 -> strip
Z1_4_cell = cell(1, Nvkz);
for nz = 1:Nvkz
    Z1_4_cell{1, nz} = GetScatteringMatrix_v2_GPU(nu0, Nm, Nv, Nk, Nvky, ym1_vec, ymV_vec, dyM, dyMK, DyV, dyVK, dzVK, 0, zmV4K_vec(nz), 'ImpSurf', 'Dielectric');
end
Z1_4 = Cell2MatGPU(Z1_4_cell);

Z1_G   = GetScatteringMatrix_v2_GPU(nu0, Nm, Ng, Nk, Ngk, ym1_vec, ymG_vec, dyM, dyMK, DyG, dyGK, 0, 0, z_ground, 'ImpSurf', 'ImpSurf');
Z1_LW  = GetScatteringMatrix_v2_GPU(nu0, Nm, Nw, Nk, Nwk, ym1_vec, yLW_loc, dyM, dyMK, DzW, dzWK, 0, 0, zmW_vec, 'ImpSurf', 'GWall');
Z1_RW  = GetScatteringMatrix_v2_GPU(nu0, Nm, Nw, Nk, Nwk, ym1_vec, yRW_loc, dyM, dyMK, DzW, dzWK, 0, 0, zmW_vec, 'ImpSurf', 'GWall');

% (3) strip -> 고정 구조체
Z4_1_cell = cell(Nvkz, 1);
for mz = 1:Nvkz
    Z4_1_cell{mz, 1} = GetScatteringMatrix_v2_GPU(nu0, Nv, Nm, Nvky, Nk, ymV_vec, ym1_vec, DyV, dyVK, dyM, dyMK, 0, zmV4K_vec(mz), 0, 'Dielectric', 'ImpSurf');
end
Z4_1 = Cell2MatGPU(Z4_1_cell);

ZG_1   = GetScatteringMatrix_v2_GPU(nu0, Ng, Nm, Ngk, Nk, ymG_vec, ym1_vec, DyG, dyGK, dyM, dyMK, 0, z_ground, 0, 'ImpSurf', 'ImpSurf');

ZLW_1  = GetScatteringMatrix_v2_GPU(nu0, Nw, Nm, Nwk, Nk, yLW_loc, ym1_vec, DzW, dzWK, dyM, dyMK, 0, zmW_vec, 0, 'GWall', 'ImpSurf');

ZRW_1  = GetScatteringMatrix_v2_GPU(nu0, Nw, Nm, Nwk, Nk, yRW_loc, ym1_vec, DzW, dzWK, dyM, dyMK, 0, zmW_vec, 0, 'GWall', 'ImpSurf');

% (4) strip 자기 상호작용
Z1_1   = GetScatteringMatrix_v2_GPU(nu0, Nm, Nm, Nk, Nk, ym1_vec, ym1_vec, dyM, dyMK, dyM, dyMK, 0, 0, 0, 'ImpSurf', 'ImpSurf');

% 블록 조립:  [A B; C D] [Is; Ii] = [Vs; Vi]
Z_cav = [Z4_4   Z4_G   Z4_LW   Z4_RW; ...
         ZG_4   ZG_G   ZG_LW   ZG_RW; ...
         ZLW_4  ZLW_G  ZLW_LW  ZLW_RW; ...
         ZRW_4  ZRW_G  ZRW_LW  ZRW_RW];

Bmtx = [Z4_1; ...
        ZG_1; ...
        ZLW_1; ...
        ZRW_1];

Cmtx = [Z1_4  Z1_G  Z1_LW  Z1_RW];

ZDmtx = Z1_1;

n4  = Nv * Nvky * Nvkz;
nG  = Ng * Ngk;
nLW = Nw * Nwk;
nRW = Nw * Nwk;

P_4_gpu   = (1 / (1j * omega0 * (eps_sub - 1) * epsilon0)) * eye(n4,'double','gpuArray');
EtaG_gpu  = zeros(nG,  nG,  'double','gpuArray');
EtaLW_gpu = zeros(nLW, nLW, 'double','gpuArray');
EtaRW_gpu = zeros(nRW, nRW, 'double','gpuArray');
Eta_cavity_gpu = blkdiag(P_4_gpu, EtaG_gpu, EtaLW_gpu, EtaRW_gpu);

Amtx_gpu = Eta_cavity_gpu + Z_cav;

% 소스 목록 (위치, 계수):  고정 구조체 / 능동 strip
n_src_V4 = Nv * Nvky * Nvkz;
r0_V = sqrt(dyVK * dzVK / pi);
y_fx = [repmat(ymVK_vec(:), Nvkz, 1); ymGK_vec(:); yLW_loc*ones(Nw*Nwk,1); yRW_loc*ones(Nw*Nwk,1)];
z_fx = [repelem(zmV4K_vec(:), Nv*Nvky, 1); z_ground*ones(Ng*Ngk,1); zmWK_vec(:); zmWK_vec(:)];
c_fx = [-(eta0*pi*r0_V/2) * besselj(1, k0*r0_V) * ones(n_src_V4,1); ...   % 유전체 (Richmond)
        -(eta0*k0/4) * dyGK * ones(Ng*Ngk,1); ...                         % 바닥
        -(eta0*k0/4) * dzWK * ones(2*Nw*Nwk,1)];                          % 좌우 벽 (안쪽 면)
y_ac = ymMK_vec(:);
z_ac = zeros(Nm*Nk,1);
c_ac = -(eta0*k0/4) * dyMK * ones(Nm*Nk,1);                               % strip (z = 0)

%% =========================================================================
% [AIE] 막힌 개구면 (z = z_ap PEC) image 보정 + 개구면 미지수 추가
% =========================================================================
tAIE = tic;
n_fx = numel(y_fx);  n_ac = numel(y_ac);
assert(numel(Vs) == n_fx && numel(Vi) == n_ac, 'source list / RHS size mismatch');
y_all_ = [y_fx; y_ac];  z_all_ = [z_fx; z_ac];  c_all_ = [c_fx; c_ac];

% (0) 규약 검사: 기존 MoM 블록의 먼 원소 = -c_n H0(kR) 인지 (바닥, 유전체)
nGk = Ng*Ngk;  nVk = Nv*Nvky;
Zg_ref = gather(ZG_G);
Zg_new = -green2D_points(ymGK_vec(:), z_ground*ones(nGk,1), ymGK_vec(:), z_ground*ones(nGk,1), c_fx(n_src_V4+1)*ones(nGk,1), k0);
farG   = abs(ymGK_vec(:) - ymGK_vec(:).') > 3*dyGK;
Zd_ref = gather(Z4_4(1:nVk, 1:nVk));
Zd_new = -green2D_points(ymVK_vec(:), zmV4K_vec(1)*ones(nVk,1), ymVK_vec(:), zmV4K_vec(1)*ones(nVk,1), c_fx(1)*ones(nVk,1), k0);
farD   = abs(ymVK_vec(:) - ymVK_vec(:).') > 3*dyVK;
errG = norm(Zg_new(farG) - Zg_ref(farG))/norm(Zg_ref(farG));
errD = norm(Zd_new(farD) - Zd_ref(farD))/norm(Zd_ref(farD));
fprintf('[AIE] MoM convention check: ground %.2e, dielectric %.2e\n', errG, errD);
assert(errG < 1e-2 && errD < 1e-2, 'MoM 규약 불일치: image 보정을 기존 행렬에 더할 수 없습니다.');
clear Zg_ref Zg_new Zd_ref Zd_new farG farD

% (1) 모든 블록에 image 보정: Z_int = Z_free + c_n H0(k |r_m - r_n,img|)
Zimg = green2D_points(y_all_, z_all_, y_all_, 2*z_ap - z_all_, c_all_, k0);
Amtx_gpu = Amtx_gpu + gpuArray(Zimg(1:n_fx, 1:n_fx));
Bmtx     = Bmtx     + gpuArray(Zimg(1:n_fx, n_fx+1:end));
Cmtx     = Cmtx     + gpuArray(Zimg(n_fx+1:end, 1:n_fx));
ZDmtx    = ZDmtx    + gpuArray(Zimg(n_fx+1:end, n_fx+1:end));
clear Zimg
Vs = Vs - const_inc * besselh(0, 2, k0 * hypot(y_fx - yi, z_fx - (2*z_ap - z_feed)));   % 피드 image
Vi = Vi - const_inc * besselh(0, 2, k0 * hypot(y_ac - yi, z_ac - (2*z_ap - z_feed)));

% (2) 개구면 펄스 -> 안쪽 모든 점의 E (2D Rayleigh-Sommerfeld, 아래로)
[tq, wq] = gauss_legendre(24, -D_ap/2, D_ap/2);
Dm_all = zeros(n_fx + n_ac, Na_ap);
for n = 1:Na_ap
    yp  = y_ap(n) + tq.';
    Zr  = z_ap - z_all_;
    rho = max(hypot(y_all_ - yp, Zr), 1e-30);
    Dm_all(:, n) = sum(wq.' .* (-1i*k0/2) .* (Zr./rho) .* besselh(1, 2, k0*rho), 2);
end
% (3) 안쪽 전류 / 피드 -> 개구면 H_y (image 로 2 배)
rho    = hypot(y_ap(:) - y_all_.', z_ap - z_all_.');
Hs_all = 2*(1i/(omega0*mu0)) * c_all_.' .* (-k0) .* besselh(1, 2, k0*rho) .* (z_ap - z_all_.')./rho;
rf     = hypot(y_ap(:) - yi, z_ap - z_feed);
Hf_ap  = 2*(1i/(omega0*mu0)) * const_inc * (-k0) .* besselh(1, 2, k0*rf) .* (z_ap - z_feed)./rf;
% (4) 개구면 자기 연산자: 바깥 (3D 플랜지, x 평균) - 안쪽 (2D)
lags    = (0:Na_ap).' * D_ap;
Hext_ap = -(2/(1i*omega0*mu0)) * a_operator(pk3(lags, D_ap, Lx_3D, k0), Na_ap, D_ap, k0);
Hint_ap =  (2/(1i*omega0*mu0)) * a_operator(pk2(lags, D_ap, k0, gam_E), Na_ap, D_ap, k0);
% 유한 플랜지 보정:  Hext = Hext(무한 플랜지) + (유한 - 무한)
Hb_ff = z_ap - z_ground;                       % 단면 금속 몸체 높이 (개구면 ~ 바닥)
FF    = ff_prepare(Lx_3D, t_flange, Hb_ff, ff_nslot, ff_dl, ff_beta_Y, k0, omega0, mu0, epsilon0, z_f - z_ap, ff_cache);
dH_ff = ff_dH_matrix(FF, Na_ap, D_ap);
fprintf('[AIE-FF] flange correction: ||dH|| / ||Hext_inf|| = %.3f\n', norm(dH_ff)/norm(Hext_ap));
Hext_ap = Hext_ap + dH_ff;
Yap     = Hext_ap - Hint_ap;
% (5) 확장: 고정 미지수 = [구조 전류 (n_fx); 개구면 e (Na_ap)].  H 식은 eta0 를 곱해 크기를 맞춤
Amtx_gpu = [Amtx_gpu, gpuArray(-Dm_all(1:n_fx, :)); gpuArray(-eta0*Hs_all(:, 1:n_fx)), gpuArray(eta0*Yap)];
Bmtx     = [Bmtx; gpuArray(-eta0*Hs_all(:, n_fx+1:end))];
Cmtx     = [Cmtx, gpuArray(-Dm_all(n_fx+1:end, :))];
Vs       = [Vs; eta0*Hf_ap];
ie_ap    = n_fx + (1:Na_ap);                    % 확장 고정 벡터 안의 개구면 위치
fprintf('[AIE] image + aperture blocks: %.1f s\n', toc(tAIE));

y_src = [y_fx; y_ac];  z_src = [z_fx; z_ac];  c_src = [c_fx; c_ac];

%% =========================================================================
% 6. Schur 축약 (strip 전류 Ii 만 남김)
% =========================================================================
%   Is = Is0 - T_Is*Ii,   (S0 + diag(Z)) Ii = rhs_Ii
Vi_gpu = gpuArray(Vi);
Vs_gpu = gpuArray(Vs);

T_Is_gpu = Amtx_gpu \ Bmtx;
Is0_gpu  = Amtx_gpu \ Vs_gpu;

S0_gpu     = ZDmtx - Cmtx * T_Is_gpu;
rhs_Ii_gpu = Vi_gpu - Cmtx * Is0_gpu;

% 개구면 e = H_ap*Ii + E0_ap
H_ap_gpu  = -T_Is_gpu(ie_ap, :);
E0_ap_gpu = Is0_gpu(ie_ap);

% 초점 field 는 개구면 e 에서만 나옴:  E_focus = Hred*Ii + Ebase
P_focus_gpu = gpuArray(ff_focal_matrix(x_obs(:), FF.Fzf, FF, y_ap, D_ap, ff_beta_F*k0));
Hred_gpu    = P_focus_gpu * H_ap_gpu;
Ebase_gpu   = P_focus_gpu * E0_ap_gpu;
clear P_focus_gpu

fprintf('Section 5-6 GPU precompute elapsed: %.3f s\n', toc(tPrecompute));

% 가장자리 벌점용 마스크
edge_mask     = abs(y_ap(:)) > Width/2 - edge_width;
edge_mask_gpu = gpuArray(double(edge_mask));
fprintf('[edge] w_edge = %.3g, %d of %d aperture points in edge zone\n', ...
    w_edge, nnz(edge_mask), numel(edge_mask));

% 피드 위치 산란 field (손실 벌점의 P_in 계산용):  E_sc(feed) = a_feed.'*Ii + b_feed
%   구조 전류 (image 포함) + 개구면 (아래로 RS) + 피드 image
g_s   = @(ys_, zs_, cs_) cs_(:).' .* (besselh(0, 2, k0*hypot(yi - ys_(:).', z_feed - zs_(:).')) ...
                                    - besselh(0, 2, k0*hypot(yi - ys_(:).', z_feed - (2*z_ap - zs_(:).'))));
rrf   = hypot(yi - (y_ap(:) + tq.'), z_ap - z_feed);
g_e   = sum(wq.' .* (-1i*k0/2) .* ((z_ap - z_feed)./rrf) .* besselh(1, 2, k0*rrf), 2).';
g_fxe = [g_s(y_fx, z_fx, c_fx), g_e];
g_acf = g_s(y_ac, z_ac, c_ac);
b_img = const_inc * (-besselh(0, 2, k0*2*abs(z_ap - z_feed)));
a_feed_gpu = (gpuArray(g_acf) - gpuArray(g_fxe) * T_Is_gpu).';
b_feed     = gather(gpuArray(g_fxe) * Is0_gpu) + b_img;
fprintf('[loss] w_loss = %.3g\n', w_loss);

%% =========================================================================
% 7. Adjoint 최적화 (연속 커패시턴스)
% =========================================================================
nvars = Nm;
target_gpu = gpuArray(target_PSF_lin(:));

% (C, R, X) 테이블 -> pchip 보간 (C^1) 과 그 도함수
cap_data = double(x_fit(:));
R_data = double(y_res_fit(:));
X_data = double(y_imp_fit(:));
[cap_data, order] = sort(cap_data);
R_data = R_data(order);
X_data = X_data(order);
[cap_data, uniq] = unique(cap_data, 'stable');
R_data = R_data(uniq);
X_data = X_data(uniq);
assert(numel(cap_data)>=2, 'Capacitance fitting needs >=2 unique points.');

if ~optimize_R
    R_data(:) = 0;
end

ppR = pchip(cap_data, R_data);
ppX = pchip(cap_data, X_data);
ppRd = ppR; ppXd = ppX;
orderR = ppR.order;
orderX = ppX.order;
ppRd.coefs = ppR.coefs(:,1:end-1) .* (orderR-1:-1:1);
ppXd.coefs = ppX.coefs(:,1:end-1) .* (orderX-1:-1:1);
ppRd.order = orderR-1;
ppXd.order = orderX-1;

Cmin = cap_data(1);
Cspan = cap_data(end)-Cmin;

% 초기값: 저장된 Nm 개 커패시턴스 인덱스 (""  또는 파일 없음 = 범위 중앙)
initial_file = "";
if isfile(initial_file)
    init_data = load(initial_file,'X');
    X_seed = round(double(init_data.X(:)));
    X_seed = max(1,min(numel(x_fit),X_seed));
    C_initial = double(x_fit(X_seed(:)));
    u0 = (C_initial(:)-Cmin)/Cspan;
    fprintf('Adjoint start: %s\n',initial_file);
else
    u0 = 0.5*ones(nvars,1);
    fprintf('Adjoint start: center of capacitance range (initial file absent).\n');
end
u0 = max(0,min(1,u0));

w_shape = 1.0;
w_centroid = 0.05;
y_norm_gpu = gpuArray((x_obs(:)-delta_y)/lambda0);

obj_adj = @(u) Cost_Adjoint_1Layer_GPU(u, ...
    S0_gpu, rhs_Ii_gpu, Hred_gpu, Ebase_gpu, ...
    target_gpu, y_norm_gpu, Nm, Nk, ...
    ppR, ppX, ppRd, ppXd, Cmin, Cspan, w_shape, w_centroid, ...
    H_ap_gpu, E0_ap_gpu, edge_mask_gpu, w_edge, ...
    a_feed_gpu, b_feed, Winc, I0, dyMK, w_loss);

if strlength(replay_X_file) == 0
    % 초기점에서 adjoint gradient 를 중앙차분과 비교
    [f0,g0] = obj_adj(u0);
    assert(isfinite(f0) && all(isfinite(g0)), 'Initial objective/gradient is nonfinite.');
    idx_test = unique(round(linspace(1,nvars,min(nvars,6))));
    h_fd = 1e-5;
    fprintf('Adjoint gradient check (central differences):\n');
    for j = idx_test
        up=u0; um=u0;
        up(j)=min(1,up(j)+h_fd);
        um(j)=max(0,um(j)-h_fd);
        fp = obj_adj(up); fm = obj_adj(um);
        g_fd = (fp-fm)/(up(j)-um(j));
        relerr = abs(g_fd-g0(j))/max([1,abs(g_fd),abs(g0(j))]);
        fprintf('  var %02d: adj %+.5e | FD %+.5e | relative %.3e\n',j,g0(j),g_fd,relerr);
    end

    adj_options = optimoptions('fmincon', ...
        'Algorithm','interior-point', ...
        'SpecifyObjectiveGradient',true, ...
        'Display','iter', ...
        'MaxIterations',300, ...
        'MaxFunctionEvaluations',3000, ...
        'OptimalityTolerance',1e-7, ...
        'StepTolerance',1e-10, ...
        'UseParallel',false);

    tAdj = tic;
    [u_opt,fval,exitflag,output] = fmincon(obj_adj,u0,[],[],[],[], ...
        zeros(nvars,1),ones(nvars,1),[],adj_options);
    optimization_elapsed=toc(tAdj);
    fprintf('Adjoint complete: cost %.6e, exitflag %d, time %.2f s\n', ...
        fval,exitflag,optimization_elapsed);

    % 연속 C -> 가장 가까운 테이블 인덱스
    C_optimized = Cmin + Cspan*u_opt(:);
    X = zeros(1,nvars);
    for j = 1:nvars
        [~,X(j)] = min(abs(double(x_fit(:))-C_optimized(j)));
    end
else
    % 최적화 생략: 저장된 X 로 결과만 계산
    rp = load(replay_X_file, 'X');
    X  = round(double(rp.X(:).'));
    fprintf('[replay] X loaded from %s\n', replay_X_file);
end

%% =========================================================================
% 8. 최적 결과 복원 및 시각화 (양자화된 인덱스 X 사용)
% =========================================================================
X_indices = round(X);
X_indices = max(1, min(length(x_fit), X_indices));

idx_layer1 = X_indices(1:Nm);

if optimize_R
    R_opt1 = y_res_fit(idx_layer1);
else
    R_opt1 = zeros(1, Nm);
end

X_opt1 = y_imp_fit(idx_layer1);
x_ctrl_opt1 = x_fit(idx_layer1);

Z_strip_opt1 = R_opt1 + 1j * X_opt1;
Z_s_expanded = repelem(Z_strip_opt1(:), Nk);

% Schur complement 로 전류 복원
Smtx_opt_gpu    = S0_gpu + diag(gpuArray(Z_s_expanded));
Ii_opt_gpu      = Smtx_opt_gpu \ rhs_Ii_gpu;
Is_opt_gpu      = Is0_gpu - T_Is_gpu * Ii_opt_gpu;
I_total_opt_gpu = [Is_opt_gpu; Ii_opt_gpu];

Ii_opt      = gather(Ii_opt_gpu);
I_total_opt = gather(I_total_opt_gpu);

% 오믹 손실
P_ohmic = sum(0.5 * real(Z_s_expanded) .* abs(Ii_opt).^2 * dyMK);

% 방사 효율 (Winc 기준)
eR    = 1 - P_ohmic / Winc;
eR_dB = 10 * log10(max(eR, 1e-12));

% 실제 입력 전력 기준 효율 + 에너지 검사 (P_in = P_ohmic + P_ap)
E_ap_opt_gpu = Is_opt_gpu(ie_ap);
E_ap_opt     = gather(E_ap_opt_gpu);
E_sc_feed = [g_fxe, g_acf] * I_total_opt + b_img;
P_in      = Winc - 0.5 * real(E_sc_feed * conj(I0));
eR_true   = 1 - P_ohmic / P_in;
P_ap      = 0.5 * real(sum(E_ap_opt .* conj(Hext_ap * E_ap_opt)) * D_ap);
fprintf(' - P_in / W_inc                          : %.4f\n', P_in / Winc);
fprintf(' - Radiation efficiency (실제 입력 기준)  : %.2f %%\n', 100 * eR_true);
fprintf(' - 검증: P_ap (개구면 통과)/W_inc = %.4f  vs  (P_in - P_ohmic)/W_inc = %.4f\n', ...
    P_ap / Winc, (P_in - P_ohmic) / Winc);

% 넓은 범위 초점면 field
x_obs_wide = linspace(-30*lambda0, 30*lambda0, N_obs);
E_total_focal = ff_focal_matrix(x_obs_wide(:), FF.Fzf, FF, y_ap, D_ap, ff_beta_F*k0) * E_ap_opt;

target_PSF_ideal = interp1(y_coords, E_focus_ideal_abs, x_obs_wide, 'pchip', 0);

PSF_opt_abs = abs(E_total_focal);
P_dB = 20 * log10(PSF_opt_abs / max(PSF_opt_abs));

% 초점 효율 및 빔폭
P_focal_density = abs(E_total_focal(:)).^2;
P_peak = max(P_focal_density);
P_th_3dB = 0.5 * P_peak;
[~, idx_peak] = max(P_focal_density);
idx_left = find(P_focal_density(1:idx_peak) <= P_th_3dB, 1, 'last');
if isempty(idx_left), idx_left = 1; end
idx_right_rel = find(P_focal_density(idx_peak:end) <= P_th_3dB, 1, 'first');
if isempty(idx_right_rel)
    idx_right = length(P_focal_density);
else
    idx_right = idx_peak + idx_right_rel - 1;
end
P_spot   = trapz(x_obs_wide(idx_left:idx_right), P_focal_density(idx_left:idx_right));
P_total  = trapz(x_obs_wide, P_focal_density);
eta_foc  = (P_spot / P_total) * 100;
FWHM_val = (x_obs_wide(idx_right) - x_obs_wide(idx_left)) / lambda0;
fprintf(' - 3 dB Spot Width (FWHM) : %.3f \\lambda_0\n', FWHM_val);
fprintf(' - Focusing Efficiency (\\eta_foc): %.2f %%\n', eta_foc);
fprintf(' - Radiation Efficiency (eR)      : %.2f %% (%.2f dB)\n', eR*100, eR_dB);

dx_exact_theory      = lambda0 / ( (delta_y + Width/2)/sqrt((delta_y + Width/2)^2 + z_f^2) - (delta_y - Width/2)/sqrt((delta_y - Width/2)^2 + z_f^2) );
dx_exact_FWHM_theory = 0.4429 * lambda0 / NA_num;
title_str = sprintf("DL FWHM: %.3f mm, Optimized FWHM: %.3f mm, Ratio: %.2f %%", ...
    dx_exact_FWHM_theory * 1e3, FWHM_val * lambda0 * 1e3, FWHM_val*lambda0/dx_exact_FWHM_theory*100);

% Figure 1: 초점면 패턴 (dB)
figure(1); clf;
plot(x_obs_wide / lambda0, 20*log10(target_PSF_ideal) - max(20*log10(target_PSF_ideal)), '--k', 'LineWidth', LW, 'DisplayName', 'Target (Ideal Source)');
hold on;
plot(x_obs_wide / lambda0, P_dB, '-b', 'LineWidth', LW, 'DisplayName', 'Cavity 1-Layer MS Opt');
grid on;
xline(-y_ROI / lambda0 + delta_y / lambda0, '--r', "LineWidth", 3, 'HandleVisibility', 'off');
xline( y_ROI / lambda0 + delta_y / lambda0, '--r', "LineWidth", 3, 'HandleVisibility', 'off');
xlim([-10 10]); ylim([-30 0]);
xlabel('${y}~/~\lambda_0$', 'Interpreter', 'latex');
ylabel('$|\mathbf{E}|$ (dB)', 'Interpreter', 'latex', 'FontSize', FS);
title(title_str, "FontName", Fonts, "FontSize", FS);
legend('Interpreter', 'latex', 'FontSize', FS-2, 'Location', 'SouthWest');
ax = gca;
ax.XAxis.FontName = Fonts; ax.YAxis.FontName = Fonts;
ax.XAxis.FontSize = FS;   ax.YAxis.FontSize = FS;
ax.Box = 'on'; ax.LineWidth = 1.5; ax.XColor = 'k'; ax.YColor = 'k';

% Figure 2: 초점면 |E| (V/m)
figure(2); clf;
hold on;
plot(x_obs_wide / lambda0, PSF_opt_abs, '-b', 'LineWidth', LW, 'DisplayName', 'Cavity 1-Layer MS Opt');
grid on;
xline(-y_ROI / lambda0 + delta_y / lambda0, '--r', "LineWidth", 3, 'HandleVisibility', 'off');
xline( y_ROI / lambda0 + delta_y / lambda0, '--r', "LineWidth", 3, 'HandleVisibility', 'off');
xlim([-10 10]); ylim([0 max(PSF_opt_abs)*1.1]);
xlabel('${y}~/~\lambda_0$', 'Interpreter', 'latex');
ylabel('$|\mathbf{E}|$ (V/m)', 'Interpreter', 'latex', 'FontSize', FS);
title(title_str, "FontName", Fonts, "FontSize", FS);
legend('Interpreter', 'latex', 'FontSize', FS-2, 'Location', 'NorthWest');
ax = gca;
ax.XAxis.FontName = Fonts; ax.YAxis.FontName = Fonts;
ax.XAxis.FontSize = FS;   ax.YAxis.FontSize = FS;
ax.Box = 'on'; ax.LineWidth = 1.5; ax.XColor = 'k'; ax.YColor = 'k';

% Figure 90: 초점면 세기 (선형, 정규화)
figure(90); clf;
plot(x_obs_wide / lambda0, abs(target_PSF_ideal).^2 / max(abs(target_PSF_ideal).^2), '--k', 'LineWidth', LW, 'DisplayName', 'Target (Ideal Source)');
hold on;
plot(x_obs_wide / lambda0, PSF_opt_abs.^2 / max(PSF_opt_abs.^2), '-b', 'LineWidth', LW, 'DisplayName', 'Cavity 1-Layer MS Opt');
grid on;
xline(-y_ROI / lambda0 + delta_y / lambda0, '--r', "LineWidth", 3, 'HandleVisibility', 'off');
xline( y_ROI / lambda0 + delta_y / lambda0, '--r', "LineWidth", 3, 'HandleVisibility', 'off');
xlim([-10 10]); ylim([0 1]);
ax = gca;
ax.XAxis.FontName = Fonts; ax.YAxis.FontName = Fonts;
ax.XAxis.FontSize = FS;   ax.YAxis.FontSize = FS;
ax.Box = 'on'; ax.LineWidth = 1.5; ax.XColor = 'k'; ax.YColor = 'k';
xlabel('${y}~/~\lambda_0$', 'Interpreter', 'latex');
ylabel('$|\mathbf{E}|^2$ (norm.)', 'Interpreter', 'latex', 'FontSize', FS);
legend('Interpreter', 'latex', 'FontSize', FS-2, 'Location', 'northeast');

fprintf('\n==========================================================\n');
fprintf(' - 엄밀한 이론적 해상도 (Exact NA)   : %8.3f mm (%.3f lambda0)\n', dx_exact_theory * 1e3, dx_exact_theory / lambda0);
fprintf(' - 엄밀한 이론적 FWHM                : %8.3f mm (%.3f lambda0)\n', dx_exact_FWHM_theory * 1e3, dx_exact_FWHM_theory / lambda0);
fprintf(' - 최적화된 메타표면 PSF FWHM (D_dB)  : %8.3f mm (%.3f lambda0), 비율 : %.2f %%\n', FWHM_val * lambda0 * 1e3, FWHM_val, FWHM_val*lambda0/dx_exact_FWHM_theory*100);
fprintf('----------------------------------------------------------\n');

% 최적 R / X / C 출력 (HFSS 입력용)
R_ov = R_opt1;
X_ov = X_opt1;
C_ov = x_ctrl_opt1;

fprintf('str_res = "%s"\n', strjoin(arrayfun(@(x) sprintf('%.2f', x), R_ov, 'UniformOutput', false), ', '));
fprintf('str_imp = "%s"\n', strjoin(arrayfun(@(x) sprintf('%.2f', x), X_ov, 'UniformOutput', false), ', '));
fprintf('str_cap = "%s"\n', strjoin(arrayfun(@(x) sprintf('%.3f', x), C_ov, 'UniformOutput', false), ' '));

all_vals = [R_ov(:).', X_ov(:).'];
all_str = strjoin(arrayfun(@(x) sprintf('%.2f', x), all_vals, 'UniformOutput', false), '\t');
clipboard('copy', all_str);

%% 결과 저장
edge_frac = sum(abs(E_ap_opt(edge_mask)).^2) / sum(abs(E_ap_opt).^2);
fprintf('[edge] aperture energy in edge zone (|y| > W/2 - %.2f lam0): %.1f %%\n', edge_width/lambda0, 100*edge_frac);

if ~exist(save_folder, 'dir'), mkdir(save_folder); end
if abs(delta_y) < 1e-12, dy_str = '0p00';
elseif delta_y < 0,       dy_str = ['m', strrep(sprintf('%.2f', abs(delta_y/lambda0)), '.', 'p')];
else,                     dy_str = strrep(sprintf('%.2f', delta_y/lambda0), '.', 'p');
end

file_name = sprintf('X_1L_H_%slam.mat', dy_str);

if strlength(replay_X_file) == 0
    save(fullfile(save_folder, file_name), 'X');
end

[~, base_name] = fileparts(file_name);
save(fullfile(save_folder, sprintf('PSF_for_HFSS_compare_%s.mat', base_name)), 'x_obs_wide', 'E_total_focal', 'delta_y', ...
     'lambda0', 'z_f', 'z_ap', 'Lx_3D', 'H_baffle', 'R_ov', 'X_ov', 'X', 'w_edge', 'w_loss', 'edge_frac', 'eR_true');
fprintf('Saved to folder: %s\n', save_folder);

%% =========================================================================
% 9. yz 평면 (x = 0) |E| 지도 (HFSS yz 평면 field plot 과 같은 단면)
%   z <= z_ap (캐비티 안)     : 피드 + 구조 전류 (image 포함) + 개구면 RS (아래로)
%   z_ap < z < z_ap + lam0/4  : 개구면 -> 3D (무한 플랜지 근사)
%   그 위                     : 유한 플랜지 ky 적분
%   캐비티 바깥 아래쪽은 계산하지 않음 (NaN)
% =========================================================================

%{

yz_y_half  = Width/2 + 1.0*lambda0;      % y 범위 +-
yz_z_min   = z_ground - 0.5*lambda0;     % z 아래 끝
yz_z_max   = z_f + 3.0*lambda0;          % z 위 끝
yz_Ny      = 241;
yz_Nz      = 321;
yz_clim_k  = 1.5;                        % 색 상한 = yz_clim_k * (z_f 선 위 최대 |E|)

yz_y = linspace(-yz_y_half, yz_y_half, yz_Ny);
yz_z = linspace(yz_z_min, yz_z_max, yz_Nz);
[YY, ZZ] = meshgrid(yz_y, yz_z);
Yf = YY(:); Zf = ZZ(:);
E_yz = nan(numel(Yf), 1);

idx_2d = find(Zf <= z_ap & Zf >= z_ground & abs(Yf) <= Width/2);
idx_3d = find(Zf >  z_ap & Zf < z_ap + lambda0/4);
idx_ff = find(Zf >= z_ap + lambda0/4);
I_str  = [I_total_opt(1:n_fx); I_total_opt(n_fx+Na_ap+1:end)];
for b = 1:3000:numel(idx_2d)
    ii = idx_2d(b:min(b+2999, numel(idx_2d)));
    E_yz(ii) = interior_field(Yf(ii), Zf(ii), y_src, z_src, c_src, I_str, y_ap, z_ap, E_ap_opt, tq, wq, ...
                              yi, z_feed, const_inc, k0);
end
E_yz(idx_3d) = gather(ap_matrix(Yf(idx_3d), Zf(idx_3d), y_ap, D_ap, z_ap, Lx_3D, k0, 2048) * E_ap_opt_gpu);
[zr, ~, iz_r] = unique(Zf(idx_ff));
Fz_yz = ff_fields_z(FF, zr - z_ap, ff_beta_Y*k0);                         % Nb x Nr
b_all = FF.beta;  c_all = FF.wts .* (D_ap*ff_sinc(b_all*D_ap/(2*pi))) / pi;
Cy    = cos(yz_y(:)*b_all).*(cos(b_all.'*y_ap)*E_ap_opt).' + sin(yz_y(:)*b_all).*(sin(b_all.'*y_ap)*E_ap_opt).';
E_rows = (Cy .* c_all) * Fz_yz;                                           % Ny x Nr
[~, iy_c] = min(abs(Yf(idx_ff) - yz_y), [], 2);
E_yz(idx_ff) = E_rows(sub2ind(size(E_rows), iy_c, iz_r));
A_yz = reshape(abs(E_yz), yz_Nz, yz_Ny);

[~, iz_f] = min(abs(yz_z - z_f));
c_top = yz_clim_k * max(A_yz(iz_f, :));

figure(4); clf;
h_img = imagesc(yz_y*1e3, yz_z*1e3, A_yz);
set(h_img, 'AlphaData', ~isnan(A_yz));
set(gca, 'YDir', 'normal', 'Color', [1 1 1]);
colormap(jet(256)); clim([0 c_top]); colorbar;
hold on;
line([yLW_loc yRW_loc]*1e3, [z_ground z_ground]*1e3, 'Color', 'k', 'LineWidth', 2);
line([yLW_loc yLW_loc]*1e3, [z_ground H_baffle]*1e3, 'Color', 'k', 'LineWidth', 2);
line([yRW_loc yRW_loc]*1e3, [z_ground H_baffle]*1e3, 'Color', 'k', 'LineWidth', 2);
yline(z_f*1e3, '-k', 'LineWidth', 1.2);
yline(z_ap*1e3, ':m', 'LineWidth', 1.2);
plot(yi*1e3, z_feed*1e3, 'wo', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
axis equal tight;
xlabel('y (mm)'); ylabel('z (mm)');
title(sprintf('|E| on x = 0 plane  (color max = %.1f x max on z_f line)', yz_clim_k));
set(gca, 'FontName', Fonts, 'FontSize', FS-4);

% 축 위 |E| (y = delta_y): 초점이 z_f 에 맺히는지 확인
[~, iy0] = min(abs(yz_y - delta_y));
z_axis_ok = yz_z > 2*lambda0;        % 개구면 근접장 제외
[~, im] = max(A_yz(z_axis_ok, iy0));
z_tmp = yz_z(z_axis_ok);
fprintf('yz map: |E| max on axis (above aperture) at z = %.1f mm (z_f = %.1f mm)\n', z_tmp(im)*1e3, z_f*1e3);

save(fullfile(save_folder, 'YZ_map.mat'), 'yz_y', 'yz_z', 'A_yz', 'z_f', 'z_ap', 'Lx_3D');

%}

%% =========================================================================
% Local functions
% =========================================================================
function [f, grad] = Cost_Adjoint_1Layer_GPU(u, S0, rhs, Hred, Ebase, ...
    target, y_centered, Nm, Nk, ppR, ppX, ppRd, ppXd, ...
    Cmin, Cspan, w_shape, w_centroid, Hap, E0ap, emask, w_edge, ...
    a_feed, b_feed, Winc, I0, dyMK, w_loss)
% 1층 메타표면 Schur complement 에 대한 discrete adjoint
%   State:    S(u)*Ii = rhs
%   Observer: E = Hred*Ii + Ebase
%   Adjoint:  S(u)'*lambda = Hred'*gE
%   Gradient: -2*real(lambda'*(dS/du_j)*Ii)     (규약 df = 2*real(g'*dE))
%   w_edge > 0 : f += w_edge * sum(emask.*|E_ap|^2)/sum(|E_ap|^2),  E_ap = Hap*Ii + E0ap
%   w_loss > 0 : f += w_loss * P_ohmic/P_in
%                P_ohmic = sum 0.5*Re(Z)|Ii|^2 dyMK,  P_in = Winc - 0.5*Re[(a_feed.'*Ii + b_feed) conj(I0)]

u = double(u(:));
C = Cmin + Cspan*u;  % pF
Z = ppval(ppR,C) + 1i*ppval(ppX,C);
dZdu = Cspan*(ppval(ppRd,C) + 1i*ppval(ppXd,C));

Zsub = repelem(Z,Nk);
S = S0 + diag(gpuArray(Zsub));
Ii = S \ rhs;
E = Hred*Ii + Ebase;

% PSF 모양: RMS 로 정규화한 |E| 와 target 의 차이
N = numel(E);
a = abs(E);
q = sqrt(sum(a.^2)/N + 1e-20);
m = a/q;
r = m-target;

fshape = sum(r.^2)/N;
ga = (r/(N*q)).*(E./max(a,1e-12));
proj = sum(r.*a)/(N*N*q^3);
ga = ga - proj*E;

% 세기 중심 (목표 초점 위치 기준)
P = sum(a.^2) + 1e-20;
center = sum(y_centered.*a.^2)/P;
fcenter = center^2;
gcenter = center*E.*(y_centered-center)/P;

f_gpu = real(w_shape*fshape + w_centroid*fcenter);
gE = w_shape*ga + w_centroid*gcenter;

% 가장자리 벌점
gA = [];
if w_edge > 0
    Eap = Hap*Ii + E0ap;
    pa  = abs(Eap).^2;
    Pt  = sum(pa) + 1e-30;
    Pe  = sum(emask.*pa);
    f_gpu = f_gpu + w_edge*real(Pe/Pt);
    gA  = w_edge*(emask.*Eap/Pt - (Pe/Pt^2)*Eap);
end

% 손실 벌점 (Ii 에 대한 gradient gI)
gI = [];
if w_loss > 0
    Rg    = gpuArray(real(Zsub(:)));
    P_ohm = sum(0.5*Rg.*abs(Ii).^2*dyMK);
    P_in  = Winc - 0.5*real((a_feed.'*Ii + b_feed)*conj(I0));
    f_gpu = f_gpu + w_loss*real(P_ohm/P_in);
    g_Po  = 0.5*dyMK*Rg.*Ii;
    g_Pi  = -0.25*I0*conj(a_feed);
    gI    = w_loss*(g_Po/P_in - (P_ohm/P_in^2)*g_Pi);
end

if nargout>1
    rhs_adj = Hred'*gE;
    if ~isempty(gA), rhs_adj = rhs_adj + Hap'*gA; end
    if ~isempty(gI), rhs_adj = rhs_adj + gI; end
    lambda = S' \ rhs_adj;
    sensitivities = reshape(conj(lambda).*Ii,Nk,Nm);
    grad_gpu = -2*real(sum(sensitivities,1).' .* gpuArray(dZdu(:)));
    if w_loss > 0
        % P_ohmic 이 Re(Z) 를 직접 포함하는 항 (전류 고정일 때의 미분)
        Isq = sum(reshape(abs(Ii).^2, Nk, Nm), 1).';
        grad_gpu = grad_gpu + (w_loss/P_in) * 0.5*dyMK * real(gpuArray(dZdu(:))) .* Isq;
    end
    grad = double(gather(grad_gpu));
end
f = double(gather(f_gpu));
end

function G = green2D_points(y_o, z_o, y_s, z_s, c_s, k0)
    % 2D: 소스 (y_s, z_s, 계수 c_s) -> 관측점 (y_o, z_o) 의 E_x.   G(o, s) = c_s * H0^(2)(k0 R)
    %   출력 No x Ns (CPU)
    R = sqrt((y_o(:) - y_s(:).').^2 + (z_o(:) - z_s(:).').^2);
    G = besselh(0, 2, k0 * R) .* c_s(:).';
end

%% ======================= [AIE] local functions =======================
function E = interior_field(Y, Z, y_s, z_s, c_s, I_s, y_ap, z_ap, e, tq, wq, yi, z_feed, ci, k0)
    % 막힌 개구면 아래 (캐비티 안) 의 E_x: 피드 (image 포함) + 구조 전류 (image 포함) + 개구면 RS
    Y = Y(:); Z = Z(:);
    E = ci*(besselh(0,2,k0*hypot(Y - yi, Z - z_feed)) - besselh(0,2,k0*hypot(Y - yi, Z - (2*z_ap - z_feed)))) ...
      + (green2D_points(Y, Z, y_s, z_s, c_s, k0) - green2D_points(Y, Z, y_s, 2*z_ap - z_s, c_s, k0)) * I_s;
    Zr = z_ap - Z;  Dm = zeros(numel(Y), numel(y_ap));
    for q = 1:numel(tq)
        rho = max(hypot(Y - (y_ap(:).' + tq(q)), Zr), 1e-30);
        Dm  = Dm + wq(q) * (-1i*k0/2) * (Zr./rho) .* besselh(1, 2, k0*rho);
    end
    E = E + Dm * e(:);
end

function P = ap_matrix(yo, zo, y_ap, D, z_ap, Lx, k0, chunkRows)
    % 플랜지 개구면 펄스 (중심 y_ap, 폭 D, |x| <= Lx/2) -> (x = 0, yo, zo) 의 E_x.  3D Rayleigh-Sommerfeld
    %   출력 numel(yo) x numel(y_ap) gpuArray
    yo = double(yo(:)); zo = double(zo(:));
    assert(all(zo > z_ap), 'ap_matrix: observation points must be above z_ap.');
    lam = 2*pi/k0;  nx = max(8, ceil(Lx/(lam/16)));
    xq = ((1:nx) - 0.5)/nx*Lx - Lx/2;  wx = Lx/nx;
    [tq, wq] = gauss_legendre(6, -D/2, D/2);
    ya = gpuArray(y_ap(:).');
    P  = complex(zeros(numel(yo), numel(y_ap), 'double', 'gpuArray'));
    for i0 = 1:chunkRows:numel(yo)
        ii = i0:min(i0+chunkRows-1, numel(yo));
        Zr = gpuArray(zo(ii) - z_ap);  yy = gpuArray(yo(ii));
        acc = complex(zeros(numel(ii), numel(y_ap), 'double', 'gpuArray'));
        for q = 1:numel(tq)
            Y2 = (yy - (ya + tq(q))).^2 + Zr.^2;
            for p = 1:nx
                R = sqrt(Y2 + xq(p)^2);
                acc = acc + wq(q) * (Zr./R) .* (1i*k0 + 1./R) .* exp(-1i*k0*R) ./ R;
            end
        end
        P(ii,:) = (wx/(2*pi)) * acc;
    end
end

function A = a_operator(Pv, Na, D, k)
    % (k^2 + d^2/dy^2) int K(y-y') e(y') dy',  펄스 기저, 펄스 중심 점정합 (경계 점프 -> 스칼라 퍼텐셜 차분)
    idx = abs((1:Na).' - (1:Na));
    V   = k^2 * D * Pv(idx + 1);
    B   = 0:Na;
    Pb  = Pv(abs(B.' - B) + 1);
    J   = zeros(Na+1, Na);
    for b = 0:Na
        if b < Na, J(b+1, b+1) = J(b+1, b+1) + 1; end
        if b > 0,  J(b+1, b)   = J(b+1, b)   - 1; end
    end
    Phi = Pb * J;
    A   = V + (Phi(2:end,:) - Phi(1:end-1,:)) / D;
end

function P = pk2(s, D, k, gam)
    % (1/D) int (-j/4) H0(k|s+t|) dt  (2D Green 의 펄스 평균)
    P = zeros(size(s));
    [t, wt] = gauss_legendre(16, -D/2, D/2);
    for i = 1:numel(s)
        if abs(s(i)) < 1e-12
            P(i) = (-1i/4)*(1 - 1i*(2/pi)*(log(gam*k*D/4) - 1));
        else
            P(i) = sum(wt .* (-1i/4) .* besselh(0,2,k*abs(s(i) + t))) / D;
        end
    end
end

function P = pk3(s, D, a, k)
    % 개구면 폭 a 로 x 평균한 3D Green 을 펄스 폭 D 로 y 평균 (정적 1/R 부분은 해석적으로)
    P = zeros(size(s));
    [u, wu] = gauss_legendre(64, 0, a);
    for i = 1:numel(s)
        edges = [-D/2, D/2];
        if -D/2 < -s(i) && -s(i) < D/2, edges = [-D/2, -s(i), D/2]; end
        tot = 0;
        for j = 1:numel(edges)-1
            [t, wt] = gauss_legendre(48, edges(j), edges(j+1));
            v  = max(abs(s(i) + t), 1e-30);
            st = 2*(a*asinh(a./v) - (sqrt(a^2 + v.^2) - v)) / (4*pi);
            R  = max(sqrt(u.'.^2 + v.^2), 1e-30);
            dy = 2*sum(wu.' .* (a - u.') .* (exp(-1i*k*R) - 1)./(4*pi*R), 2);
            tot = tot + sum(wt .* (st + dy));
        end
        P(i) = tot / (a*D);
    end
end

function [x, w] = gauss_legendre(n, a, b)
    i = (1:n-1).';  bta = i./sqrt(4*i.^2 - 1);
    [V, L] = eig(diag(bta, 1) + diag(bta, -1));
    [x, o] = sort(diag(L));  w = 2*V(1, o).'.^2;
    x = 0.5*(b - a)*x + 0.5*(b + a);  w = 0.5*(b - a)*w;
end

%% ======================= [AIE-FF] 유한 플랜지 바깥 계산 =======================
%  단면 (x-z) 금속 몸체: 슬롯 (폭 a) + 양옆 플랜지 (폭 t) + 바깥 옆면 + 바닥 (높이 Hb), y 방향으로는 일정
%  ky = beta 성분마다: 슬롯 위 자기전류 (E_ap = 1, x 로 일정) + 몸체 유도전류 (J_t, J_y) 를 2.5D MoM 으로 풂
%    Y(beta) : 슬롯 평균 H_y (바깥 쪽)  = 바깥이 개구면으로 되미는 양
%    F(beta) : (x = 0, z) 의 E_x         = 바깥으로 나가는 field

function FF = ff_prepare(a, t, Hb, nslot, dl, bYmax, k0, w, mu0, eps0, zf_rel, cache)
    if strlength(cache) == 0
        cache = sprintf('ffcache_a%.3f_t%.3f_H%.3f_n%d_dl%.3f_b%.1f_zf%.3f.mat', ...
            a*1e3, t*1e3, Hb*1e3, nslot, dl*1e3, bYmax, zf_rel*1e3);
    end
    if isfile(cache)
        S = load(cache, 'FF');  FF = S.FF;
        fprintf('[AIE-FF] exterior data loaded from %s (%d ky values)\n', cache, numel(FF.beta));
        return;
    end
    assert(t >= dl, '[AIE-FF] t_flange 는 ff_dl 이상이어야 합니다.');
    b1 = 0:(k0/60):(0.9*k0);
    n2 = 60;  b2 = 0.9*k0 + ((1:n2) - 0.5)*(0.2*k0/n2);          % k0 근처 촘촘히 (k0 자체는 피함)
    b3 = (1.1*k0 + k0/40):(k0/40):(bYmax*k0);
    beta = unique([b1, b2, b3]);
    wts  = zeros(size(beta));                                        % 사다리꼴 가중치
    db   = diff(beta);
    wts(1:end-1) = wts(1:end-1) + db/2;  wts(2:end) = wts(2:end) + db/2;

    geom = ff_contour(a, t, Hb, nslot, dl);
    Nb = numel(beta);  Ns = geom.Ns;
    FF.beta = beta;  FF.wts = wts;  FF.geom = geom;  FF.a = a;  FF.t = t;  FF.Hb = Hb;
    FF.k0 = k0;  FF.w = w;  FF.mu0 = mu0;  FF.eps0 = eps0;
    FF.Y = zeros(1, Nb);  FF.Yinf = zeros(1, Nb);  FF.Fzf = zeros(1, Nb);
    FF.Jt = zeros(Ns, Nb);  FF.Jy = zeros(Ns, Nb);
    tff = tic;
    for j = 1:Nb
        [Y, Jt, Jy] = ff_solve(beta(j), geom, a, k0, w, mu0, eps0);
        FF.Y(j) = Y;  FF.Jt(:,j) = Jt;  FF.Jy(:,j) = Jy;
        FF.Yinf(j) = ff_Yinf(beta(j), a, k0, w, mu0);
        FF.Fzf(j)  = ff_obs(geom, beta(j), Jt, Jy, [0, zf_rel], k0, w, mu0, eps0);
        if mod(j, 40) == 0 || j == Nb
            fprintf('[AIE-FF] exterior ky %d/%d (%.0f s)\n', j, Nb, toc(tff));
        end
    end
    save(cache, 'FF');
    fprintf('[AIE-FF] exterior data saved to %s\n', cache);
end

function geom = ff_contour(a, t, Hb, nslot, dl)
    % 시계 방향: 윗면 (왼 -> 오) , 오른 옆면 (아래로), 바닥 (오 -> 왼), 왼 옆면 (위로).  z = 0 = 개구면
    X = a/2 + t;  pts = zeros(0, 2);
    nt = max(1, ceil(t/dl));  nh = ceil(Hb/dl);  nb = ceil(2*X/dl);
    pts = [pts; line_pts([-X 0], [-a/2 0], nt)];
    pts = [pts; line_pts([-a/2 0], [a/2 0], nslot)];
    pts = [pts; line_pts([a/2 0], [X 0], nt)];
    pts = [pts; line_pts([X 0], [X -Hb], nh)];
    pts = [pts; line_pts([X -Hb], [-X -Hb], nb)];
    pts = [pts; line_pts([-X -Hb], [-X 0], nh)];
    geom.A0 = pts;  geom.A1 = circshift(pts, -1, 1);
    geom.Ns = size(pts, 1);
    geom.C  = 0.5*(geom.A0 + geom.A1);
    D = geom.A1 - geom.A0;  geom.L = sqrt(sum(D.^2, 2));  geom.T = D ./ geom.L;
    geom.slot = false(geom.Ns, 1);  geom.slot(nt+1:nt+nslot) = true;
    function P = line_pts(p0, p1, n)
        P = p0 + ((0:n-1).'/n) .* (p1 - p0);
    end
end

function kp = ff_kap(beta, k0)
    if beta < k0, kp = sqrt(k0^2 - beta^2); else, kp = -1i*sqrt(beta^2 - k0^2); end
end

function g = ff_G(R, kp)
    % 2.5D Green:  (1/4j) H0^(2)(kp R)  (kp 가 허수면 K0/(2 pi))
    if isreal(kp)
        x = kp*R;  g = (besselj(0, x) - 1i*bessely(0, x)) / (4i);
    else
        g = besselk(0, abs(imag(kp))*R) / (2*pi);
    end
end

function g = ff_dGdR(R, kp)
    if isreal(kp)
        x = kp*R;  g = -kp*(besselj(1, x) - 1i*bessely(1, x)) / (4i);
    else
        p = abs(imag(kp));  g = -p*besselk(1, p*R) / (2*pi);
    end
end

function M = ff_segint(P, A0, A1, kp, nq)
    % M(i,s) = int_{seg s} G(|P_i - r'|) dl'   (로그 특이 부분은 해석적으로)
    D = A1 - A0;  L = sqrt(sum(D.^2, 2));  T = D ./ L;
    dx = P(:,1) - A0(:,1).';  dz = P(:,2) - A0(:,2).';
    sp = dx.*T(:,1).' + dz.*T(:,2).';
    h  = abs(dx.*T(:,2).' - dz.*T(:,1).');
    u1 = -sp;  u2 = L.' - sp;
    st = -(1/(2*pi)) * (Fint(u2, h) - Fint(u1, h));
    [x, wq] = gauss_legendre(nq, 0, 1);
    M = complex(st);
    for q = 1:nq
        rp = A0 + (x(q)*L).*T;
        R  = max(hypot(P(:,1) - rp(:,1).', P(:,2) - rp(:,2).'), 1e-14);
        M  = M + (wq(q)*L.') .* (ff_G(R, kp) + log(R)/(2*pi));
    end
    function v = Fint(u, h)
        r = sqrt(u.^2 + h.^2);
        v = -u;
        nz = r > 0;   v(nz) = v(nz) + u(nz).*log(r(nz));
        hz = h > 0;   v(hz) = v(hz) + h(hz).*atan(u(hz)./h(hz));
    end
end

function [Y, Jt, Jy] = ff_solve(beta, g, a, k0, w, mu0, eps0)
    kp = ff_kap(beta, k0);  Ns = g.Ns;  L = g.L;  T = g.T;  C = g.C;  A0 = g.A0;  A1 = g.A1;
    nxt = [2:Ns, 1].';  prv = [Ns, 1:Ns-1].';
    Ld  = 0.5*L + 0.5*L(nxt);
    nodepot = @(P) (ff_segint(P, C, A1, kp, 12) + ff_segint(P, A1, C(nxt,:), kp, 12)) ./ Ld.';
    Kc = ff_segint(C, A0, A1, kp, 12);   Kn = ff_segint(A1, A0, A1, kp, 12);
    Pc = nodepot(C);                    Pn = nodepot(A1);
    Dq = zeros(Ns);  Dq(sub2ind([Ns Ns], (1:Ns).', nxt)) = 1;  Dq = Dq - eye(Ns);
    cq = -1/(1i*w);  cr = beta/w;
    phic_Jt = (1/eps0)*Pc*(cq*Dq);  phic_Jy = (1/eps0)*Kc*cr;
    phin_Jt = (1/eps0)*Pn*(cq*Dq);  phin_Jy = (1/eps0)*Kn*cr;
    TT = T*T.';
    Et_Jt = -1i*w*mu0*TT.*Kc - (phin_Jt - phin_Jt(prv,:))./L;
    Et_Jy = -(phin_Jy - phin_Jy(prv,:))./L;
    Ey_Jt = 1i*beta*phic_Jt;
    Ey_Jy = -1i*w*mu0*Kc + 1i*beta*phic_Jy;
    Z = [Et_Jt, Et_Jy; Ey_Jt, Ey_Jy];
    % 슬롯 위 자기전류 M_y = -E_ap = -1 이 만드는 E (몸체 표면, 금속 쪽)
    My = -1;  EM = complex(zeros(Ns, 1));
    off = abs(C(:,2)) > 1e-12;  Pm = C(off,:);
    fx = complex(zeros(size(Pm,1),1));  fz = fx;
    [xq, wq] = gauss_legendre(32, 0, 1);
    for s = find(g.slot).'
        for q = 1:numel(xq)
            rp = A0(s,:) + xq(q)*L(s)*T(s,:);
            dxm = Pm(:,1) - rp(1);  dzm = Pm(:,2) - rp(2);  R = hypot(dxm, dzm);
            gg = ff_dGdR(R, kp) * (wq(q)*L(s)*My);
            fx = fx + gg.*dxm./R;  fz = fz + gg.*dzm./R;
        end
    end
    EM(off)    = T(off,1).*fz + T(off,2).*(-fx);
    EM(g.slot) = 0.5*My;
    x  = Z \ [-EM; zeros(Ns,1)];
    Jt = x(1:Ns);  Jy = x(Ns+1:end);
    Y  = sum(-Jt(g.slot).*L(g.slot)) / a;          % 슬롯 평균 H_y (윗면에서 t = +x)
end

function Ex = ff_obs(g, beta, Jt, Jy, P, k0, w, mu0, eps0)
    % (x, z) 점들의 E_x  (P: N x 2, z 는 개구면 기준)
    kp = ff_kap(beta, k0);  Ns = g.Ns;  L = g.L;  T = g.T;  A0 = g.A0;  A1 = g.A1;
    nxt = [2:Ns, 1].';
    cq = -1/(1i*w);  cr = beta/w;  My = -1;
    Q  = cq*(Jt(nxt) - Jt);
    Ex = complex(zeros(size(P,1),1));
    [xq, wq] = gauss_legendre(8, 0, 1);
    sl = g.slot;
    for q = 1:numel(xq)
        rp = A0 + (xq(q)*L).*T;
        dx = P(:,1) - rp(:,1).';  dz = P(:,2) - rp(:,2).';  R = max(hypot(dx, dz), 1e-14);
        ws = (wq(q)*L).';
        gg = ff_G(R, kp);  gr = ff_dGdR(R, kp);
        Ex = Ex - 1i*w*mu0*(gg*(ws.'.*Jt.*T(:,1)));
        Ex = Ex - (1/eps0)*((gr.*dx./R)*(ws.'.*cr.*Jy));
        Ex = Ex + (gr(:,sl).*dz(:,sl)./R(:,sl))*(ws(sl).'*My);
    end
    dx = P(:,1) - A1(:,1).';  dz = P(:,2) - A1(:,2).';  R = max(hypot(dx, dz), 1e-14);
    Ex = Ex - (1/eps0)*((ff_dGdR(R, kp).*dx./R)*Q);
end

function Y = ff_Yinf(beta, a, k0, w, mu0)
    % 무한 플랜지 슬롯 (폭 a, E_x 균일) 의 같은 정의의 어드미턴스 (해석 적분)
    persistent x wx ta wa
    if isempty(x)
        [x, wx] = gauss_legendre(2000, 0, 1);  [ta, wa] = gauss_legendre(400, 0, pi/2);
    end
    q2 = k0^2 - beta^2;  s2 = @(kx) (a*ff_sinc(kx*a/(2*pi))).^2;
    if q2 > 0
        q = sqrt(q2);
        tot = 2*sum(wa.*s2(q*sin(ta)));
        S = acosh(max(600/a, 2*q)/q);
        tot = tot + 2i*sum(S*wx.*s2(q*cosh(S*x)));
    else
        p = sqrt(-q2);  S = asinh(600/a/p);
        tot = 2i*sum(S*wx.*s2(p*sinh(S*x)));
    end
    Y = (k0^2 - beta^2)/(2*pi*a*w*mu0) * tot;
end

function dH = ff_dH_matrix(FF, Na, D)
    % 개구면 펄스 점정합 연산자의 유한 플랜지 보정 (Toeplitz):  (1/pi) int dY(b) D ff_sinc(bD/2) cos(b l D) db
    b  = FF.beta;  S = D*ff_sinc(b*D/(2*pi));
    c  = (FF.wts .* (FF.Y - FF.Yinf) .* S) / pi;
    l  = (0:Na-1).' * D;
    h  = cos(l*b) * c.';
    dH = toeplitz(h);
end

function P = ff_focal_matrix(yo, Fz, FF, y_ap, D, bmax)
    % 개구면 펄스 -> (x = 0, yo, z) 의 E_x.  Fz: 해당 z 의 F(beta) (1 x Nb)
    sel = FF.beta <= bmax;  b = FF.beta(sel);
    c = FF.wts(sel) .* (D*ff_sinc(b*D/(2*pi))) .* Fz(sel) / pi;
    yo = yo(:);  ya = y_ap(:).';
    P = (cos(yo*b).*c) * cos(b.'*ya) + (sin(yo*b).*c) * sin(b.'*ya);
end

function Fz = ff_fields_z(FF, z_rel, bmax)
    % 여러 높이 (x = 0) 의 F(beta, z):  Nb_sel x Nz
    sel = find(FF.beta <= bmax);  z_rel = z_rel(:);
    Fz = complex(zeros(numel(sel), numel(z_rel)));
    for jj = 1:numel(sel)
        j = sel(jj);
        Fz(jj,:) = ff_obs(FF.geom, FF.beta(j), FF.Jt(:,j), FF.Jy(:,j), [zeros(size(z_rel)), z_rel], ...
                          FF.k0, FF.w, FF.mu0, FF.eps0).';
    end
end

function y = ff_sinc(x)
    % sin(pi x)/(pi x)  (Signal Processing Toolbox 없이)
    y = ones(size(x));  nz = x ~= 0;
    y(nz) = sin(pi*x(nz)) ./ (pi*x(nz));
end
