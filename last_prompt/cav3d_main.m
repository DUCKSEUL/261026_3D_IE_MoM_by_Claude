% cav3d_main
%   HFSS STEP (1L_cavity.step) 전체 구조를 3D 적분방정식 + MoM 으로 푸는 주 스크립트
%
%   구조: 금속 cavity (L, R, LW, RW, gnd) + 기판 (RT5880) + E자형 구리 cell 24 쌍
%         + lumped RLC (MS1..24, 직렬 R + C) + 피드 (Rectangle1, x 방향 current source 1 A)
%   정식화 (계획): 자유공간 3D Green 함수 하나 (e^{-jkR}/4piR, e^{jwt})
%         PEC 박스 RWG + CFIE,  구리 시트 RWG + EFIE,  기판 VIE (D-형 rooftop)
%         lumped C 는 24-port 축약,  관측 z = z_ap + z_f
%
%   S0  STEP 읽기, body 분류, 표 검사
%   S1  메시 + 검사 + 그림
%         PEC 박스: 면마다 graded 삼각형 (안쪽 / rim h_in, 바깥 h_out), 변 공유로 정합, x 대칭
%         구리 시트: 구리 + MS 꼭짓점 tensor grid,  기판: tensor grid voxel
%   좌표 단위: 기하/메시 [mm] (STEP 과 같음),  전자기 상수 [SI]
clc; clear; close all;
format short e;
format compact;

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'geom'), fullfile(root, 'mesh'), fullfile(root, 'post'));

FS = 18;
Fonts = "Times New Roman";

%% =========================================================================
% 1. 파라미터
% =========================================================================
nu0      = 10e9;
omega0   = 2 * pi * nu0;
epsilon0 = 8.8541878176e-12;
mu0      = (4 * pi) * 1e-7;
c0       = 1 / sqrt(epsilon0 * mu0);
freq0    = 10e9;
lambda0  = c0 / freq0;
eta0     = sqrt(mu0 / epsilon0);
k0       = 2 * pi / lambda0;
mm       = 1e-3;

step_file   = fullfile(root, "1L_cavity.step");
save_folder = fullfile(root, "cav3d_results");

% 2D 코드 (main_1layer_AIE) 와 같은 구조 상수 -> STEP 과 대조
Width = 8.000 * lambda0;
Nm    = 24;
Dy    = Width / Nm;
Lx_3D = lambda0 / 3;
t_sub = 1.52e-3;
d     = 0.750 * lambda0;
z_ap  = lambda0 / 40;
z_f   = Width * 0.7317;           % 관측 z = z_ap + z_f (개구면 끝 기준)

% 재질 / 부하 / 피드
eps_sub   = 2.20;
tand_sub  = 0.0009;
eps_c_sub = eps_sub * (1 - 1j * tand_sub);
sigma_cu  = 5.8e7;                % HFSS copper
t_cu      = 0.018e-3;
R_lump    = 10.2;                 % MS: 직렬 R + C
C_range   = [0.051, 0.234] * 1e-12;
I0        = 1.0;                  % Rectangle1 current source [A], x 방향

% 메시 크기 [mm]
mopt.h_in        = lambda0 / 15 / mm;    % PEC 안쪽 (cavity, 기판과 맞닿은 면)
mopt.h_out       = lambda0 / 8 / mm;     % PEC 바깥
mopt.h_rim       = lambda0 / 15 / mm;    % PEC 윗면 rim (z = z_ap)
mopt.h_sub       = [1.0, 1.5, 0.76];     % 기판 voxel [x y z]
mopt.h_cu        = [0.5, 0.25];          % 구리 시트 [x y]
mopt.merge_sheet = false;                % true = 구리 선을 3D grid 에 넣어 전역 정합
mopt.tol         = 1e-5;                 % STEP 허용오차 [mm]

%% =========================================================================
% 2. 실행 환경 (GPU 필수: 행렬 채우기 / 풀이는 gpuArray)
% =========================================================================
assert(canUseGPU(), 'GPU 를 쓸 수 없습니다. 이 코드는 gpuArray 로 행렬을 채우고 풉니다.');
gdev = gpuDevice;
fprintf('GPU: %s, %.1f GB (사용 가능 %.1f GB)\n', gdev.Name, gdev.TotalMemory / 1e9, gdev.AvailableMemory / 1e9);

%% =========================================================================
% 3. S0: STEP 읽기 / 분류 / 표 검사
% =========================================================================
tic;
G = step_read_planar(step_file);
B = classify_bodies(G);
fprintf('S0 STEP 읽기: body %d 개, 면 %d 개 (%.2f s)\n', numel(G.body), G.count.ADVANCED_FACE, toc);

prm = struct('lambda0', lambda0 / mm, 't_sub', t_sub / mm, 'd', d / mm, ...
             'Width', Width / mm, 'Lx_3D', Lx_3D / mm, 'z_ap', z_ap / mm);
check_step_table(G, B, prm);

%% =========================================================================
% 4. S1: tensor-grid 메시
% =========================================================================
tic;
TG    = build_tensor_grid(G, B, mopt);
M     = mesh_occupancy(G, B, TG);
box   = pec_face_mesh(M.pec0, B, mopt);
sheet = tri_from_quads(M.sheet.P, M.sheet.Q, M.sheet.tag);
fprintf('S1 메시: grid %d x %d x %d, 시트 grid %d x %d (%.2f s)\n', ...
        numel(TG.x), numel(TG.y), numel(TG.z), numel(TG.sx), numel(TG.sy), toc);

minfo = check_mesh(box, sheet, M, B, G, TG, gdev.TotalMemory);
plot_cav3d_geom(G, B, M, TG, box, sheet, FS, Fonts);

%% =========================================================================
% 5. 저장
% =========================================================================
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end
save(fullfile(save_folder, "cav3d_S1_mesh.mat"), 'G', 'B', 'TG', 'M', 'box', 'sheet', 'mopt', 'minfo', '-v7.3');
fprintf('저장: %s\n', fullfile(save_folder, "cav3d_S1_mesh.mat"));
