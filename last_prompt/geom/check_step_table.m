function check_step_table(G, B, prm)
% check_step_table
%   STEP 분석 표 (인계 문서 2 절) 와 같은지 assert 로 검사하고 요약을 출력한다 [mm]
%   prm: lambda0, t_sub, d, Width, Lx_3D, z_ap [mm]
%
%   breakpoint 개수는 두 기준으로 검사
%     꼭짓점 기준 (메시에 쓰는 값)  : x 10, y 198, z 7
%     CARTESIAN_POINT 전체 (평면 원점 포함, 인계 문서 표의 값) : x 17, y 273, z 7

tol  = 1e-3;          % 표 값은 소수 셋째 자리
tolp = 1e-2;          % 2D 모델 상수와 비교
npass = 0;

% -------------------------------------------------------------------------
% 1. 엔티티 / body 개수
% -------------------------------------------------------------------------
kind = [G.body.kind];
chk(numel(G.body) == 79,              'body 79 개 (%d)', numel(G.body));
chk(nnz(kind == "solid") == 54,       'MANIFOLD_SOLID_BREP 54 개 (%d)', nnz(kind == "solid"));
chk(nnz(kind == "sheet") == 25,       'SHELL_BASED_SURFACE_MODEL 25 개 (%d)', nnz(kind == "sheet"));
chk(nnz([G.body.name] ~= "") == 55,   '이름 속성 55 개 (%d)', nnz([G.body.name] ~= ""));
chk(G.count.ADVANCED_FACE == 733 && G.count.PLANE == 733, 'ADVANCED_FACE = PLANE = 733 개');
chk(G.count.EDGE_CURVE == 1900 && G.count.LINE == 1900,   'EDGE_CURVE = LINE = 1900 개');

ok = true;
for b = 1:numel(G.body)
    for f = 1:numel(G.body(b).face)
        n = G.body(b).face(f).normal;
        a = find(abs(n) > 0.5);
        L = vertcat(G.body(b).face(f).loop{:});
        ok = ok && nnz(abs(abs(n) - 1) < 1e-12) == 1 && nnz(abs(n) < 1e-12) == 2 && ...
             all(abs(L(:,a) - L(1,a)) < 1e-9);
    end
end
chk(ok, '모든 면 법선이 +-x, +-y, +-z 이고 꼭짓점이 그 평면 위');

% -------------------------------------------------------------------------
% 2. breakpoint
% -------------------------------------------------------------------------
P = zeros(0, 3);
for b = 1:numel(G.body)
    for f = 1:numel(G.body(b).face)
        P = [P; vertcat(G.body(b).face(f).loop{:})]; %#ok<AGROW>
    end
end
ux = uniquetol(P(:,1), 1e-5, 'DataScale', 1)';
uy = uniquetol(P(:,2), 1e-5, 'DataScale', 1)';
uz = uniquetol(P(:,3), 1e-5, 'DataScale', 1)';
x_exp = [19.9965, 4.9965, 3.99475, 3.49475, 0.14];
x_exp = sort([-x_exp, x_exp]);
z_exp = [-39.0042, -24.0042, -12.7621, -1.52, 0, 0.018, 0.7495];
chk(numel(ux) == 10 && max(abs(ux - x_exp)) < tol, '꼭짓점 x breakpoint 10 개 %s', mat2str(ux, 6));
chk(numel(uz) == 7  && max(abs(uz - z_exp)) < tol, '꼭짓점 z breakpoint 7 개 %s', mat2str(uz, 6));
chk(numel(uy) == 198, '꼭짓점 y breakpoint 198 개 (%d)', numel(uy));

cx = uniquetol(round(G.cpt(:,1), 4), 1e-5, 'DataScale', 1)';
cy = uniquetol(round(G.cpt(:,2), 4), 1e-5, 'DataScale', 1)';
cz = uniquetol(round(G.cpt(:,3), 4), 1e-5, 'DataScale', 1)';
x17 = [19.9965, 12.4965, 4.9965, 3.9947, 3.7447, 3.4947, 1.8174, 0.14];
x17 = sort([-x17, 0, x17]);
chk(numel(cx) == 17 && max(abs(cx - x17)) < tol, 'CARTESIAN_POINT x 17 개 (표 값과 같음)');
chk(numel(cy) == 273, 'CARTESIAN_POINT y 273 개 (%d)', numel(cy));
chk(numel(cz) == 7,   'CARTESIAN_POINT z 7 개 (%d)', numel(cz));

% -------------------------------------------------------------------------
% 3. body 범위 / 재질
% -------------------------------------------------------------------------
tab = { ...
    "L",          "pec",                   [-19.997 -4.997;  -134.916  134.916; -39.004  0.749];
    "R",          "pec",                   [  4.997 19.997;  -134.916  134.916; -39.004  0.749];
    "LW",         "pec",                   [ -4.997  4.997;  -134.916 -119.916; -24.004  0.749];
    "RW",         "pec",                   [ -4.997  4.997;   119.916  134.916; -24.004  0.749];
    "gnd",        "pec",                   [ -4.997  4.997;  -134.916  134.916; -39.004 -24.004];
    "sub1",       "Rogers RT/duroid 5880", [ -4.997  4.997;  -119.916  119.916;  -1.520  0    ];
    "Rectangle1", "vacuum",                [ -4.997  4.997;    -0.150    0.150; -12.762 -12.762]};
for r = 1:size(tab, 1)
    b = find(B.name == tab{r,1});
    chk(numel(b) == 1, 'body %s 1 개', tab{r,1});
    chk(startsWith(G.body(b).material, tab{r,2}), '%s 재질 = %s', tab{r,1}, tab{r,2});
    chk(max(abs(G.body(b).bbox - tab{r,3}), [], 'all') < tol, '%s 범위 %s', tab{r,1}, mat2str(tab{r,3}));
end
chk(isequal(sort(B.name(B.idx.pec)), sort(["L", "R", "LW", "RW", "gnd"])), 'PEC 박스 = L, R, LW, RW, gnd');

% 구리 반쪽 cell
ic = B.idx.cu;
chk(numel(ic) == 48 && numel(B.cell) == 24, '구리 반쪽 48 개 = cell 24 개');
chk(nnz([G.body(ic).name] ~= "") == 24, '구리 이름 속성 24 개 (나머지 24 개는 x 대칭 짝)');
ok = true;
for b = ic
    bb = G.body(b).bbox;
    xr = sort(abs(bb(1,:)));
    ok = ok && max(abs(xr - [0.14 3.995])) < tol && max(abs(bb(3,:) - [0 0.018])) < tol && ...
         abs(diff(bb(2,:)) - 3.747) < tol && numel(G.body(b).face) == 14;
end
chk(ok, '구리: |x| in [0.14, 3.995], z in [0, 0.018], y 폭 3.747, 면 14 개');
ycell = [B.cell.yc];
chk(max(abs(diff(ycell) - prm.Width / 24)) < tol, '구리 cell y 주기 %.3f (= W/24)', prm.Width / 24);

ok = true;
for k = 1:numel(B.cell)
    for h = 1:2
        L  = B.cell(k).poly{h}{1};
        dy = diff(uniquetol(L(:,2), 1e-5, 'DataScale', 1))';
        xs = uniquetol(abs(L(:,1)), 1e-5, 'DataScale', 1)';
        ok = ok && numel(B.cell(k).poly{h}) == 1 && numel(dy) == 5 && ...
             max(abs(dy - [0.949 0.3 1.249 0.3 0.949])) < tol && max(abs(xs - [0.14 3.495 3.995])) < tol;
    end
end
chk(ok, 'E자형: finger 0.949 / slot 0.3 / finger 1.249 / slot 0.3 / finger 0.949, slot 끝 x 3.495, spine 3.995');

% lumped MS
chk(numel(B.ms) == 24, 'MS 24 개');
ok = true;
for m = 1:numel(B.ms)
    r  = B.ms(m).rect;
    ok = ok && max(abs(r(1:2) - [-0.14 0.14])) < tol && abs(diff(r(3:4)) - 0.40) < tol && abs(B.ms(m).z - 0.018) < tol;
end
chk(ok, 'MS: x in +-0.14, y 폭 0.40, z = 0.018');
chk(all(arrayfun(@(c) abs(mean(B.ms(c.ms).rect(3:4)) - c.yc) < tol, B.cell)), 'MS 가 cell 가운데 finger 중앙에 있음');
chk(isequal([B.ms.num], 1:24), 'MS 번호 = y 오름차순 cell 번호');

% -------------------------------------------------------------------------
% 4. 2D 모델 상수와 비교
% -------------------------------------------------------------------------
chk(abs(diff(B.cav.x) - prm.Lx_3D) < tolp,                  'cavity x 폭 %.3f = lambda0/3', diff(B.cav.x));
chk(abs(diff(B.cav.y) - prm.Width) < tolp,                  'cavity y 폭 %.3f = 8 lambda0', diff(B.cav.y));
chk(abs(B.cav.z_floor - (-prm.t_sub - prm.d)) < tolp,       '바닥 z %.4f = -t_sub - d', B.cav.z_floor);
chk(abs(B.cav.z_ap - prm.z_ap) < tolp,                      '개구면 z_ap %.4f = lambda0/40', B.cav.z_ap);
chk(abs(B.feed.z - (-prm.t_sub - prm.d/2)) < tolp,          '피드 z %.4f = -t_sub - d/2', B.feed.z);
chk(max(abs(B.feed.rect(1:2) - B.cav.x)) < tol,             '피드가 x 벽에서 벽까지');
chk(abs(diff(B.cav.z_sub) - prm.t_sub) < tol,               '기판 두께 %.3f', diff(B.cav.z_sub));
chk(abs(B.cav.x(2) - 3.99475 - 1.0) < 2e-3,                 '구리 - x 벽 간격 %.4f (약 1.0)', B.cav.x(2) - 3.99475);

fprintf('S0 STEP 표 검사: %d 항목 통과\n', npass);
fprintf('  cavity  x [%.4f, %.4f]  y [%.3f, %.3f]  바닥 %.4f  기판 [%.2f, %.2f]  z_ap %.4f\n', ...
        B.cav.x, B.cav.y, B.cav.z_floor, B.cav.z_sub, B.cav.z_ap);

% =========================================================================
    function chk(cond, fmt, varargin)
        assert(cond, ['S0 검사 실패: ', fmt], varargin{:});
        npass = npass + 1;
    end
end
