function info = check_mesh(box, sheet, M, B, G, TG, gpu_mem)
% check_mesh
%   S1 메시를 assert 로 검사하고 미지수 / 메모리 추정을 출력한다 [mm]
%   PEC 박스: 닫힌 면 (모든 edge 에 삼각형 2 개), 방향 일관 (directed edge 중복 없음),
%             Euler 수 2, 연결 성분 1, 둘러싼 부피 = PEC cell 부피 = body 부피 합, 최소 내각 > 10 도
%   대칭: PEC / 구리 꼭짓점, 기판 voxel 이 x = 0 에 대해 대칭
%   구리 시트: edge 당 삼각형 2 개 이하, 법선 +z, 연결 성분 24 (cell 마다 반쪽 2 개 + MS),
%             면적 = 구리 발자국 + MS 면적
%   기판 voxel: 부피 = sub1 부피
%   gpu_mem: GPU 메모리 [byte] (dense 행렬이 들어가는지 출력)

tol = 1e-6;
npass = 0;

% -------------------------------------------------------------------------
% PEC 박스
% -------------------------------------------------------------------------
[Eb, cb, Edb] = tri_edges(box.T);
chk(all(cb == 2), 'PEC: 닫힌 면 (edge %d 개 모두 삼각형 2 개)', numel(cb));
chk(size(unique(Edb, 'rows'), 1) == size(Edb, 1), 'PEC: 삼각형 방향 일관');
chi = size(box.P, 1) - size(Eb, 1) + size(box.T, 1);
chk(chi == 2, 'PEC: Euler 수 V - E + F = %d', chi);
chk(max(n_comp(box.T, size(box.P, 1))) == 1, 'PEC: 연결 성분 1 개');
chk(all(abs(abs(box.n) - round(abs(box.n))) < tol, 'all'), 'PEC: 삼각형 법선 축 정렬');

V_surf = sum(dot(box.P(box.T(:,1),:), cross(box.P(box.T(:,2),:), box.P(box.T(:,3),:), 2), 2)) / 6;
[dx, dy, dz] = cell_size(TG);
dV     = dx(:) .* reshape(dy, 1, []) .* reshape(dz, 1, 1, []);
V_cell = sum(dV(M.occ == 1));
V_body = sum(arrayfun(@(b) prod(diff(G.body(b).bbox, 1, 2)), B.idx.pec));
chk(abs(V_surf - V_cell) < 1e-9 * V_cell && abs(V_cell - V_body) < 1e-9 * V_body, ...
    'PEC: 바깥 법선 부피 %.6g = cell 부피 = body 부피 [mm^3]', V_surf);
n_sub_face = nnz(box.tag(:,3) == 2);
chk(n_sub_face > 0, 'PEC: 기판과 맞닿은 삼각형 %d 개', n_sub_face);
ang_b = min_angle(box.P, box.T);
chk(min(ang_b) > 10, 'PEC: 최소 내각 %.1f 도', min(ang_b));
chk(sym_ok(box.P, tol), 'PEC: 꼭짓점이 x = 0 대칭');
chk(sym_ok(sheet.P, tol), '구리: 꼭짓점이 x = 0 대칭');
chk(sym_ok((M.sub.lo + M.sub.hi) / 2, tol), '기판: voxel 이 x = 0 대칭');

% -------------------------------------------------------------------------
% 구리 시트
% -------------------------------------------------------------------------
[Es, cs, Eds] = tri_edges(sheet.T);
chk(all(cs <= 2), '구리: edge 당 삼각형 2 개 이하');
chk(size(unique(Eds, 'rows'), 1) == size(Eds, 1), '구리: 삼각형 방향 일관');
chk(all(abs(sheet.n(:,3) - 1) < tol), '구리: 법선 +z');
comp = n_comp(sheet.T, size(sheet.P, 1));
tri_comp = comp(sheet.T(:,1));
chk(max(comp(unique(sheet.T))) == 24, '구리: 연결 성분 24 개 (cell 마다 반쪽 2 개 + MS)');
ok = true;
for c = 1:24
    tg = sheet.tag(tri_comp == c, :);
    ok = ok && numel(unique(tg(tg(:,1) == 1, 2))) == 1 && nnz(unique(tg(tg(:,1) == 1, 3))) == 2 && ...
         numel(unique(tg(tg(:,1) == 2, 2))) == 1;
end
chk(ok, '구리: 성분마다 cell 1 개 (반쪽 2 개) + MS 1 개');

A_poly = 0;
for k = 1:numel(B.cell)
    for h = 1:2
        L = B.cell(k).poly{h};
        A_poly = A_poly + polyarea(L{1}(:,1), L{1}(:,2)) - sum(cellfun(@(p) polyarea(p(:,1), p(:,2)), L(2:end)));
    end
end
A_ms = sum(arrayfun(@(m) diff(m.rect(1:2)) * diff(m.rect(3:4)), B.ms));
A_cu = sum(sheet.area(sheet.tag(:,1) == 1));
A_m  = sum(sheet.area(sheet.tag(:,1) == 2));
chk(abs(A_cu - A_poly) < 1e-9 * A_poly && abs(A_m - A_ms) < 1e-9 * A_ms, ...
    '구리: 면적 %.4f + MS %.4f [mm^2] = 발자국 면적', A_cu, A_m);

% port edge: MS 안, x = 0 위의 edge
mid  = (sheet.P(Es(:,1),:) + sheet.P(Es(:,2),:)) / 2;
on0  = abs(sheet.P(Es(:,1),1)) < tol & abs(sheet.P(Es(:,2),1)) < tol & cs == 2;
n_port = zeros(1, numel(B.ms));
for m = 1:numel(B.ms)
    r = B.ms(m).rect;
    n_port(m) = nnz(on0 & mid(:,2) > r(3) & mid(:,2) < r(4));
end
chk(all(n_port >= 1), 'port: MS 마다 x = 0 edge %d 개', n_port(1));
chk(nnz(on0) == sum(n_port), 'port: x = 0 위 내부 edge 는 모두 MS 안');

% -------------------------------------------------------------------------
% 기판 voxel / 피드
% -------------------------------------------------------------------------
V_vox = sum(prod(M.sub.hi - M.sub.lo, 2));
V_sub = prod(diff(G.body(B.idx.sub).bbox, 1, 2));
chk(abs(V_vox - V_sub) < 1e-9 * V_sub, '기판: voxel %d 개, 부피 %.4f [mm^3]', size(M.sub.lo, 1), V_vox);
chk(all(M.sub.tagb(:,3) ~= 2), '기판: 경계면 바깥은 공기 또는 PEC');

fc  = [mean(B.feed.rect(1:2)), mean(B.feed.rect(3:4)), B.feed.z];
ifc = [find(TG.x < fc(1), 1, 'last'), find(TG.y < fc(2), 1, 'last'), find(TG.z < fc(3), 1, 'last')];
chk(M.occ(ifc(1), ifc(2), ifc(3)) == 0, '피드 중심이 cavity 공기 안');

fprintf('S1 메시 검사: %d 항목 통과\n', npass);

% -------------------------------------------------------------------------
% 미지수 / 메모리 추정
% -------------------------------------------------------------------------
N_box  = size(Eb, 1);
N_cu   = nnz(cs == 2);
N_roof = sum(M.sub.nface);

xb   = box.P(:,1);
N_box_s = nnz((xb(Eb(:,1)) + xb(Eb(:,2))) / 2 > tol) + nnz(abs(xb(Eb(:,1))) < tol & abs(xb(Eb(:,2))) < tol);
xs   = sheet.P(:,1);
int_s = cs == 2;
N_cu_s = nnz(int_s & (xs(Es(:,1)) + xs(Es(:,2))) / 2 > tol) + nnz(on0);
[N_roof_s] = roof_sym(M, TG, tol);

N   = N_box + N_cu + N_roof;
N_s = N_box_s + N_cu_s + N_roof_s;
fprintf('\n  미지수        전체      x = 0 대칭\n');
fprintf('  PEC RWG   %8d   %8d\n', N_box, N_box_s);
fprintf('  구리 RWG  %8d   %8d\n', N_cu, N_cu_s);
fprintf('  rooftop   %8d   %8d   (voxel %d, 면 내부 %d / 경계 %d)\n', N_roof, N_roof_s, ...
        size(M.sub.lo, 1), M.sub.nface(1), M.sub.nface(2));
fprintf('  합계      %8d   %8d\n', N, N_s);
fprintf('  dense Z (complex double)  %7.1f GB   %7.1f GB\n', 16 * N^2 / 1e9, 16 * N_s^2 / 1e9);
fprintf('  dense Z (complex single)  %7.1f GB   %7.1f GB\n', 8 * N^2 / 1e9, 8 * N_s^2 / 1e9);
if nargin > 6 && ~isempty(gpu_mem)
    fprintf('  GPU 메모리 %.1f GB: 대칭 double LU %s, single LU %s\n', gpu_mem / 1e9, ...
            yesno(16 * N_s^2 < 0.9 * gpu_mem), yesno(8 * N_s^2 < 0.9 * gpu_mem));
end
fprintf('  삼각형: PEC %d, 구리 %d (MS %d)\n', size(box.T, 1), size(sheet.T, 1), nnz(sheet.tag(:,1) == 2));
fprintf('  edge 길이 [mm]: PEC %.3f ~ %.3f, 구리 %.3f ~ %.3f\n', edge_len(box.P, Eb), edge_len(sheet.P, Es));
ang_s = sort(ang_b);
fprintf('  최소 내각 [도]: PEC %.1f (5 %% 분위 %.1f), 구리 %.1f\n', ang_s(1), ang_s(ceil(0.05 * end)), ...
        min(min_angle(sheet.P, sheet.T)));

info = struct('N_box', N_box, 'N_cu', N_cu, 'N_roof', N_roof, 'N', N, ...
              'N_box_s', N_box_s, 'N_cu_s', N_cu_s, 'N_roof_s', N_roof_s, 'N_s', N_s, 'n_port', n_port);

% =========================================================================
    function chk(cond, fmt, varargin)
        assert(cond, ['S1 검사 실패: ', fmt], varargin{:});
        npass = npass + 1;
    end
end

% =========================================================================
function [E, cnt, Ed] = tri_edges(T)
Ed = [T(:,[1 2]); T(:,[2 3]); T(:,[3 1])];
[E, ~, ie] = unique(sort(Ed, 2), 'rows');
cnt = accumarray(ie, 1);
end

function comp = n_comp(T, np)
g = graph([T(:,1); T(:,2); T(:,3)], [T(:,2); T(:,3); T(:,1)], [], np);
comp = conncomp(g);
end

function [dx, dy, dz] = cell_size(TG)
dx = diff(TG.x(:));
dy = diff(TG.y(:));
dz = diff(TG.z(:));
end

function n = roof_sym(M, TG, tol)
% x 중심 > 0 인 voxel 면 + x = 0 위의 x 법선 면 (D_x 가 대칭면을 지남)
in = M.occ == 2;
x  = TG.x(:);
xc = (x(1:end-1) + x(2:end)) / 2;
pos = reshape(xc > tol, [], 1, 1);
n = 0;
% x 법선 면: node 평면 x_i
Rp = cat(1, false(1, size(in,2), size(in,3)), in, false(1, size(in,2), size(in,3)));
F  = Rp(1:end-1,:,:) | Rp(2:end,:,:);
n  = n + nnz(F & reshape(x > tol | abs(x) < tol, [], 1, 1));
% y, z 법선 면
for a = 2:3
    pad = size(in);
    pad(a) = 1;
    Rp = cat(a, false(pad), in, false(pad));
    if a == 2
        F = Rp(:,1:end-1,:) | Rp(:,2:end,:);
    else
        F = Rp(:,:,1:end-1) | Rp(:,:,2:end);
    end
    n = n + nnz(F & pos);
end
end

function ang = min_angle(P, T)
A = P(T(:,1),:);  Bv = P(T(:,2),:);  C = P(T(:,3),:);
la = vecnorm(Bv - C, 2, 2);  lb = vecnorm(C - A, 2, 2);  lc = vecnorm(A - Bv, 2, 2);
ga = acosd(max(-1, min(1, (lb.^2 + lc.^2 - la.^2) ./ (2 * lb .* lc))));
gb = acosd(max(-1, min(1, (la.^2 + lc.^2 - lb.^2) ./ (2 * la .* lc))));
ang = min([ga, gb, 180 - ga - gb], [], 2);
end

function ok = sym_ok(P, tol)
ok = all(ismembertol(P .* [-1 1 1], P, tol, 'ByRows', true, 'DataScale', 1));
end

function v = edge_len(P, E)
L = vecnorm(P(E(:,1),:) - P(E(:,2),:), 2, 2);
v = [min(L), max(L)];
end

function s = yesno(b)
if b
    s = "가능";
else
    s = "불가";
end
end
