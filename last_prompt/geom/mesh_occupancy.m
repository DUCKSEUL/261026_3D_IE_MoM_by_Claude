function M = mesh_occupancy(G, B, TG)
% mesh_occupancy
%   tensor grid cell 을 PEC / 기판 / 공기로 분류하고 면 사각형을 만든다 [mm]
%   - 안/밖 판정: cell 중심에서 +z 방향 ray 와 z 법선 면의 교차 수 (축 정렬 다면체)
%   - PEC 경계 = PEC cell 과 비PEC cell (격자 밖 = 공기) 사이의 면, breakpoint grid (TG.bp) 위에서
%       -> 박스 body 사이의 맞닿은 면은 자동으로 없어지고 닫힌 면이 됨
%       -> 사각형 꼭짓점 순서는 바깥 법선 (PEC -> 비PEC) 기준 반시계
%       -> 삼각형은 pec_face_mesh 가 면마다 graded 로 만듦
%   - 기판 voxel = 세분 grid 의 기판 cell, 기판 경계면 (그림/검사용)
%   - 구리 시트: 시트 grid 사각형 중심이 구리 발자국 또는 MS 안이면 선택 (법선 +z)
%
%   M.occ           세분 grid (nx-1, ny-1, nz-1) int8: 0 공기, 1 PEC, 2 기판
%   M.node.P        세분 grid 꼭짓점 (nx*ny*nz x 3),  id = i + (j-1)*nx + (k-1)*nx*ny
%   M.pec0.P/Q/tag  breakpoint grid 꼭짓점, PEC 사각형 (Nq x 4), [축, 법선 부호, 맞닿은 매질 (0 공기 / 2 기판)]
%   M.sub.ijk/lo/hi 기판 voxel,  M.sub.Qb/tagb = 기판 경계면 사각형 [축, 부호, 맞닿은 매질 (0 / 1 PEC)]
%   M.sub.nface     voxel 면 개수 [내부, 경계]
%   M.sheet.P/Q/tag 구리 시트 꼭짓점, 사각형, [종류 (1 구리, 2 MS), cell 또는 MS 번호, 반쪽 (+1 / -1 / 0)]

x = TG.x(:);  y = TG.y(:);  z = TG.z(:);
nx = numel(x); ny = numel(y); nz = numel(z);
xc = (x(1:end-1) + x(2:end)) / 2;
yc = (y(1:end-1) + y(2:end)) / 2;
zc = (z(1:end-1) + z(2:end)) / 2;
[XC, YC] = ndgrid(xc, yc);

inPEC = false(nx-1, ny-1, nz-1);
for b = B.idx.pec
    inPEC = inPEC | solid_inside(G.body(b), XC, YC, zc);
end
inSUB = solid_inside(G.body(B.idx.sub), XC, YC, zc);
assert(~any(inPEC & inSUB, 'all'), 'PEC 와 기판이 겹치는 cell 이 있습니다.');

occ = zeros(nx-1, ny-1, nz-1, 'int8');
occ(inPEC) = 1;
occ(inSUB) = 2;
M.occ = occ;

[X3, Y3, Z3] = ndgrid(x, y, z);
M.node.P  = [X3(:), Y3(:), Z3(:)];
M.node.n  = [nx, ny, nz];
nid = @(i, j, k) i + (j-1) * nx + (k-1) * nx * ny;

% -------------------------------------------------------------------------
% PEC 경계 사각형: breakpoint grid (TG.bp) 위에서 만들고 pec_face_mesh 가 면마다 다시 나눔
% -------------------------------------------------------------------------
x0 = TG.bp{1}(:);  y0 = TG.bp{2}(:);  z0 = TG.bp{3}(:);
n0 = [numel(x0), numel(y0), numel(z0)];
[XC0, YC0] = ndgrid((x0(1:end-1) + x0(2:end)) / 2, (y0(1:end-1) + y0(2:end)) / 2);
zc0 = (z0(1:end-1) + z0(2:end)) / 2;
inPEC0 = false(n0 - 1);
for b = B.idx.pec
    inPEC0 = inPEC0 | solid_inside(G.body(b), XC0, YC0, zc0);
end
inSUB0 = solid_inside(G.body(B.idx.sub), XC0, YC0, zc0);
occ0 = zeros(n0 - 1, 'int8');
occ0(inPEC0) = 1;
occ0(inSUB0) = 2;
assert(max(abs(x0 + flipud(x0))) < TG.opt.tol && isequal(occ0, flip(occ0, 1)), ...
       'PEC / 기판 배치가 x = 0 에 대해 대칭이 아닙니다.');
[X0, Y0, Z0] = ndgrid(x0, y0, z0);
nid0 = @(i, j, k) i + (j-1) * n0(1) + (k-1) * n0(1) * n0(2);
M.pec0.P = [X0(:), Y0(:), Z0(:)];
[M.pec0.Q, M.pec0.tag] = boundary_quads(inPEC0, occ0, 1, nid0, n0);

% -------------------------------------------------------------------------
% 기판 voxel
% -------------------------------------------------------------------------
[iv, jv, kv] = ind2sub(size(occ), find(inSUB));
M.sub.ijk = [iv, jv, kv];
M.sub.lo  = [x(iv), y(jv), z(kv)];
M.sub.hi  = [x(iv+1), y(jv+1), z(kv+1)];
[M.sub.Qb, M.sub.tagb] = boundary_quads(inSUB, occ, 2, nid, [nx ny nz]);
n_int = nnz(inSUB(1:end-1,:,:) & inSUB(2:end,:,:)) + ...
        nnz(inSUB(:,1:end-1,:) & inSUB(:,2:end,:)) + ...
        nnz(inSUB(:,:,1:end-1) & inSUB(:,:,2:end));
M.sub.nface = [n_int, size(M.sub.Qb, 1)];

% -------------------------------------------------------------------------
% 구리 시트 (z = TG.sz)
% -------------------------------------------------------------------------
sx = TG.sx(:);  sy = TG.sy(:);
nsx = numel(sx); nsy = numel(sy);
[SXC, SYC] = ndgrid((sx(1:end-1) + sx(2:end)) / 2, (sy(1:end-1) + sy(2:end)) / 2);
kind = zeros(size(SXC));
num  = zeros(size(SXC));
half = zeros(size(SXC));
for k = 1:numel(B.cell)
    for h = 1:2
        in = false(size(SXC));
        for l = 1:numel(B.cell(k).poly{h})
            L  = B.cell(k).poly{h}{l};
            in = xor(in, inpolygon(SXC, SYC, L(:,1), L(:,2)));
        end
        assert(~any(kind(in) ~= 0), '구리 발자국이 겹칩니다 (cell %d).', k);
        kind(in) = 1;
        num(in)  = k;
        half(in) = 3 - 2*h;
    end
end
for m = 1:numel(B.ms)
    r  = B.ms(m).rect;
    in = SXC > r(1) & SXC < r(2) & SYC > r(3) & SYC < r(4);
    assert(~any(kind(in) ~= 0), 'MS%d 가 구리와 겹칩니다.', B.ms(m).num);
    kind(in) = 2;
    num(in)  = m;
end
[is, js] = find(kind > 0);
sid = @(i, j) i + (j-1) * nsx;
M.sheet.Q   = [sid(is, js), sid(is+1, js), sid(is+1, js+1), sid(is, js+1)];
li          = sub2ind(size(kind), is, js);
M.sheet.tag = [kind(li), num(li), half(li)];
[SX, SY] = ndgrid(sx, sy);
M.sheet.P = [SX(:), SY(:), TG.sz * ones(nsx*nsy, 1)];
M.sheet.n = [nsx, nsy];
end

% =========================================================================
function in = solid_inside(bd, XC, YC, zc)
% 축 정렬 다면체: cell 중심 위쪽 (+z) 의 z 법선 면 교차 수가 홀수이면 안
cnt = zeros([size(XC), numel(zc)]);
for f = 1:numel(bd.face)
    F = bd.face(f);
    if abs(F.normal(3)) < 0.5
        continue;
    end
    zf = F.loop{1}(1,3);
    in2 = false(size(XC));
    for l = 1:numel(F.loop)
        in2 = xor(in2, inpolygon(XC, YC, F.loop{l}(:,1), F.loop{l}(:,2)));
    end
    cnt = cnt + double(in2) .* reshape(double(zc(:) < zf), 1, 1, []);
end
in = mod(cnt, 2) == 1;
end

function [Q, tag] = boundary_quads(inR, occ, self, nid, n)
% 영역 inR 과 바깥 (격자 밖 포함) 사이 면.  tag = [축, 법선 부호 (영역 -> 바깥), 바깥 cell 매질]
%   self: 바깥 매질을 셀 때 쓰는 영역 자신의 occ 값 (검사용)
Q = zeros(0, 4);
tag = zeros(0, 3);
for a = 1:3
    pad = size(inR);
    pad(a) = 1;
    Rp = cat(a, false(pad), inR, false(pad));
    Op = cat(a, zeros(pad, 'int8'), occ, zeros(pad, 'int8'));
    idx_lo = repmat({':'}, 1, 3);
    idx_hi = idx_lo;
    idx_lo{a} = 1:n(a);
    idx_hi{a} = 2:n(a)+1;
    lo = Rp(idx_lo{:});
    hi = Rp(idx_hi{:});
    olo = Op(idx_lo{:});
    ohi = Op(idx_hi{:});
    for s = [1, -1]
        if s == 1
            F = lo & ~hi;  med = ohi;   % 영역이 아래쪽 cell -> 법선 +축
        else
            F = ~lo & hi;  med = olo;
        end
        [i1, i2, i3] = ind2sub(size(F), find(F));
        md = double(med(F));
        assert(~any(md == self), '경계 면 바깥 매질 오류');
        ijk = [i1, i2, i3];        % 면이 있는 node 평면 = ijk(a)
        q = face_nodes(ijk, a, nid);
        if s == -1
            q = q(:, [1 4 3 2]);
        end
        Q   = [Q; q]; %#ok<AGROW>
        tag = [tag; [a * ones(numel(i1), 1), s * ones(numel(i1), 1), md]]; %#ok<AGROW>
    end
end
end

function q = face_nodes(ijk, a, nid)
% 축 a 에 수직인 면, node (i,j,k) 에서 시작해 +a 법선 기준 반시계
i = ijk(:,1); j = ijk(:,2); k = ijk(:,3);
switch a
    case 1   % (y, z):  e_y x e_z = +x
        q = [nid(i, j, k), nid(i, j+1, k), nid(i, j+1, k+1), nid(i, j, k+1)];
    case 2   % (z, x):  e_z x e_x = +y
        q = [nid(i, j, k), nid(i, j, k+1), nid(i+1, j, k+1), nid(i+1, j, k)];
    case 3   % (x, y):  e_x x e_y = +z
        q = [nid(i, j, k), nid(i+1, j, k), nid(i+1, j+1, k), nid(i, j+1, k)];
end
end
