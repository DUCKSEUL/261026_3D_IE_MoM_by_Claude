function S = pec_face_mesh(R, B, opt)
% pec_face_mesh
%   PEC 박스의 닫힌 면을 면 (region) 마다 graded 삼각형으로 나눈다 [mm]
%   입력 R = breakpoint grid 위의 PEC 경계 사각형 (mesh_occupancy 의 M.pec0)
%     R.P 꼭짓점, R.Q 사각형 (바깥 법선 기준 반시계), R.tag = [축, 법선 부호, 맞닿은 매질 (0 공기 / 2 기판)]
%
%   (1) region = 같은 평면 / 법선 / x 쪽 / 매질 / 등급 의 사각형 묶음
%       등급: 1 안쪽 (cavity 공기 또는 기판과 맞닿음) h_in,  2 윗면 rim (z = z_ap, +z) h_rim,  3 바깥 h_out
%       매질이 다르면 다른 region -> 벽의 기판 / 공기 경계선 (z = -t_sub, 0) 이 삼각형 변이 됨
%   (2) 사각형 변 (breakpoint 사이 선분) 을 양쪽 사각형 h 의 최솟값으로 한 번만 등분
%       -> 이웃 region 이 같은 점을 공유 (정합, 닫힌 면 유지)
%   (3) 사각형 안쪽 격자점 (ceil(L/h) 등분) + 변 위 점 -> region 경계를 constraint 로 constrained Delaunay
%   (4) x > 0 region 만 만들고 x = 0 에 대해 거울 복사 (꼭짓점 순서를 뒤집어 바깥 법선 유지)
%       -> 메시가 정확히 x 대칭, x = 0 위의 점은 공유
%
%   S.P, S.T, S.tag = [축, 법선 부호, 맞닿은 매질, 등급, x 쪽 (+1 / -1)], S.n, S.area, S.ctr

tol = opt.tol;
P0  = R.P;
Q   = R.Q;
tg  = R.tag;
nq  = size(Q, 1);
Nn  = size(P0, 1);
uv_of = [2 3; 3 1; 1 2];          % 축 a 평면의 (u, v): e_u x e_v = +e_a

ax   = tg(:,1);
sg   = tg(:,2);
med  = tg(:,3);
qc   = (P0(Q(:,1), :) + P0(Q(:,3), :)) / 2;
side = sign(qc(:,1));
assert(all(side ~= 0), 'x = 0 을 지나는 PEC 사각형이 있습니다 (x = 0 breakpoint 누락).');
cpl  = P0(sub2ind(size(P0), Q(:,1), ax));

% 등급 / h
probe = qc + 1e-3 * sg .* (ax == 1:3);
inner = probe(:,1) > B.cav.x(1) & probe(:,1) < B.cav.x(2) & ...
        probe(:,2) > B.cav.y(1) & probe(:,2) < B.cav.y(2) & ...
        probe(:,3) > B.cav.z_floor & probe(:,3) < B.cav.z_ap;
rim  = ax == 3 & sg == 1 & abs(cpl - B.cav.z_ap) < tol;
cls  = 3 * ones(nq, 1);
cls(rim)   = 2;
cls(inner) = 1;
hq = opt.h_out * ones(nq, 1);
hq(cls == 2) = opt.h_rim;
hq(cls == 1) = opt.h_in;

[~, ~, rid] = unique([ax, sg, round(cpl / tol), side, med, cls], 'rows');
nr = max(rid);

% -------------------------------------------------------------------------
% (2) 사각형 변 등분 (변마다 사각형 2 개: 닫힌 면)
% -------------------------------------------------------------------------
qe   = [Q(:,[1 2]); Q(:,[2 3]); Q(:,[3 4]); Q(:,[4 1])];     % (4 nq) x 2, 행 = 사각형 k 의 변 j -> (j-1) nq + k
qe_s = sort(qe, 2);
[E, ~, ie] = unique(qe_s(:,1) * (Nn + 1) + qe_s(:,2));
E    = [floor(E / (Nn + 1)), mod(E, Nn + 1)];
ecnt = accumarray(ie, 1);
assert(all(ecnt == 2), 'PEC 사각형 면이 닫혀 있지 않습니다.');
eh   = accumarray(ie, repmat(hq, 4, 1), [], @min);
ne   = size(E, 1);

PG    = P0;                         % 전체 꼭짓점 (등분점, 격자점을 뒤에 붙임)
epts  = cell(ne, 1);                % 변 위 점 id (E(e,1) -> E(e,2) 순서)
for e = 1:ne
    p1 = P0(E(e,1), :);
    p2 = P0(E(e,2), :);
    n  = max(1, ceil(norm(p2 - p1) / eh(e) - 1e-9));
    t  = (1:n-1)' / n;
    id = size(PG, 1) + (1:n-1)';
    PG = [PG; p1 + t .* (p2 - p1)]; %#ok<AGROW>
    epts{e} = [E(e,1); id; E(e,2)];
end
qedge = reshape(ie, nq, 4);         % 사각형 k 의 변 j -> 변 번호

% -------------------------------------------------------------------------
% (3) region 마다 constrained Delaunay (x > 0 쪽만)
% -------------------------------------------------------------------------
T   = zeros(0, 3);
tag = zeros(0, 5);
for r = 1:nr
    qs = find(rid == r);
    if side(qs(1)) < 0
        continue;
    end
    a  = ax(qs(1));
    s  = sg(qs(1));
    h  = hq(qs(1));
    uv = uv_of(a, :);

    % 변 위 점 + 경계 constraint
    eq   = qedge(qs, :);
    [eu, ~, je] = unique(eq(:));
    bnd  = eu(accumarray(je, 1) == 1);            % region 안에서 한 번만 쓰인 변 = 경계
    ids  = unique(vertcat(epts{eu}));
    C    = zeros(0, 2);
    for e = bnd'
        p = epts{e};
        C = [C; p(1:end-1), p(2:end)]; %#ok<AGROW>
    end

    % 사각형 안쪽 격자점
    for k = qs'
        pk = P0(Q(k,:), :);
        lo = min(pk(:, uv));
        hi = max(pk(:, uv));
        nu = max(1, ceil((hi(1) - lo(1)) / h - 1e-9));
        nv = max(1, ceil((hi(2) - lo(2)) / h - 1e-9));
        if nu < 2 || nv < 2
            continue;
        end
        [U, V] = ndgrid(lo(1) + (1:nu-1) * (hi(1) - lo(1)) / nu, lo(2) + (1:nv-1) * (hi(2) - lo(2)) / nv);
        pn = zeros(numel(U), 3);
        pn(:, a)     = cpl(k);
        pn(:, uv(1)) = U(:);
        pn(:, uv(2)) = V(:);
        ids = [ids; size(PG, 1) + (1:numel(U))']; %#ok<AGROW>
        PG  = [PG; pn]; %#ok<AGROW>
    end

    loc = zeros(size(PG, 1), 1);
    loc(ids) = 1:numel(ids);
    dt  = delaunayTriangulation(PG(ids, uv), loc(C));
    assert(size(dt.Points, 1) == numel(ids), 'region %d: Delaunay 가 점을 추가했습니다 (constraint 교차).', r);
    tl  = dt.ConnectivityList(isInterior(dt), :);
    p2  = dt.Points;
    sa  = (p2(tl(:,2),1) - p2(tl(:,1),1)) .* (p2(tl(:,3),2) - p2(tl(:,1),2)) - ...
          (p2(tl(:,3),1) - p2(tl(:,1),1)) .* (p2(tl(:,2),2) - p2(tl(:,1),2));
    tl(sa < 0, :) = tl(sa < 0, [1 3 2]);          % (u, v) 반시계 -> 법선 +a
    if s < 0
        tl = tl(:, [1 3 2]);
    end
    T   = [T; ids(tl)]; %#ok<AGROW>
    tag = [tag; repmat([a, s, med(qs(1)), cls(qs(1)), 1], size(tl, 1), 1)]; %#ok<AGROW>
end

% -------------------------------------------------------------------------
% (4) x = 0 거울 복사
% -------------------------------------------------------------------------
used = unique(T(:));
x    = PG(used, 1);
assert(all(x > -tol), 'x < 0 꼭짓점이 있습니다.');
on0  = abs(x) < tol;
PG(used(on0), 1) = 0;
mir  = zeros(size(PG, 1), 1);
mir(used(on0))  = used(on0);
mid  = used(~on0);
mir(mid) = size(PG, 1) + (1:numel(mid))';
PG   = [PG; PG(mid, :) .* [-1 1 1]];
Tm   = mir(T(:, [1 3 2]));
tagm = tag;
tagm(tag(:,1) == 1, 2) = -tagm(tag(:,1) == 1, 2);     % 거울 복사로 x 법선 부호만 바뀜
tagm(:, 5) = -1;

T   = [T; Tm];
tag = [tag; tagm];
[used, ~, inv] = unique(T(:));
S.P   = PG(used, :);
S.T   = reshape(inv, [], 3);
S.tag = tag;
S = tri_geom(S);
end
