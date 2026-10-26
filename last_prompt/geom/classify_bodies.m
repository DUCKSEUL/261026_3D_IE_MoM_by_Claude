function B = classify_bodies(G)
% classify_bodies
%   body 를 역할별로 나누고 구리 cell / lumped (MS) / 피드 정보를 만든다
%   역할  "pec"  : 박스 (L, R, LW, RW, gnd)            -> 닫힌 PEC 면 (CFIE)
%         "sub"  : 기판 (Rogers RT/duroid 5880)         -> VIE voxel
%         "cu"   : E자형 구리 반쪽 cell (48 개)         -> z = 0 시트 (EFIE)
%         "lump" : MS1..MS24 lumped RLC 시트            -> port
%         "feed" : Rectangle1 current source            -> impressed 전류 (미지수 아님)
%   이름/재질 속성이 없는 구리 짝 (x 대칭) 은 같은 label 의 body 에서 재질을 물려받고
%   name = label + "_m" 으로 둔다
%
%   B.role(b)          body 역할
%   B.idx.(role)       역할별 body 번호
%   B.cell(k)          구리 cell (y 오름차순), .body = [x>0, x<0], .poly{1|2} = z = 0 발자국 loop,
%                      .ybox, .yc, .ms (이 cell 의 MS 번호)
%   B.ms(k)            lumped 시트 (y 오름차순), .body, .name, .num (MS 번호), .rect [x0 x1 y0 y1], .z
%   B.feed             .body, .rect [x0 x1 y0 y1], .z, .dir = "x"
%   B.cav              안쪽 cavity 범위 (x, y, 바닥, 기판, 개구면 z)

tol = 1e-4;
nb  = numel(G.body);

% 재질 물려받기 (이름 없는 x 대칭 짝)
for b = 1:nb
    if G.body(b).material == ""
        src = find([G.body.label] == G.body(b).label & [G.body.material] ~= "", 1);
        assert(~isempty(src), 'body #%d (%s) 의 재질을 정할 수 없습니다.', G.body(b).id, G.body(b).label);
        G.body(b).material = G.body(src).material;
        G.body(b).name     = G.body(b).label + "_m";
    end
end
B.name     = [G.body.name];
B.material = [G.body.material];

role = strings(1, nb);
for b = 1:nb
    bd = G.body(b);
    if bd.kind == "solid" && bd.material == "pec"
        role(b) = "pec";
    elseif bd.kind == "solid" && startsWith(bd.material, "Rogers RT/duroid 5880")
        role(b) = "sub";
    elseif bd.kind == "solid" && bd.material == "copper"
        role(b) = "cu";
    elseif bd.kind == "sheet" && startsWith(bd.label, "MS")
        role(b) = "lump";
    elseif bd.kind == "sheet" && bd.label == "Rectangle1"
        role(b) = "feed";
    else
        error('분류할 수 없는 body: %s (%s, %s)', bd.label, bd.kind, bd.material);
    end
end
B.role = role;
for r = ["pec", "sub", "cu", "lump", "feed"]
    B.idx.(r) = find(role == r);
end
assert(numel(B.idx.sub) == 1 && numel(B.idx.feed) == 1, '기판/피드가 하나가 아닙니다.');

% -------------------------------------------------------------------------
% lumped 시트 (MS)
% -------------------------------------------------------------------------
ib = B.idx.lump;
yc = arrayfun(@(b) mean(G.body(b).bbox(2,:)), ib);
[~, ord] = sort(yc);
ib = ib(ord);
B.ms = struct('body', {}, 'name', {}, 'num', {}, 'rect', {}, 'z', {});
for k = 1:numel(ib)
    bb = G.body(ib(k)).bbox;
    assert(abs(bb(3,2) - bb(3,1)) < tol, '%s 가 z 평면 시트가 아닙니다.', G.body(ib(k)).label);
    B.ms(k).body = ib(k);
    B.ms(k).name = G.body(ib(k)).label;
    B.ms(k).num  = sscanf(G.body(ib(k)).label, 'MS%d');
    B.ms(k).rect = [bb(1,:), bb(2,:)];
    B.ms(k).z    = bb(3,1);
end

% -------------------------------------------------------------------------
% 구리 cell: 반쪽 2 개 (x > 0, x < 0) 를 y 범위로 묶음, z = zmin 의 -z 면 = 발자국
% -------------------------------------------------------------------------
ib  = B.idx.cu;
ycu = arrayfun(@(b) mean(G.body(b).bbox(2,:)), ib);
xcu = arrayfun(@(b) mean(G.body(b).bbox(1,:)), ib);
[~, ord] = sortrows([round(ycu(:)/tol), -sign(xcu(:))]);
ib  = ib(ord);
assert(mod(numel(ib), 2) == 0, '구리 반쪽 개수가 홀수입니다.');
B.cell = struct('body', {}, 'poly', {}, 'ybox', {}, 'yc', {}, 'ms', {});
for k = 1:numel(ib)/2
    bp = ib(2*k-1);
    bm = ib(2*k);
    assert(mean(G.body(bp).bbox(1,:)) > 0 && mean(G.body(bm).bbox(1,:)) < 0, 'cell %d 반쪽 짝이 맞지 않습니다.', k);
    assert(max(abs(G.body(bp).bbox(2,:) - G.body(bm).bbox(2,:))) < tol, 'cell %d 반쪽 y 범위가 다릅니다.', k);
    B.cell(k).body = [bp, bm];
    B.cell(k).poly = {footprint(G.body(bp)), footprint(G.body(bm))};
    B.cell(k).ybox = G.body(bp).bbox(2,:);
    B.cell(k).yc   = mean(B.cell(k).ybox);
    in_cell = arrayfun(@(m) mean(m.rect(3:4)) > B.cell(k).ybox(1) && mean(m.rect(3:4)) < B.cell(k).ybox(2), B.ms);
    assert(nnz(in_cell) == 1, 'cell %d 에 MS 가 %d 개 있습니다.', k, nnz(in_cell));
    B.cell(k).ms = find(in_cell);
end

% -------------------------------------------------------------------------
% 피드
% -------------------------------------------------------------------------
bb = G.body(B.idx.feed).bbox;
assert(abs(bb(3,2) - bb(3,1)) < tol, '피드가 z 평면 시트가 아닙니다.');
B.feed.body = B.idx.feed;
B.feed.rect = [bb(1,:), bb(2,:)];
B.feed.z    = bb(3,1);
B.feed.dir  = "x";

% -------------------------------------------------------------------------
% cavity 범위
% -------------------------------------------------------------------------
sb   = G.body(B.idx.sub).bbox;
pecb = cat(3, G.body(B.idx.pec).bbox);
gnd  = G.body(B.idx.pec(B.name(B.idx.pec) == "gnd")).bbox;
B.cav.x       = sb(1,:);
B.cav.y       = sb(2,:);
B.cav.z_floor = gnd(3,2);
B.cav.z_sub   = sb(3,:);
B.cav.z_ap    = max(pecb(3,2,:));
B.cav.z_cu    = sb(3,2);
end

% =========================================================================
function poly = footprint(bd)
% 구리 solid 의 아래면 (법선 -z, z = zmin) loop 들 (xy, 바깥 경계 먼저)
zmin = bd.bbox(3,1);
poly = {};
for f = 1:numel(bd.face)
    F = bd.face(f);
    if F.normal(3) < -0.5 && all(abs(F.loop{1}(:,3) - zmin) < 1e-6)
        poly = [poly, cellfun(@(L) L(:,1:2), F.loop, 'UniformOutput', false)]; %#ok<AGROW>
    end
end
assert(~isempty(poly), '%s 의 아래면을 찾지 못했습니다.', bd.label);
end
