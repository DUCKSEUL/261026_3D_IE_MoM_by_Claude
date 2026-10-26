function G = step_read_planar(fname)
% step_read_planar
%   STEP (AP242) 텍스트를 직접 파싱해서 평면 다각형 body 목록을 만든다 (MATLAB 기본 함수만 사용)
%   면은 PLANE, 모서리는 LINE 만 허용 (곡면/곡선이 있으면 error)
%   가정: 문자열 안에 ';' 가 없음 (엔티티 개수로 검사)
%
%   G.body(b)  .id        body 엔티티 번호
%              .label     MANIFOLD_SOLID_BREP / SHELL_BASED_SURFACE_MODEL 의 이름
%              .name      AEDT_EntityName_V1   (속성이 없으면 "")
%              .material  AEDT_MaterialName_V1 (속성이 없으면 "")
%              .kind      "solid" | "sheet"
%              .face(f)   .loop{l}  꼭짓점 (Nv x 3) [mm], loop{1} = 바깥 경계
%                         .normal   PLANE 축 방향 x same_sense (1 x 3)
%              .bbox      [xmin xmax; ymin ymax; zmin zmax] [mm]
%   G.count    엔티티 종류별 개수 (struct, 복합 엔티티는 COMPLEX)
%   G.cpt      모든 CARTESIAN_POINT 좌표 (평면 원점 포함) [mm]
%   G.unit     "mm"

txt  = fileread(fname);
iD   = regexp(txt, '(?m)^DATA;', 'once');
assert(~isempty(iD), 'DATA 섹션이 없습니다: %s', fname);
data = txt(iD+5:end);
data = regexprep(data, '/\*.*?\*/', '');
data = regexprep(data, '[\r\n]', '');

tok = regexp(data, '#(\d+)\s*=\s*([^;]*);', 'tokens');
assert(numel(tok) == numel(regexp(data, '#\d+\s*=', 'start')), ...
       'STEP 엔티티 분리 실패 (문자열 안에 ; 가 있음)');
tok = vertcat(tok{:});
eid = str2double(tok(:,1));
rhs = string(tok(:,2));
N   = max(eid);

typ  = regexp(rhs, '^([A-Z0-9_]+)\s*\(', 'tokens', 'once');
isS  = ~cellfun(@isempty, typ);
TYPE = strings(N, 1);
ARG  = strings(N, 1);
TYPE(eid(isS))  = string(cellfun(@(c) c{1}, typ(isS), 'UniformOutput', false));
TYPE(eid(~isS)) = "COMPLEX";
ARG(eid(isS))   = regexprep(rhs(isS), '^[A-Z0-9_]+\s*\((.*)\)\s*$', '$1');

% 종류별 개수
[ut, ~, ic] = unique(TYPE(eid));
cnt = accumarray(ic, 1);
G.count = struct();
for t = 1:numel(ut)
    G.count.(ut(t)) = cnt(t);
end

ids_of = @(name) eid(TYPE(eid) == name);

% 곡면 / 곡선 종류 검사
surf_types = ["CYLINDRICAL_SURFACE", "CONICAL_SURFACE", "SPHERICAL_SURFACE", "TOROIDAL_SURFACE", ...
              "B_SPLINE_SURFACE_WITH_KNOTS", "SURFACE_OF_REVOLUTION", "SURFACE_OF_LINEAR_EXTRUSION"];
curv_types = ["CIRCLE", "ELLIPSE", "B_SPLINE_CURVE_WITH_KNOTS", "TRIMMED_CURVE", "POLYLINE"];
bad = intersect([surf_types, curv_types], string(ut));
assert(isempty(bad), '평면/직선이 아닌 엔티티가 있습니다: %s', strjoin(bad, ', '));

% -------------------------------------------------------------------------
% 기하 엔티티 표
% -------------------------------------------------------------------------
PT = nan(N, 3);
for id = ids_of("CARTESIAN_POINT")'
    PT(id,:) = paren_nums(ARG(id));
end
DIR = nan(N, 3);
for id = ids_of("DIRECTION")'
    DIR(id,:) = paren_nums(ARG(id));
end
VPT = zeros(N, 1);
for id = ids_of("VERTEX_POINT")'
    VPT(id) = refs(ARG(id));
end
EC = zeros(N, 3);                         % [시작 vertex, 끝 vertex, curve]
for id = ids_of("EDGE_CURVE")'
    r = refs(ARG(id));
    EC(id,:) = r(1:3);
end
assert(all(TYPE(EC(ids_of("EDGE_CURVE"), 3)) == "LINE"), 'LINE 이 아닌 EDGE_CURVE 가 있습니다.');

OE = zeros(N, 2);                         % [edge, 방향 (1 = .T.)]
for id = ids_of("ORIENTED_EDGE")'
    OE(id,:) = [refs(ARG(id)), last_bool(ARG(id))];
end
EL = cell(N, 1);
for id = ids_of("EDGE_LOOP")'
    EL{id} = refs(ARG(id));
end
FB = zeros(N, 3);                         % [loop, 방향, outer 표시]
for id = ids_of("FACE_BOUND")'
    FB(id,:) = [refs(ARG(id)), last_bool(ARG(id)), 0];
end
for id = ids_of("FACE_OUTER_BOUND")'
    FB(id,:) = [refs(ARG(id)), last_bool(ARG(id)), 1];
end
AXP = zeros(N, 3);                        % [원점, 축, 기준 방향]
for id = ids_of("AXIS2_PLACEMENT_3D")'
    r = refs(ARG(id));
    AXP(id, 1:numel(r)) = r;
end
PLN = zeros(N, 1);
for id = ids_of("PLANE")'
    PLN(id) = refs(ARG(id));
end
SH = cell(N, 1);
for id = [ids_of("CLOSED_SHELL"); ids_of("OPEN_SHELL")]'
    SH{id} = refs(ARG(id));
end

% -------------------------------------------------------------------------
% AEDT 속성 (이름, 재질)
%   GEOMETRIC_ITEM_SPECIFIC_USAGE('', 'properties', #sa, #ctx, #body)
%   PROPERTY_DEFINITION(key, '', #sa) <- PROPERTY_DEFINITION_REPRESENTATION(#pd, #rep)
%   REPRESENTATION('', (items), #ctx),  DESCRIPTIVE_REPRESENTATION_ITEM(key, value)
% -------------------------------------------------------------------------
sa2body = zeros(N, 1);
for id = ids_of("GEOMETRIC_ITEM_SPECIFIC_USAGE")'
    r = refs(ARG(id));
    sa2body(r(1)) = r(3);
end
bname = strings(N, 1);
bmat  = strings(N, 1);
for id = ids_of("PROPERTY_DEFINITION_REPRESENTATION")'
    r   = refs(ARG(id));
    pd  = r(1);
    rep = r(2);
    key = first_str(ARG(pd));
    if key ~= "AEDT_EntityName_V1" && key ~= "AEDT_MaterialName_V1"
        continue;
    end
    sa   = refs(ARG(pd));
    body = sa2body(sa(1));
    if body == 0
        continue;
    end
    items = refs(ARG(rep));
    for it = items(1:end-1)
        if TYPE(it) ~= "DESCRIPTIVE_REPRESENTATION_ITEM"
            continue;
        end
        s = all_str(ARG(it));
        if s(1) == key && key == "AEDT_EntityName_V1"
            bname(body) = s(2);
        elseif s(1) == key
            bmat(body) = s(2);
        end
    end
end

% -------------------------------------------------------------------------
% body
% -------------------------------------------------------------------------
bid  = [ids_of("MANIFOLD_SOLID_BREP"); ids_of("SHELL_BASED_SURFACE_MODEL")];
kind = [repmat("solid", numel(ids_of("MANIFOLD_SOLID_BREP")), 1); ...
        repmat("sheet", numel(ids_of("SHELL_BASED_SURFACE_MODEL")), 1)];
[bid, ord] = sort(bid);
kind = kind(ord);

G.body = struct('id', {}, 'label', {}, 'name', {}, 'material', {}, 'kind', {}, 'face', {}, 'bbox', {});
for b = 1:numel(bid)
    id     = bid(b);
    shells = refs(ARG(id));
    faces  = [SH{shells}];
    face   = struct('loop', {}, 'normal', {});
    allP   = zeros(0, 3);
    for f = 1:numel(faces)
        face(f) = read_face(faces(f));
        allP = [allP; vertcat(face(f).loop{:})]; %#ok<AGROW>
    end
    G.body(b).id       = id;
    G.body(b).label    = first_str(ARG(id));
    G.body(b).name     = bname(id);
    G.body(b).material = bmat(id);
    G.body(b).kind     = kind(b);
    G.body(b).face     = face;
    G.body(b).bbox     = [min(allP)', max(allP)'];
end
G.cpt  = PT(ids_of("CARTESIAN_POINT"), :);
G.unit = "mm";

% =========================================================================
    function F = read_face(af)
        assert(TYPE(af) == "ADVANCED_FACE", '#%d 는 ADVANCED_FACE 가 아닙니다.', af);
        r     = refs(ARG(af));
        surf  = r(end);
        bnds  = r(1:end-1);
        assert(TYPE(surf) == "PLANE", '#%d 면이 PLANE 이 아닙니다 (%s).', af, TYPE(surf));
        n = DIR(AXP(PLN(surf), 2), :);
        if ~last_bool(ARG(af))
            n = -n;
        end
        loops = cell(1, numel(bnds));
        area  = zeros(1, numel(bnds));
        outer = false(1, numel(bnds));
        for l = 1:numel(bnds)
            lp  = FB(bnds(l), 1);
            oes = EL{lp};
            v   = zeros(numel(oes), 2);
            for t = 1:numel(oes)
                e = OE(oes(t), 1);
                if OE(oes(t), 2)
                    v(t,:) = EC(e, [1 2]);
                else
                    v(t,:) = EC(e, [2 1]);
                end
            end
            assert(all(v(:,2) == circshift(v(:,1), -1)), 'EDGE_LOOP #%d 가 닫혀 있지 않습니다.', lp);
            P = PT(VPT(v(:,1)), :);
            if ~FB(bnds(l), 2)
                P = flipud(P);
            end
            loops{l} = P;
            area(l)  = abs(dot(newell(P), n)) / 2;
            outer(l) = FB(bnds(l), 3) == 1;
        end
        if ~any(outer)
            [~, io] = max(area);
        else
            io = find(outer, 1);
        end
        F.loop   = loops([io, setdiff(1:numel(bnds), io)]);
        F.normal = n;
    end
end

% =========================================================================
function r = refs(s)
r = str2double(regexp(char(s), '(?<=#)\d+', 'match'));
end

function v = paren_nums(s)
inner = regexp(char(s), '\(([^()]*)\)', 'tokens', 'once');
v = sscanf(inner{1}, '%f,')';
end

function b = last_bool(s)
m = regexp(char(s), '\.(T|F)\.', 'tokens');
b = double(m{end}{1} == 'T');
end

function s = first_str(s0)
t = regexp(char(s0), '^''([^'']*)''', 'tokens', 'once');
s = string(t{1});
end

function s = all_str(s0)
t = regexp(char(s0), '''([^'']*)''', 'tokens');
s = string(cellfun(@(c) c{1}, t, 'UniformOutput', false));
end

function n = newell(P)
Q = circshift(P, -1);
n = [sum((P(:,2) - Q(:,2)) .* (P(:,3) + Q(:,3))), ...
     sum((P(:,3) - Q(:,3)) .* (P(:,1) + Q(:,1))), ...
     sum((P(:,1) - Q(:,1)) .* (P(:,2) + Q(:,2)))];
end
