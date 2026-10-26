function TG = build_tensor_grid(G, B, opt)
% build_tensor_grid
%   3D tensor grid (x, y, z: 기판 voxel, cell 분류) 와 구리 시트용 2D grid (sx, sy) 를 만든다 [mm]
%   PEC 면 삼각형은 이 grid 가 아니라 breakpoint grid (TG.bp) + pec_face_mesh 로 만든다
%   (1) breakpoint = PEC + 기판 꼭짓점 좌표 (+ x = 0 대칭면)
%   (2) 이웃 breakpoint 사이 구간을 목표 h 이하로 등분
%       구간의 h = min(opt.h_out, 구간 중점을 포함하는 영역들의 h)
%         영역: cavity 안 (h_in), 기판 (h_sub), 윗면 rim (z = z_ap 면, x·y 에만 h_rim)
%   (3) 구리 시트 grid: 구리 + MS 꼭짓점, cell 안 y 구간은 h_cu(2), 구리 x 구간은 h_cu(1)
%       opt.merge_sheet = true 이면 구리 breakpoint/영역을 3D grid 에 넣고 시트 grid 를 3D grid 에서 잘라 씀
%       (전역 정합, 대신 구리의 가는 선이 박스 면 전체로 퍼짐)
%   피드 (Rectangle1) 는 미지수가 아니므로 breakpoint 에 넣지 않음
%
%   opt: h_in, h_out, h_rim, h_sub [hx hy hz], h_cu [hx hy], merge_sheet, tol

tol = opt.tol;

% -------------------------------------------------------------------------
% breakpoint
% -------------------------------------------------------------------------
P3 = body_points(G, [B.idx.pec, B.idx.sub]);
bp = cell(1, 3);
for a = 1:3
    bp{a} = uniquetol(P3(:,a), tol, 'DataScale', 1)';
end
bp{1} = uniquetol([bp{1}, 0], tol, 'DataScale', 1);

Ps  = body_points(G, [B.idx.cu, B.idx.lump]);
sbp = {uniquetol([Ps(:,1); 0], tol, 'DataScale', 1)', uniquetol(Ps(:,2), tol, 'DataScale', 1)'};

% -------------------------------------------------------------------------
% 영역별 목표 h   (lo, hi: 1 x 3,  h: 1 x 3)
% -------------------------------------------------------------------------
sb   = G.body(B.idx.sub).bbox;
allb = [min(P3)', max(P3)'];
reg  = struct('lo', {}, 'hi', {}, 'h', {});
reg(end+1) = struct('lo', [B.cav.x(1), B.cav.y(1), B.cav.z_floor], ...
                    'hi', [B.cav.x(2), B.cav.y(2), B.cav.z_ap], 'h', opt.h_in * [1 1 1]);
reg(end+1) = struct('lo', sb(:,1)', 'hi', sb(:,2)', 'h', opt.h_sub);
reg(end+1) = struct('lo', [allb(1:2,1)', B.cav.z_ap], 'hi', [allb(1:2,2)', B.cav.z_ap], ...
                    'h', [opt.h_rim, opt.h_rim, inf]);

cu_reg = struct('lo', {}, 'hi', {}, 'h', {});
for k = 1:numel(B.cell)
    xb = [min(cellfun(@(p) min(p(:,1)), [B.cell(k).poly{:}])), max(cellfun(@(p) max(p(:,1)), [B.cell(k).poly{:}]))];
    cu_reg(end+1) = struct('lo', [xb(1), B.cell(k).ybox(1), B.cav.z_cu], ...
                           'hi', [xb(2), B.cell(k).ybox(2), B.cav.z_cu], 'h', [opt.h_cu, inf]); %#ok<AGROW>
end

if opt.merge_sheet
    for a = 1:2
        bp{a} = uniquetol([bp{a}, sbp{a}], tol, 'DataScale', 1);
    end
    reg = [reg, cu_reg];
end

TG.x = refine_axis(bp{1}, @(m) h_at(reg, 1, m, opt.h_out));
TG.y = refine_axis(bp{2}, @(m) h_at(reg, 2, m, opt.h_out));
TG.z = refine_axis(bp{3}, @(m) h_at(reg, 3, m, opt.h_out));
TG.bp = bp;

% -------------------------------------------------------------------------
% 구리 시트 grid (z = 기판 윗면)
% -------------------------------------------------------------------------
xs = [min(Ps(:,1)), max(Ps(:,1))];
ys = [min(Ps(:,2)), max(Ps(:,2))];
if opt.merge_sheet
    TG.sx = TG.x(TG.x >= xs(1) - tol & TG.x <= xs(2) + tol);
    TG.sy = TG.y(TG.y >= ys(1) - tol & TG.y <= ys(2) + tol);
else
    TG.sx = refine_axis(sbp{1}, @(m) h_at(cu_reg, 1, m, inf));
    TG.sy = refine_axis(sbp{2}, @(m) h_at(cu_reg, 2, m, inf));
end
TG.sz  = B.cav.z_cu;
TG.sbp = sbp;
TG.opt = opt;
end

% =========================================================================
function P = body_points(G, ib)
P = zeros(0, 3);
for b = ib
    for f = 1:numel(G.body(b).face)
        P = [P; vertcat(G.body(b).face(f).loop{:})]; %#ok<AGROW>
    end
end
end

function h = h_at(reg, a, m, h_out)
h = h_out;
for r = 1:numel(reg)
    if m >= reg(r).lo(a) && m <= reg(r).hi(a)
        h = min(h, reg(r).h(a));
    end
end
end

function g = refine_axis(bp, hfun)
% 이웃 breakpoint 사이를 ceil(L/h) 등분 (breakpoint 는 정확히 유지)
g = bp(1);
for i = 1:numel(bp)-1
    L = bp(i+1) - bp(i);
    h = hfun((bp(i) + bp(i+1)) / 2);
    n = max(1, ceil(L / h - 1e-9));
    g = [g, bp(i) + (1:n-1) * L / n, bp(i+1)]; %#ok<AGROW>
end
end
