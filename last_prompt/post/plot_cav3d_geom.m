function plot_cav3d_geom(G, B, M, TG, box, sheet, FS, Fonts)
% plot_cav3d_geom
%   재질별 색으로 STEP 원본과 S1 메시를 그려 HFSS 화면과 대조한다 [mm]
%   Fig 1  STEP 원본 면 (역할별 색)
%   Fig 2  S1 메시 (PEC 삼각형, 기판 경계면, 구리 시트, MS, 피드)
%   Fig 3  z = 0 구리 시트 확대 (cell 1, 2): 삼각형, MS, port edge (x = 0), 기판 voxel 선
%   Fig 4  단면: x-z (y = 0),  y-z (x = 0+, y 끝 벽 근처)

col.pec  = [0.60 0.60 0.62];
col.sub  = [0.902 0.882 0.863];
col.cu   = [0.949 0.549 0.400];
col.lump = [0.85 0.10 0.10];
col.feed = [0.10 0.35 0.90];

% -------------------------------------------------------------------------
% Fig 1: STEP 원본
% -------------------------------------------------------------------------
figure('Name', 'STEP 원본', 'Color', 'w');
hold on;
alpha_of = struct('pec', 0.15, 'sub', 0.5, 'cu', 1, 'lump', 1, 'feed', 0.8);
for b = 1:numel(G.body)
    r = B.role(b);
    for f = 1:numel(G.body(b).face)
        L = G.body(b).face(f).loop{1};
        patch(L(:,1), L(:,2), L(:,3), col.(r), 'FaceAlpha', alpha_of.(r), ...
              'EdgeColor', 0.4 * col.(r), 'EdgeAlpha', 0.5);
    end
end
dress('STEP 원본 (회색 PEC, 베이지 기판, 주황 구리, 빨강 MS, 파랑 피드)');

% -------------------------------------------------------------------------
% Fig 2: S1 메시
% -------------------------------------------------------------------------
figure('Name', 'S1 메시', 'Color', 'w');
hold on;
patch('Faces', box.T, 'Vertices', box.P, 'FaceColor', col.pec, 'FaceAlpha', 0.12, ...
      'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.15);
patch('Faces', M.sub.Qb, 'Vertices', M.node.P, 'FaceColor', col.sub, 'FaceAlpha', 0.35, ...
      'EdgeColor', [0.6 0.5 0.4], 'EdgeAlpha', 0.2);
icu = sheet.tag(:,1) == 1;
patch('Faces', sheet.T(icu,:), 'Vertices', sheet.P, 'FaceColor', col.cu, 'EdgeColor', 'none');
patch('Faces', sheet.T(~icu,:), 'Vertices', sheet.P, 'FaceColor', col.lump, 'EdgeColor', 'none');
fr = B.feed.rect;
patch(fr([1 2 2 1]), fr([3 3 4 4]), B.feed.z * [1 1 1 1], col.feed, 'EdgeColor', col.feed, 'LineWidth', 2);
dress(sprintf('S1 메시: PEC %d, 구리 %d 삼각형, 기판 voxel %d', size(box.T,1), size(sheet.T,1), size(M.sub.lo,1)));

% -------------------------------------------------------------------------
% Fig 3: 구리 시트 확대
% -------------------------------------------------------------------------
figure('Name', '구리 시트 확대', 'Color', 'w');
hold on;
yr = [B.cell(1).ybox(1) - 1, B.cell(2).ybox(2) + 1];
xg = TG.x(TG.x >= B.cav.x(1) - 1e-9 & TG.x <= B.cav.x(2) + 1e-9);
yg = TG.y(TG.y >= yr(1) & TG.y <= yr(2));
plot([xg; xg], repmat(yr', 1, numel(xg)), ':', 'Color', [0.7 0.6 0.5]);
plot(repmat(B.cav.x', 1, numel(yg)), [yg; yg], ':', 'Color', [0.7 0.6 0.5]);
sel = any(reshape(sheet.P(sheet.T, 2), [], 3) > yr(1) & reshape(sheet.P(sheet.T, 2), [], 3) < yr(2), 2);
patch('Faces', sheet.T(sel & icu,  :), 'Vertices', sheet.P(:,1:2), 'FaceColor', col.cu, 'EdgeColor', 'k', 'LineWidth', 0.3);
patch('Faces', sheet.T(sel & ~icu, :), 'Vertices', sheet.P(:,1:2), 'FaceColor', col.lump, 'EdgeColor', 'k', 'LineWidth', 0.3);
for k = 1:2
    for h = 1:2
        L = B.cell(k).poly{h}{1};
        plot(L([1:end 1], 1), L([1:end 1], 2), 'b-', 'LineWidth', 1.2);
    end
end
xline(0, 'r--', 'x = 0 (대칭면 / port)', 'LineWidth', 1.2, 'FontSize', FS - 6);
xline(B.cav.x, 'k-', 'LineWidth', 2);
axis equal;
xlim(B.cav.x + [-0.5 0.5]);
ylim(yr);
xlabel('x [mm]');
ylabel('y [mm]');
title('z = 0 구리 시트 (cell 1, 2): 파랑 = STEP 발자국, 점선 = 기판 voxel 선');
set(gca, 'FontSize', FS - 4, 'FontName', Fonts, 'Box', 'on');

% -------------------------------------------------------------------------
% Fig 4: 단면
% -------------------------------------------------------------------------
figure('Name', '단면', 'Color', 'w');
cmap = [1 1 1; col.pec; col.sub];
xc = (TG.x(1:end-1) + TG.x(2:end)) / 2;
yc = (TG.y(1:end-1) + TG.y(2:end)) / 2;
[~, j0] = min(abs(yc));
[~, i0] = min(abs(xc - 1e-3) + 1e9 * (xc < 0));
subplot(1, 2, 1);
draw_slice(TG.x, TG.z, squeeze(M.occ(:, j0, :)), cmap);
hold on;
plot(B.feed.rect(1:2), B.feed.z * [1 1], '-', 'Color', col.feed, 'LineWidth', 3);
plot([0.14 3.995; -0.14 -3.995]', B.cav.z_cu * [1 1; 1 1]', '-', 'Color', col.cu, 'LineWidth', 3);
xlabel('x [mm]');
ylabel('z [mm]');
title(sprintf('x-z 단면 (y = %.2f)', yc(j0)));
set(gca, 'FontSize', FS - 4, 'FontName', Fonts);
subplot(1, 2, 2);
draw_slice(TG.y, TG.z, squeeze(M.occ(i0, :, :)), cmap);
hold on;
for k = 1:numel(B.cell)
    plot(B.cell(k).ybox, B.cav.z_cu * [1 1], '-', 'Color', col.cu, 'LineWidth', 3);
end
xlim([B.cav.y(2) - 25, max(TG.y)]);
xlabel('y [mm]');
ylabel('z [mm]');
title(sprintf('y-z 단면 (x = %.2f), y 끝 벽 근처', xc(i0)));
set(gca, 'FontSize', FS - 4, 'FontName', Fonts);

% =========================================================================
    function dress(ttl)
        axis equal tight;
        view(35, 25);
        xlabel('x [mm]');
        ylabel('y [mm]');
        zlabel('z [mm]');
        title(ttl);
        set(gca, 'FontSize', FS - 4, 'FontName', Fonts);
        grid on;
        camlight headlight;
        lighting none;
    end
end

% =========================================================================
function draw_slice(u, v, O, cmap)
% 비등간격 grid 단면: cell 마다 사각형
u = u(:);
v = v(:);
[I, J] = ndgrid(1:numel(u)-1, 1:numel(v)-1);
I = I(:);
J = J(:);
nf = numel(I);
X = [u(I), u(I+1), u(I+1), u(I)]';
Y = [v(J), v(J), v(J+1), v(J+1)]';
patch('Faces', reshape(1:4*nf, 4, nf)', 'Vertices', [X(:), Y(:)], ...
      'FaceVertexCData', cmap(double(O(:)) + 1, :), 'FaceColor', 'flat', ...
      'EdgeColor', [0.85 0.85 0.85], 'LineWidth', 0.2);
axis equal tight;
set(gca, 'Box', 'on');
end
