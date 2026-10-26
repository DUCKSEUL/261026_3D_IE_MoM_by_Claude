function S = tri_geom(S)
% tri_geom  삼각형 법선 (꼭짓점 순서 기준), 면적, 중심을 S.n, S.area, S.ctr 에 넣는다
e1 = S.P(S.T(:,2), :) - S.P(S.T(:,1), :);
e2 = S.P(S.T(:,3), :) - S.P(S.T(:,1), :);
c  = cross(e1, e2, 2);
S.area = vecnorm(c, 2, 2) / 2;
S.n    = c ./ (2 * S.area);
S.ctr  = (S.P(S.T(:,1), :) + S.P(S.T(:,2), :) + S.P(S.T(:,3), :)) / 3;
end
