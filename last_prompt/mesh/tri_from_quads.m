function S = tri_from_quads(P, Q, qtag)
% tri_from_quads
%   사각형 (꼭짓점 반시계, 법선 방향 유지) 을 대각선 1-3 으로 나눠 삼각형 2 개로 만들고
%   쓰지 않는 꼭짓점을 지워 번호를 다시 매긴다
%
%   S.P (Np x 3), S.T (Nt x 3), S.tag (Nt x ?) = 사각형 tag,  S.quad (Nt x 1) = 원래 사각형 번호
%   S.n (Nt x 3) 단위 법선,  S.area (Nt x 1),  S.ctr (Nt x 3)

nq = size(Q, 1);
T  = zeros(2*nq, 3);
T(1:2:end, :) = Q(:, [1 2 3]);
T(2:2:end, :) = Q(:, [1 3 4]);

[used, ~, inv] = unique(T(:));
S.P    = P(used, :);
S.T    = reshape(inv, [], 3);
S.tag  = repelem(qtag, 2, 1);
S.quad = repelem((1:nq)', 2, 1);
S = tri_geom(S);
end
