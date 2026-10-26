# 작업 인계: HFSS STEP 기반 1-layer cavity 메타표면의 3D 적분방정식 + MoM MATLAB 코드 작성

## 0. 이 세션에 함께 첨부할 파일
- `1L_cavity_step.step` : HFSS에서 export한 전체 구조 (필수)
- `main_1layer_AIE.m` : 지금 쓰고 있는 2D MoM + 개구면 적분방정식(AIE) 설계/최적화 코드 (참고용)
- `cavity3D_IEMoM_pipeline.pdf` : 이전 세션에서 정리한 3D 파이프라인 (참고용)

첨부가 없으면 먼저 요청해 줘. 답변과 코드 주석은 한국어로, 변수 이름과 코드 스타일은 `main_1layer_AIE.m`과 비슷하게 맞춰 줘.

## 1. 배경
- 10 GHz (λ₀ = 29.98 mm), 금속 cavity 안에 1층 메타표면(strip 24개)을 두고, cavity 안의 x 방향 선전류 피드로 개구면 위 초점면에 원하는 PSF를 만드는 설계야.
- 기존 설계 코드(`main_1layer_AIE.m`)는 2D y-z MoM(캐비티 안, x 방향 PEC 판 사이)과 개구면 적분방정식, 유한 플랜지 ky 스펙트럼 보정으로 바깥을 계산해. strip은 커패시턴스 C → (R, X) 테이블을 쓴 **impedance strip**이고, strip 24개의 C를 adjoint + fmincon으로 최적화해.
- 문제: HFSS에서 impedance strip으로 시뮬레이션하면 잘 맞는데, **물리 unit cell(구리 패턴 + lumped C)로 바꾸면 빔이 완전히 깨짐.**
- 목표: STEP 구조를 근사 없이 **전부 3D 적분방정식 + MoM**으로 풀고, **개구면 끝(z_ap)에서 z_f 위**의 E-field를 계산하는 MATLAB 코드를 만드는 것. 이 계산으로 impedance strip 근사가 맞는지 판정하고, 나중에는 C 최적화까지 3D로 하고 싶어.

## 2. STEP 분석 결과 (이전 세션에서 Python으로 직접 파싱해서 확인한 값)
- 형식: AP242, HOOPS Exchange, 단위 mm. 엔티티는 **PLANE 734개, LINE 1900개뿐이고 곡면/곡선이 없어.** 모든 면 법선이 ±x, ±y, ±z라서 **완전히 축 정렬된 직육면체/직선 다각형 구조**야.
- body 이름은 `DESCRIPTIVE_REPRESENTATION_ITEM('AEDT_EntityName_V1', …)`, 재질은 `'AEDT_MaterialName_V1'`에 있어. 연결 경로: `PROPERTY_DEFINITION(name,'',#shape_aspect)` → `PROPERTY_DEFINITION_REPRESENTATION` → `REPRESENTATION(items)`이고, `GEOMETRIC_ITEM_SPECIFIC_USAGE('', 'properties', #shape_aspect, #ctx, #body)`가 body를 가리켜.
- body 79개 = `MANIFOLD_SOLID_BREP` 54개 + `SHELL_BASED_SURFACE_MODEL`(시트) 25개. 이름이 붙은 것은 55개이고, 구리 반쪽 24개는 이름 속성이 없는 x 대칭 짝이야.
- 좌표 breakpoint: x 17개 `[±19.9965, ±12.4965, ±4.9965, ±3.9947, ±3.7447, ±3.4947, ±1.8174, ±0.14, 0]`, z 7개 `[-39.0042, -24.0042, -12.7621, -1.52, 0, 0.018, 0.7495]`, y 약 273개.

| Body | 종류 / 재질 | 범위 x / y / z (mm) | 역할 |
|---|---|---|---|
| L, R | solid / pec | ±[4.997, 19.997] / ±134.916 / [−39.004, 0.749] | x 방향 두꺼운 벽 (플랜지 15 mm) |
| LW, RW | solid / pec | ±4.997 / ±[119.916, 134.916] / [−24.004, 0.749] | y 끝 벽 (15 mm) |
| gnd | solid / pec | ±4.997 / ±134.916 / [−39.004, −24.004] | 바닥 (15 mm) |
| sub1 | solid / Rogers RT5880 (εr 2.2, tanδ 0.0009) | ±4.997 / ±119.916 / [−1.52, 0] | 기판. x 벽과 y 끝 벽에 모두 닿아 있음 |
| cop, cop_1…cop_23 + 이름 없는 대칭 짝 24개 | solid / copper, 두께 18 µm | x ∈ [0.14, 3.995] (짝은 [−3.995, −0.14]) / cell마다 폭 3.747 / [0, 0.018] | E자형 반쪽 cell 48개, y 주기 9.993 (= W/24) |
| MS1…MS24 | sheet / vacuum | x ∈ ±0.14 / 폭 0.40 / z = 0.018 | **HFSS lumped RLC (C)**, 가운데 finger gap에 걸침 |
| Rectangle1 | sheet / vacuum | x ∈ ±4.997 / y ∈ ±0.15 / z = −12.762 | **피드** (벽 ↔ 벽, x 방향) |

- cavity 안쪽: x 폭 9.993 (= λ₀/3 = Lx), y 폭 239.83 (= 8λ₀ = W), 바닥 z = −24.004 (= −t_sub − d, d = 0.75λ₀), 개구면 z_ap = 0.749 (= λ₀/40). 피드 z = −t_sub − d/2.
- **E자형 반쪽 cell** (예: cop, 윗면 꼭짓점): 바깥 사각형 x [0.14, 3.995] × y [−116.793, −113.046]에서 x = 0.14 쪽부터 x = 3.495까지 폭 0.3 mm slot 2개를 판 모양이야. y 방향으로 finger 0.949 / slot 0.3 / finger 1.249 / slot 0.3 / finger 0.949이고, x = 3.495–3.995가 finger를 잇는 spine이야. x 대칭 짝과 x = 0에서 0.28 mm gap으로 마주 보고, **가운데 finger gap에만 MS(C)**가 있고 위아래 finger gap은 열려 있어.
- **2D/AIE 모델과 다른 점** (빔이 깨진 원인 후보):
  1. 구리가 x 방향 벽과 **1.0 mm 떨어져 있음** (impedance strip은 벽에서 벽까지 이어짐).
  2. 단일 strip이 아니라 E자형 3-finger이고, C는 가운데에만 있음.
  3. C → (R, X) 테이블은 3점(0.051, 0.078, 0.234 pF)을 2차 피팅한 값이고, 범위 안에 공진이 있음 (X: −175 → +34 Ω).

## 3. 이미 정한 설계 결정 (이전 세션)
1. **Green 함수는 자유공간 3D 하나만** 씀 (e^{−jkR}/4πR, e^{jωt} 규약, 기존 코드의 H0^(2)와 같은 규약). 박스·기판·구리를 모두 미지수로 넣어서 AIE/플랜지 근사를 없앰.
2. **PEC 박스**(L, R, LW, RW, gnd 합집합 = 닫힌 면): RWG + **CFIE (α ≈ 0.5)**, 금속 덩어리 내부의 가짜 공진을 막기 위해.
3. **구리 cell**: 0 두께 시트(z = 0)로 근사, RWG + **EFIE**. 손실은 표면 임피던스로 선택적으로 포함.
4. **기판**: **VIE** (voxel rooftop, D-형, Zwamborn–van den Berg). 바깥 경계면(공기·PEC·구리와 닿는 면)은 half-rooftop으로 표면 전하를 표현. ε = ε₀·2.2(1 − j0.0009).
5. **Lumped C**: gap 중앙 edge 묶음을 port p로 두고, port 벡터 a_p(i) = edge 길이. 부하는 Z_L = 1/(jωC) (+ ESR).
6. **피드**: 미지수 아님. current source면 x 방향 impressed 전류 시트, lumped port면 PEC strip + delta-gap (아래 미확인 항목 1).
7. **메시**: gmsh 없이 MATLAB만으로 **tensor-grid 정합 메시**. breakpoint를 목표 h 이하로 세분 → cell을 PEC / 기판 / 공기로 분류 → PEC 경계면 = PEC cell과 비PEC cell 사이의 면 (겹치는 면이 자동으로 사라지고 닫힌 면이 됨). 구리와 MS는 z = 0 격자 사각형에서 고르고, 사각형마다 삼각형 2개로 나눔. 박스·구리·voxel 꼭짓점이 자동으로 일치함.
   - 권장 h: 박스 안쪽 + 윗면 rim 2 mm (λ/15), 박스 바깥 3.75 mm (λ/8), 구리 y 0.25 / x 0.5 mm (slot 0.3 해상), 기판 voxel 1.0 × 1.5 × 0.76 mm.
   - 미지수 추정: RWG 박스 약 25k + 구리 약 14k + rooftop 약 10k = **약 49k (dense 38 GB)**. x = 0 대칭을 쓰면 **약 25k (10 GB)**.
8. **x = 0 대칭(선택)**: 피드가 x 방향 전류이므로 x = 0은 **전기벽(PEC 대칭면)** (J_x 짝, J_y·J_z 홀). x ≥ 0만 메시하고 image 항을 성분별 부호로 더함. 대칭면 위의 gap은 반쪽 모델에서 Z_L/2.
9. **특이 적분**: 같은/이웃 삼각형 1/R, ∇(1/R)는 해석식(Graglia, Wilton). 구리가 기판 윗면 z = 0에, 기판이 벽면에 붙어 있으므로 **같은 평면 전하끼리의 결합**은 사각형 1/R 해석식으로 처리. 자기 voxel은 직육면체 해석 적분.
10. **풀이**: 1단계는 x 대칭 + GPU dense LU, 2단계는 ACA 압축 + GMRES (근거리 블록 LU 전처리).
11. **24-port 축약 (C 분리)**: (Z + A Z_L(C) Aᵀ) I = V. 사전에 I₀ = Z⁻¹V, Y = Z⁻¹A (RHS 25개)를 구해 두면, C마다 (Z_L⁻¹ + AᵀY) w = AᵀI₀, I = I₀ − Y w (24×24)만 풀면 됨. 관측은 E = G_obs I₀ − (G_obs Y) w = E₀ + H w. 기존 `Cost_Adjoint_…`의 구조(상태식 + 관측 행렬)를 그대로 재사용해서 C 24개를 3D로 최적화할 수 있음.
12. **관측**: x = 0, **z = z_ap + z_f** (개구면 끝 기준), y ∈ [−30λ₀, 30λ₀], 20001점. 세 성분을 모두 계산. z_f = 0.7317·W ≈ 175.5 mm. 기존 2D 코드는 z_f를 z = 0 기준으로 잡았으므로 옵션으로 둘 것.

## 4. 계획한 MATLAB 파일 구성
```
cav3d_main.m
geom/  step_read_planar.m, classify_bodies.m, build_tensor_grid.m, mesh_occupancy.m
mesh/  tri_from_quads.m, rwg_build.m, vie_rooftop_build.m, apply_symmetry_x.m
em/    quad_rules.m, singular_tri.m, singular_rect.m, singular_box.m,
       fill_SS_efie.m, fill_SS_mfie.m, fill_SV.m, fill_VS.m, fill_VV.m,
       near_list.m, excitation_feed.m, ports_lumped.m
solve/ solve_dense.m, solve_aca_gmres.m, port_reduce.m
post/  field_obs.m, psf_metrics.m, power_check.m
opt/   cost_ports_adjoint.m
test/  t_sphere_pec.m, t_sphere_diel.m, t_dipole.m, t_empty_box.m
```

## 5. 검증 순서
- V1 PEC 구 Mie (EFIE/MFIE/CFIE)
- V2 유전체 구·정육면체 VIE
- V3 얇은 strip 다이폴 + lumped C, Z 행렬 대칭성(상반성)
- V4 빈 박스 + 피드 → HFSS와 초점선 |E| 비교
- V5 + 기판
- V6 + 구리, 24개 모두 같은 C (3~4가지 값) → HFSS와 2D impedance-strip 모델 둘 다와 비교 (테이블 문제 판정)
- V7 최적 C 배치 + 메시 수렴
- 상시: 전력 보존 P_in = P_rad + P_ohm (1~2 % 이내)

## 6. 아직 확인 안 된 것 (사용자에게 먼저 물어볼 것)
1. Rectangle1 여기 종류: current source인가, lumped port(50 Ω 등)인가?
2. MS1–MS24 lumped RLC: 각 C 값(24개), R·L 포함 여부, 직렬/병렬 여부
3. HFSS 바깥 경계(radiation / PML, 공기 박스 크기), copper 손실 포함 여부
4. 관측 높이 정의가 z_ap + z_f가 맞는지
5. 사용 가능한 GPU VRAM / RAM (풀이 방식 선택)

## 7. 요청 (첫 단계)
위 미확인 항목을 확인한 다음 **S0–S1부터** 작성해 줘:
1. `step_read_planar.m`: STEP 텍스트를 직접 파싱해서 body별 이름, 재질, 면 다각형을 struct로 반환. MATLAB 기본 함수만 사용.
2. `classify_bodies.m`, `build_tensor_grid.m`, `mesh_occupancy.m`, `tri_from_quads.m`: tensor-grid 정합 삼각형 메시와 기판 voxel 생성.
3. 결과를 3D 그림(재질별 색)으로 그려서 HFSS 화면과 대조할 수 있게 해 주고, 표에 있는 범위와 개수가 맞는지 자동으로 검사하는 assert를 넣어 줘.

MATLAB이 없는 환경이면 Python으로 같은 로직을 먼저 검증해도 되지만, 최종 결과물은 MATLAB 코드야.
