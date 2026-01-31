# ROCHE

Roche 楕円体（回転・自己重力・潮汐）に関する平衡形状と線形安定性の計算コードである。  
主に軸比 `(a1,a2,a3)=(1,b,c)` の楕円体について、ポテンシャル係数 A*・B* を数値積分で求め、  
Roche 系列の平衡条件を解いて `μ/(πGρ)` と `Ω^2/(πGρ)` を計算する。  
さらに、奇数・偶数モードの固有値方程式を解いて `σ^2` スペクトルを求め、安定性図を作成する。

## できること

- Roche 楕円体の平衡系列（`c` を走査して `b(c)` を解く）を計算する
- `μ/(πGρ)` と `Ω^2/(πGρ)` の系列を算出する
- 偶数/奇数モードの固有値問題から `σ^2` を算出する
- 図の生成（Fig.1/2/3 相当の曲線）を行う
- `Data/modes.dat` を生成し、臨界角度（`σ₃²=0`）を推定する

## リポジトリ構成

- `src/coeffs.jl`：楕円体ポテンシャル係数 A* の数値積分
- `src/potential.jl`：A*・Aij*・B*、体積規格化軸 `ā_i`
- `src/equilibrium.jl`：平衡条件 (⋆) を解いて Roche 系列を生成
- `src/derived.jl`：`μ/(πGρ)` と `Ω^2/(πGρ)` の計算、`ā1,ā2` 変換
- `src/stability.jl`：偶数/奇数モードの固有値計算（`σ^2`）
- `src/RocheEllipsoid.jl`：上記をまとめたモジュール
- `examples/`：図作成と安定性解析のスクリプト
- `Plots/`：生成済みの図（png/pdf）
- `Data/modes.dat`：`σ^2` スペクトル出力

## セットアップ

Julia で依存関係を導入する。

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## 使い方（例）

### 1) Roche 系列の平衡図（Fig.1/2 相当）

```sh
julia --project=. examples/make_equilibrium_figures.jl
```

出力：
- `Plots/roche_fig1_like.png`
- `Plots/roche_fig2_like.png`

### 2) 安定性解析（Fig.3 相当）とモード表

```sh
julia --project=. examples/stability_analysis.jl
```

出力：
- `Plots/fig3_sigma3sq.png`
- `Data/modes.dat`

### 3) 臨界角度（`σ₃²=0`）の推定

```sh
julia --project=. examples/stability_criteria.jl
```

`Data/modes.dat` から補間して `σ₃²=0` となる角度 `φDI` を出力する。

## 主要 API（モジュール）

```julia
include("src/RocheEllipsoid.jl")
using .RocheEllipsoid

seq = roche_sequence(p; cmin=0.01, nc=160, rtolA=5e-9)
μ, Ω2 = mu_omega_over(b, c, p)
```

### 主な関数

- `A_coeffs(b, c)`：ポテンシャル係数 A*（数値積分）
- `f_star(b, c, p)`：平衡条件 (⋆) の左辺
- `solve_b_for_c(c, p)`：与えた `c` に対して `b` を解く
- `roche_sequence(p)`：`c` を走査して Roche 系列を生成する
- `mu_omega_over(b, c, p)`：`μ/(πGρ)` と `Ω^2/(πGρ)`
- `abar12_from_seq(seq)`：体積規格化 `ā1, ā2`

## 注意
- `src/ROCHE.jl` は簡単なモジュール雛形で、実際の計算コードは `src/RocheEllipsoid.jl` にある。
