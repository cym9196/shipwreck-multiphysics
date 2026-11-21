# Shipwreck Multiphysics · 沉船冲淤模拟

用 OpenFOAM 计算固定船体周围的水流、输沙和河床冲淤，用 Rust 计算悬沙颗粒对船体的冲击载荷。工程接受闭合的二进制 STL 和来流 CSV，可输出冲淤时间序列、水动力、颗粒冲击力、泥沙收支及 ParaView 回放。

## 用通用船体跑短算例

环境需要 Linux、OpenFOAM v2412、支持 C++17 的 G++、Rust 1.75 或更新版本、Python 3.10 或更新版本，以及 Bash 和 util-linux。六核示例需要 OpenMPI。

在仓库根目录运行以下命令。若 OpenFOAM 尚未加载，把示例路径换成实际安装位置。

```bash
export OPENFOAM_BASHRC=/opt/openfoam2412/etc/bashrc
bash movable_bed_baseline/scripts/build.sh
bash shipwreck_multiphysics/build.sh movable_bed_baseline
python3 shipwreck_multiphysics/run_batch.py \
  shipwreck_multiphysics/examples/impact_serial.json \
  --base movable_bed_baseline --output runs/smoke
```

`runs/smoke/batch_summary.json` 汇总运行和数值检查。加上 `--prepare-only` 可只生成 OpenFOAM 工况。`impact_batch.json` 演示角度、材料和开关对照，`impact_extended.json` 演示来流变化与减半时间步。

来流 CSV 三列为 `time_s,speed_m_s,Cs_volume_fraction`。时间从 0 开始，严格递增并覆盖计算终点。速度和体积含沙量在时间点之间线性插值。案例 JSON 使用米、秒和角度。

## 源码入口

- [shipwreck_multiphysics](shipwreck_multiphysics) 包含 C++ 工况生成器、船体掩膜、Rust 碰撞模型、求解器组装与数值检查。
- [movable_bed_baseline](movable_bed_baseline) 提供水沙方程、床面运动、边界条件库和示例字典。
- [campaigns](campaigns) 保存三种迎流姿态的工况和绘图脚本。正式工况需要对应的船体几何和采样数据。绘图依赖 NumPy、SciPy、Matplotlib、VTK 和 Pillow。

颗粒冲击采用单向耦合。水沙场驱动冲击计算，反弹结果进入载荷和能量诊断。短算例用于检查计算链路，定量预测仍需要网格收敛和材料参数标定。

## 许可证

本项目的组合求解器与项目新增代码按 [GPL-3.0-or-later](LICENSE) 发布。仓库保留了 OpenFOAM 派生文件的版权与许可声明，以及上游目录附带的许可证原文。上游来源为 [sedExnerFoam](https://github.com/SedFoam/sedExnerFoam) 和 [wallFunctions](https://github.com/SedFoam/wallFunctions)。

本仓库包含基于上述上游源码的修改版本。2026 年 11 月 22 日整理的发布包含动床求解器适配、船体掩膜和 Rust 冲击诊断，以及工况生成、构建、数值校验和回放改动。
