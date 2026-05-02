# -*- coding: utf-8 -*-
"""
model_double.py —— 双点重建神经网络

96 → hidden → hidden → 6 的 MLP（hidden 默认 256）。
输出经 Tanh 压到 (-1, 1)，推理时乘 y_scale_mm（默认 12.5 = 晶体半边长）
还原为 mm。前 3 维 = A（ax, ay, az），后 3 维 = B（bx, by, bz）。

为什么模型本身不内嵌 permutation invariance？
  我们靠 **训练 loss** 来处理 A↔B 不可区分性（min(identity, swap)）。
  这样保留了通用 MLP 结构、便于推理后人工绑定标签。
  v1 不预测 fraction_a：信息量低，先把位置做准。
"""

from __future__ import annotations

import torch
import torch.nn as nn


class DoublePointMLP(nn.Module):
    """96 → hidden → hidden → 6，输出 Tanh ∈ (-1, 1)。"""

    def __init__(
        self,
        in_dim: int = 96,
        hidden: int = 256,
        out_dim: int = 6,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, hidden),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden, hidden),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden, out_dim),
            nn.Tanh(),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """x: (B, 96)  →  (B, 6)"""
        return self.net(x)


def count_parameters(model: nn.Module) -> int:
    return sum(p.numel() for p in model.parameters() if p.requires_grad)


if __name__ == "__main__":
    m = DoublePointMLP()
    print(m)
    print(f"\n可训练参数总数：{count_parameters(m):,}")
    x = torch.randn(8, 96)
    y = m(x)
    print(f"试跑：输入 {x.shape} → 输出 {y.shape}")
