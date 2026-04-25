# -*- coding: utf-8 -*-
"""
model.py —— 神经网络模型定义

目前提供一个 MLP（多层感知机）。CNN 留作进阶扩展。

PyTorch 模型的最小写法：
  1. 继承 torch.nn.Module
  2. 在 __init__ 里声明"有哪些层"
  3. 在 forward 里写"前向传播怎么走"
  4. 反向传播 PyTorch 会自动帮你算（autograd）

为什么最后一层要用 Tanh？
  我们把标签 y 归一化到 [-1, 1]（晶体半长）了，Tanh 的输出正好 ∈ (-1, 1)，
  网络就不用硬学"输出要落在某个区间"—— 减轻学习负担、提升稳定性。

为什么用 Dropout？
  训练阶段随机丢弃一部分神经元（例如 10%），强迫网络学到"分散的、冗余的"特征，
  避免它记住某一两个 SiPM 的特定值（过拟合）。推理时 Dropout 自动关闭。
"""

from __future__ import annotations

import torch
import torch.nn as nn


class PositionMLP(nn.Module):
    """
    96 → 256 → 256 → 3 的全连接网络。
    输入：标准化后的 96 维 SiPM 光子计数
    输出：归一化后的 3 维位置 ∈ [-1, 1]，用的时候乘回 y_scale 就是 cm。
    """

    def __init__(
        self,
        in_dim: int = 96,
        hidden: int = 256,
        out_dim: int = 3,
        dropout: float = 0.1,
    ):
        super().__init__()
        # nn.Sequential 按顺序把若干层串起来，调用时数据从上到下依次流过
        self.net = nn.Sequential(
            nn.Linear(in_dim, hidden),   # 线性层 1：96 维 → 256 维
            nn.ReLU(),                    # 非线性激活，让网络有表达非线性关系的能力
            nn.Dropout(dropout),          # 训练时随机丢弃
            nn.Linear(hidden, hidden),   # 线性层 2：256 → 256
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden, out_dim),  # 线性层 3（输出头）：256 → 3
            nn.Tanh(),                    # 把输出压到 (-1, 1)
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        x: shape (batch_size, 96)
        返回: shape (batch_size, 3)
        """
        return self.net(x)


def count_parameters(model: nn.Module) -> int:
    """工具函数：统计模型里一共有多少可训练参数。"""
    return sum(p.numel() for p in model.parameters() if p.requires_grad)


if __name__ == "__main__":
    # 直接 `python model.py` 可以快速看模型结构
    m = PositionMLP()
    print(m)
    print(f"\n可训练参数总数：{count_parameters(m):,}")
    # 造一个随机 batch 试跑前向
    x = torch.randn(8, 96)
    y = m(x)
    print(f"试跑：输入 {x.shape} → 输出 {y.shape}")
