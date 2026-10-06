from __future__ import annotations

import numpy as np


class Adam:
    def __init__(
        self, size: int, alpha: float = 1e-4, beta1: float = 0.9, beta2: float = 0.999, epsilon: float = 1e-8
    ) -> None:
        self.alpha = alpha
        self.beta1 = beta1
        self.beta2 = beta2
        self.epsilon = epsilon
        self.step = 0
        self.moment1 = np.zeros(size)
        self.moment2 = np.zeros(size)

    def update(self, weight: np.ndarray, gradient: np.ndarray) -> np.ndarray:
        self.step += 1
        self.moment1 = self.beta1 * self.moment1 + (1.0 - self.beta1) * gradient
        self.moment2 = self.beta2 * self.moment2 + (1.0 - self.beta2) * (gradient**2)
        moment1_hat = self.moment1 / (1.0 - self.beta1**self.step)
        moment2_hat = self.moment2 / (1.0 - self.beta2**self.step)
        return weight - self.alpha * moment1_hat / (np.sqrt(moment2_hat) + self.epsilon)
