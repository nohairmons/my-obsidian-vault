---
aliases:
  - 规范变换
---

$$ \mathbf{A}' = \mathbf{A} + \nabla \lambda $$
$$ V' = V - \frac{\partial \lambda}{\partial t} $$
#### 一、 核心概念提炼

1. **降维与冗余：** 在电动力学中，电场 $\mathbf{E}$ 和磁场 $\mathbf{B}$ 加起来有 6 个独立分量。通过引入标势 $V$（1 个分量）和矢势 $\mathbf{A}$（3 个分量），我们将未知数从 6 个降到了 4 个。
    
2. **多对一的映射（[[规范自由度]]）：** 势 $(V, \mathbf{A})$ 并不是唯一确定的。**多组不同的势可以对应完全相同的电磁场 $(\mathbf{E}, \mathbf{B})$。** 这种数学上的不唯一性被称为“[[规范自由度]]”。
    
3. [[规范变换]]： 我们可以在不改变物理场（$\mathbf{E}$ 和 $\mathbf{B}$）的前提下，对 $V$ 和 $\mathbf{A}$ 进行特定的数学变换。这种变换让我们拥有了按照实际需求去“定制”势函数的能力，从而极大地简化那些原本难看的微分方程。
    
总结结论：
$$
	\mathbf{A}^{\prime}=\mathbf{A} +\nabla \lambda
$$
$$
	\mathbf{V}^{\prime}=\mathbf{V}-\frac{ \partial \lambda }{ \partial t }
	
$$
不会改变物理场
#### 二、 逻辑与推导细节

我们要回答一个问题：**如果两组势 $(V, \mathbf{A})$ 和 $(V', \mathbf{A}')$ 给出相同的 $\mathbf{E}$ 和 $\mathbf{B}$，它们之间究竟差了什么？**

**第 1 步：定义差值**

设两组势的差值为 $\boldsymbol{\alpha}$ 和 $\beta$：

$$ \mathbf{A}' = \mathbf{A} + \boldsymbol{\alpha} $$

$$ V' = V + \beta $$

**第 2 步：利用磁场 $\mathbf{B}$ 相同的条件推导矢势关系**

我们知道磁场由矢势的旋度给出：$\mathbf{B} = \nabla \times \mathbf{A}$。

既然两组势给出的磁场相同，则：

$$ \nabla \times \mathbf{A}' = \nabla \times \mathbf{A} $$

将第一步的定义代入：

$$ \nabla \times (\mathbf{A} + \boldsymbol{\alpha}) = \nabla \times \mathbf{A} \implies \nabla \times \boldsymbol{\alpha} = 0 $$

根据矢量微积分基本定理：**一个旋度处处为零的矢量场，必然可以写成某个标量场的梯度。** 我们设这个标量场为 $\lambda(\mathbf{r}, t)$，则：

$$ \boldsymbol{\alpha} = \nabla \lambda $$

此时，矢势的变换关系已确定：$\mathbf{A}' = \mathbf{A} + \nabla \lambda$。

**第 3 步：利用电场 $\mathbf{E}$ 相同的条件推导标势关系**

电场由标势和矢势共同决定：$\mathbf{E} = -\nabla V - \frac{\partial \mathbf{A}}{\partial t}$。

两组势给出的电场相同，则：

$$ -\nabla V' - \frac{\partial \mathbf{A}'}{\partial t} = -\nabla V - \frac{\partial \mathbf{A}}{\partial t} $$

代入第一步的定义 $V' = V + \beta$ 和 $\mathbf{A}' = \mathbf{A} + \boldsymbol{\alpha}$：

$$ -\nabla (V + \beta) - \frac{\partial (\mathbf{A} + \boldsymbol{\alpha})}{\partial t} = -\nabla V - \frac{\partial \mathbf{A}}{\partial t} $$

消去相同项后得到：

$$ \nabla \beta + \frac{\partial \boldsymbol{\alpha}}{\partial t} = 0 $$

**第 4 步：联立求解，吸收时间常数**

将第 2 步得到的 $\boldsymbol{\alpha} = \nabla \lambda$ 代入上式：

$$ \nabla \beta + \frac{\partial}{\partial t}(\nabla \lambda) = 0 \implies \nabla \left( \beta + \frac{\partial \lambda}{\partial t} \right) = 0 $$

如果一个量的空间梯度（$\nabla$）为零，说明这个量在空间中各处相等，它仅仅依赖于时间。我们称这个量为 $k(t)$：

$$ \beta + \frac{\partial \lambda}{\partial t} = k(t) \implies \beta = -\frac{\partial \lambda}{\partial t} + k(t) $$

为了让公式更漂亮，我们可以把 $k(t)$ “吸收”进标量场 $\lambda$ 中。定义一个新的标量场 $\lambda_{\text{new}} = \lambda_{\text{old}} + \int k(t) dt$。

注意，对时间积分的常数项求空间梯度依然为零（$\nabla \lambda_{\text{new}} = \nabla \lambda_{\text{old}}$），所以这不会破坏我们第 2 步的推导。吸收后，常数项消失：

$$ \beta = -\frac{\partial \lambda}{\partial t} $$

**第 5 步：得出最终结论（公式 10.7）**

我们将 $\boldsymbol{\alpha}$ 和 $\beta$ 代回最开始的定义，得到了著名的[[规范变换方程]]：

$$ \mathbf{A}' = \mathbf{A} + \nabla \lambda $$

$$ V' = V - \frac{\partial \lambda}{\partial t} $$

只要 $\mathbf{A}$ 和 $V$ 按照上述规则**同步**发生改变，物理世界（$\mathbf{E}$ 和 $\mathbf{B}$）就完全感知不到这种变化。


#### 三、 具体的变换
其具体的变换形式包括[[库仑规范]]和[[洛伦兹规范]]
##### 2. [[洛伦兹规范]] ([[The Lorenz Gauge]])

**定义设定：** 令 $\nabla \cdot \mathbf{A} = -\mu_0\epsilon_0 \frac{\partial V}{\partial t}$

**推导过程：**
[[#^5a0260]]
看看关于 $\mathbf{A}$ 的未定方程，括号里那一项正是 $\nabla \cdot \mathbf{A} + \mu_0\epsilon_0 \frac{\partial V}{\partial t}$。洛伦兹规范的绝妙之处就在于**故意让括号里的这一整项等于 0**。

于是 $\mathbf{A}$ 的方程瞬间净化为：

$$ \nabla^2 \mathbf{A} - \mu_0\epsilon_0 \frac{\partial^2 \mathbf{A}}{\partial t^2} = -\mu_0 \mathbf{J} $$

同时，把这个规范条件代入关于 $V$ 的未定方程 $\nabla^2 V + \frac{\partial}{\partial t}(-\mu_0\epsilon_0 \frac{\partial V}{\partial t}) = -\frac{\rho}{\epsilon_0}$，整理得到：

$$ \nabla^2 V - \mu_0\epsilon_0 \frac{\partial^2 V}{\partial t^2} = -\frac{\rho}{\epsilon_0} $$

**数学评价：** $V$ 和 $\mathbf{A}$ 被**完全解耦**。它们满足了形式完全一致的**非齐次波动方程**。

为了体现这种数学之美，引入了**达朗贝尔算符 (d'Alembertian, $\square^2$)**：

$$ \square^2 \equiv \nabla^2 - \mu_0\epsilon_0 \frac{\partial^2}{\partial t^2} $$

方程最终被写成了极度对称的形式：

$$ \square^2 V = -\frac{\rho}{\epsilon_0} $$

$$ \square^2 \mathbf{A} = -\mu_0 \mathbf{J} $$


