## 1 古典格式

设$T>0$是给定的终止时刻，考虑一维热传导方程的第一边值问题(HD)
$$
\begin{cases}
u_t=ax_{xx}+f(x,t),&(x,t)\in(0,1)\times(0,T] \\
u(x,0)=u_0(x),& x\in[0,1] \\
u(0,t)=\phi_0(t),u(1,t)=\phi_1(t),&t\in(0,T]
\end{cases} \qquad (HD)
$$
### 1.1 格式构造

#### 1.1.1 计算区域的离散



记网比$\mu=\dfrac{\Delta t}{(\Delta x)^2}$



#### 1.1.2 微分方程的离散

设真解$[u]$足够光滑，将时间导数离散为一阶向前差商，空间导数离散为二阶中心差商
$$
[u_t]_j^n=\dfrac{[u]_j^{n+1}-[u]_j^n}{\Delta t}+O(\Delta t)=\dfrac{\Delta_t[u]_j^n}{\Delta t}+O(\Delta t) \\
[u_{xx}]_j^n=\dfrac{[u]^n_{j+1}-2[u]^n_j+[u]^n_{j-1}}{(\Delta x)^2}+O((\Delta x)^2)=\dfrac{\delta_x^2[u]_j^n}{(\Delta x)^2}+O((\Delta x)^2)
$$
在网格点$(x_j,t^n)$处微分方程精确成立
$$
\dfrac{\Delta_t[u]_j^n}{\Delta t}-a\dfrac{\delta_x^2[u]_j^n}{(\Delta x)^2}=f_j^n+O((\Delta x)^2+\Delta t)
$$
截去小量并用数值解$u$替换
$$
\Delta_tu_j^n=\mu a\delta_x^2u_j^n+\Delta tf_j^n
$$
称为显式离散

将时间导数离散为一阶向后差商，空间导数离散为二阶中心差商
$$
[u_t]_j^{n+1}=\dfrac{[u]_j^{n+1}-[u]_j^n}{\Delta t}+O(\Delta t)=\dfrac{\nabla_t[u]_j^{n+1}}{\Delta t}+O(\Delta t) \\
[u_{xx}]_j^{n+1}=\dfrac{[u]^{n+1}_{j+1}-2[u]^{n+1}_j+[u]^{n+1}_{j-1}}{(\Delta x)^2}+O((\Delta x)^2)=\dfrac{\delta_x^2[u]_j^{n+1}}{(\Delta x)^2}+O((\Delta x)^2)
$$
类似可得
$$
\Delta_tu_j^n=\mu a\delta_x^2u_j^{n+1}+\Delta tf_j^{n+1}
$$
称为隐式离散

#### 1.1.3 定解条件的离散

直接赋值即可

### 1.2 古典差分格式

全显格式=显式离散+数值定解条件

全隐格式=隐式离散+数值定解条件

全显格式可按照时间推进解出各格点处的函数值，全隐格式需将多个差分方程联立解出



在每一时间层内定义$u^n=[u^n_1,\dots,u^n_{J-1}]^\mathrm{T},n=0:N$

离散格式均按时间推进，单步推进过程对应线性方程组族
$$
B_1 u^{n+1} = B_0 u^n + \Delta t F^n,\quad n=0:N-1
$$
称其为差分方法的标准格式

对于全显格式
$$
B_1=I \\
B_0=\mathrm{tridiag}\{\mu a,1-2\mu a,\mu a\}=
\begin{bmatrix}
1-2\mu a & \mu a & \quad \\
\mu a & 1-2\mu a & \mu a \\
\quad & \ddots & \ddots & \ddots & \quad \\
\quad & \quad & \mu a & 1-2\mu a & \mu a \\
\quad & \quad & \quad & \mu a & 1-2\mu a \\
\end{bmatrix} 
\\
F^n=[f_1^n+\dfrac{\mu a}{\Delta t}\phi_0^n,f_2^n,\dots,f_{J-2}^n,f_{J-1}^n+\dfrac{\mu a}{\Delta t}\phi_1^n]^\mathrm{T}
$$
对于全隐格式
$$
B_1=\mathrm{tridiag}\{-\mu a,1+2\mu a,-\mu a\}=
\begin{bmatrix}
1+2\mu a & -\mu a & \quad \\
-\mu a & 1+2\mu a & -\mu a \\
\quad & \ddots & \ddots & \ddots & \quad \\
\quad & \quad & -\mu a & 1+2\mu a & -\mu a \\
\quad & \quad & \quad & -\mu a & 1+2\mu a \\
\end{bmatrix}
\\
B_0=I \\
F^n=[f_1^n+\dfrac{\mu a}{\Delta t}\phi_0^n,f_2^n,\dots,f_{J-2}^n,f_{J-1}^n+\dfrac{\mu a}{\Delta t}\phi_1^n]^\mathrm{T}
$$

### 1.3 可行性和计算效率





### 1.4 不同定解条件下的模型推广

设$T>0$是给定的终止时刻，考虑一维热传导方程的无界初值问题(HI)
$$
\begin{cases}
u_t=au_{xx}+f(x,t),&(x,t)\in\mathbb{R}\times(0,T]\\
u(x,0)=u_0(x),&x\in\mathbb{R}
\end{cases}\qquad(HI)
$$






考虑周期边值问题(HP)
$$
\begin{cases}
u_t=au_{xx}+f(x,t),&(x,t)\in\mathbb{R}\times(0,T]\\
u(x,t)=u(x+1,t),&(x,t)\in\mathbb{R}\times[0,T]
\end{cases} \qquad (HP)
$$
全隐格式可设计为
$$
\Delta_t u_j^n = \mu a \delta_x^2 u_j^{n+1} + \Delta t f_j^{n+1},\quad j=1:J,\quad n=0:N-1
$$
初值条件$u_j^0=u_0(x_j)$，周期边值条件$u^n_0=u_J^n,u^n_{J+1}=u^n_1,n=0:N$

单步时间推进
$$
B_1u^{n+1}=B_0u^n+\Delta tF^n,n=0:N
$$

$$
B_1=\mathrm{ptridiag}\{-\mu a,1+2\mu a,-\mu a\}=
\begin{bmatrix}
1+2\mu a & -\mu a & \quad & \quad & -\mu a \\
-\mu a & 1+2\mu a & -\mu a \\
\quad & \ddots & \ddots & \ddots & \quad \\
\quad & \quad & -\mu a & 1+2\mu a & -\mu a \\
-\mu a & \quad & \quad & -\mu a & 1+2\mu a \\
\end{bmatrix}
\\
B_0=I \\
F^n=[f_1^{n+1},\dots,f_J^{n+1}]^\mathrm{T}
$$


> （Sherman-Morrison）若$A\in\mathrm{GL}_n(\mathbb{R})$，$w,z\in\mathbb{R}^n$，那么$A-wz^\mathrm{T}$可逆且
> $$
> \left(A-wz^\mathrm{T}\right)^{-1}=A^{-1}+\left(1-z^\mathrm{T}A^{-1}w\right)^{-1}A^{-1}wz^\mathrm{T}A^{-1}
> $$

证明：
$$
\begin{align}
& \left( A - w z^\mathrm{T} \right) \left( A^{-1} + \left( 1 - z^\mathrm{T} A^{-1} w \right)^{-1} A^{-1} w z^\mathrm{T} A^{-1} \right) \\
=& I + \left( 1 - z^\mathrm{T} A^{-1} w \right)^{-1} w z^\mathrm{T} A^{-1} - w z^\mathrm{T} A^{-1} - \left( 1 - z^\mathrm{T} A^{-1} w \right)^{-1} w z^\mathrm{T} A^{-1} w z^\mathrm{T} A^{-1} \\
=& I + \left( 1 - z^\mathrm{T} A^{-1} w \right)^{-1} w \left( 1 - z^\mathrm{T} A^{-1} w \right) z^\mathrm{T} A^{-1} - w z^\mathrm{T} A^{-1} \\
=& I+w z^\mathrm{T} A^{-1} - w z^\mathrm{T} A^{-1} \\
=& I \qquad \mathbf{Q.E.D.}
\end{align}
$$
对于单步时间推进，注意到分解
$$
\begin{align}
	B_1 =
	\begin{bmatrix}
		1 + 3 \mu a & - \mu a & & & \\
		- \mu a & 1 + 2 \mu a & \ddots & \\
         & \ddots & \ddots & - \mu a \\
         &&-\mu a & 1 + 3 \mu a
	\end{bmatrix}
	-
	\begin{bmatrix}
		\mu a \\
		0 \\
         \vdots \\
         0 \\
         \mu a
	\end{bmatrix}
	\begin{bmatrix}
		1 & 0 & \cdots & 0 & 1
	\end{bmatrix}
	= T - wz^\mathrm{T}
\end{align}
$$
据Sherman-Morrison公式，
$$
\begin{align}
u^{n+1} &= B_1^{-1} b \\
&= (T - w z^\mathrm{T})^\mathrm{-1} b \\
&= T^{-1} b + (1 - z^\mathrm{T} T^{-1} w)^{-1} T^{-1} w z^\mathrm{T} T^{-1} b
\end{align}
$$

## 2 线性差分格式基本理论

网格点的差分方程可以汇总为差分格式
$$
B_1u^{n+1}=B_0u^n+\Delta tG^n
$$
基于可行性要求，默认$B_1$可逆，定义规范形式
$$
u^{n+1} = B u^n + \Delta t H^n
$$
其中$B=B^{-1}_1B_0,H^n=B^{-1}_1G^n$。差分格式相同当且仅当它们的规范形式是相同的。

### 2.1 相容性

2.1.1 逐点相容性

设离散对象为线性偏微分方程$L[u]=g$，其中$L$为线性微分算子，$[u]$为充分光滑的真解。

> 若差分方程在网格点$(x_j,t^n)$处具有局部描述
> $$
> L_{\Delta x,\Delta t} u^n_j = g
> $$
> 则相应的局部截断误差为
> $$
> \tau_j^n = L_{\Delta x,\Delta t}[u]_j^n - g_j^n
> $$
> 若存在时空网络加密路径$\Delta x = G(\Delta t)$，其中$G(0)=0$，使得
> $$
> \tau_j^n \rightarrow 0 \quad (\Delta x \rightarrow 0, \Delta t \rightarrow 0)
> $$
> 则称差分方程逐点相容于偏微分方程。若存在不可改善的正常数$m_1,m_2$使得
> $$
> \tau_j^n = O((\Delta x)^{m_1} + (\Delta t)^{m_2})
> $$
> 则称差分方程的局部截断误差阶为$(m_1,m_2)$。

计算时要保证差分方程同离散对象具有相同的物理量纲，计算出的误差阶才正确

> 全显格式具有(2,1)阶局部截断误差。

> 全隐格式具有(2,1)阶局部截断误差。

局部截断误差阶与推导过程无关，即在任何位置执行Taylor展开计算得到的误差阶都是相同的。

2.1.2 整体相容性

考虑差分方程$u^{n+1}=Bu^n+\Delta tH^n$。

> 若存在网格函数$\varPhi^n = \{\varPhi_j^n\}$使得
> $$
> [u]^{n+1} = B [u]^n + \Delta t H^n + \Delta t \varPhi^n
> $$
> 则称$\varPhi^n$​是差分格式的局部截断误差（网格函数）。
>
> 设$\left\| \cdot \right\|$是给定的离散范数，当$\Delta x, \Delta t \rightarrow 0$时，若
> $$
> \left\| \varPhi^n \right\| = \left\| \varPhi^n \right\|_{\Delta x} \rightarrow 0, \forall n
> $$
> 则称差分格式按$\left\| \cdot \right\|$​模整体相容于偏微分方程定解问题。若存在不可改善的正常数$m_1,m_2$使得
> $$
> \left\| \varPhi^n \right\| = O((\Delta x)^{m_1} + (\Delta t)^{m_2})
> $$
> 则称差分方程的$\left\| \cdot \right\|$模相容阶为$(m_1,m_2)$。

> 建立模型问题(HD),(HP),(HI)的古典格式在最大模度量和$L_2$模度量下的整体相容性结果。

对于全隐格式，在(HD)下，



2.1.3 导数的差商离散

待定系数法

给定等距离散模板$\{x_{j+s}\}_{s=-l:r}$和离散焦点$x^* = x_j + \theta \Delta x \quad ( \theta \in [0,1))$，利用待定系数方法确定$D^m p (x^*)$的差商公式$D^m p(x^*) = \sum_{s=-l}^r \alpha_s p(x_{j+s}) + O((\Delta x)^\sigma)$

据Taylor级数
$$
D^m p(x^*) = \sum_{k=0}^{\infty} \beta_k (\Delta x)^k D^k p(x^*) + O((\Delta x)^\sigma)
$$
其中
$$
\beta_k = \sum_{s=-l}^r \alpha_s \dfrac{(s - \theta)^k}{k!}
$$
比较系数可得
$$
\beta_0=\beta_1=\dots=\beta_{m-1}=0,\quad \beta_m=\dfrac{1}{(\Delta x)^m}
$$
据此解出$\{\alpha_s\}_{s=-l}^r$，如果方程组无解，表明需要加宽离散模板。

再计算$\beta_{m+1},\beta_{m+2},\dots$，得到导数离散的相容阶$\sigma = \min \{ k : k > m,\beta_k \neq 0 \} - m$。

> （习题2.9）设空间步长为$\Delta x$，给定等距分布网格点$x_{j-1},x_j,x_{j+1},x_{j+2}$，利用待定系数方法给出$[u_{xxx}]$在$x_j$和$x_j + \dfrac{\Delta x}{2}$处的差商离散。



函数逼近



符号演算方法

设$h$为给定的移位距离，定位移位算子乘法群$\{ \mathcal{E}^s \}_{s \in \mathbb{R}}$
$$
\mathcal{E}^s p(x) = p(x+sh)
$$
易见$\mathcal{E}^{a+b} = \mathcal{E}^a \cdot \mathcal{E}^b,\forall a,b \in \mathbb{R},\mathcal{E}^0 = 1$

- 一阶向前差分算子$\Delta_+ = \mathcal{E} - 1$
- 一阶向后差分算子$\Delta_- = 1 - \mathcal{E}^{-1}$
- 一步中心差分算子$\Delta_0 = \mathcal{E} - \mathcal{E}^{-1}$
- 半步中心差分算子$\delta = \mathcal{E}^\frac{1}{2} - \mathcal{E}^{-\frac{1}{2}}$
- 二阶中心差分算子$\delta^2 = \mathcal{E} - 2 \cdot 1 + \mathcal{E}^{-1}$

> （习题2.4）热传导方程$u_t=u_{xx}+f$有Hermite格式
> $$
> \left( 1 + \dfrac{1}{12} \delta_x^2 \right) ( u_j^{n+1} - u_j^n ) = \dfrac{1}{2} \mu \delta_x^2 ( u^{n+1}_j + u^n_j ) + \dfrac{1}{2} \Delta t \left[ f_j^{n+1} + \left( 1 + \dfrac{1}{6} \delta_x^2 \right) f_j^n \right]
> $$
> 若网比$\mu$是固定的常数，那么其具有局部截断误差$O \left( \left( \Delta x \right)^4 \right)$​。

证明：
$$
\begin{align}
\tau_j^n =& \dfrac{\left( 1 + \dfrac{1}{12} \delta_x^2 \right) \left( [u]_j^{n+1} - [u]_j^n \right)}{\Delta t}
- \dfrac{\dfrac{1}{2} \delta_x^2 \left( [u]^{n+1}_j + [u]^n_j \right)}{\left( \Delta x \right)^2}
- \dfrac{1}{2} \left[ f_j^{n+1} + \left( 1 + \dfrac{1}{6} \delta_x^2 \right) f_j^n \right] \\
=& \dfrac{\left( [u]_j^{n+1} - [u]_j^n \right) + \dfrac{1}{12} \left[ \left( [u]_{j+1}^{n+1} - [u]_{j+1}^n \right) - 2 \left( [u]_j^{n+1} - [u]_j^n \right) + \left( [u]_{j-1}^{n+1} - [u]_{j-1}^n \right) \right]}{\Delta t} \\
&- \dfrac{\dfrac{1}{2} \left[ \left( [u]^{n+1}_{j+1} + [u]^n_{j+1} \right) - 2 \left( [u]^{n+1}_j + [u]^n_j \right) + \left( [u]^{n+1}_{j-1} + [u]^n_{j-1} \right) \right]}{\left( \Delta x \right)^2} \\
&- \dfrac{1}{2} \left[ f_j^{n+1} + f_j^n + \dfrac{1}{6} \left( f_{j+1}^n - 2 f_j^n + f_{j-1}^n \right) \right] \\
=& \dfrac{ \dfrac{1}{12} \left( [u]_{j+1}^{n+1} - [u]_{j+1}^n \right) + \dfrac{5}{6} \left( [u]_j^{n+1} - [u]_j^n \right) + \dfrac{1}{12} \left( [u]_{j-1}^{n+1} - [u]_{j-1}^n \right) }{\Delta t} \\
&- \dfrac{ \left( \dfrac{1}{2} [u]_{j+1}^{n+1} - [u]_j^{n+1} + \dfrac{1}{2} [u]_{j-1}^{n+1} \right) + \left( \dfrac{1}{2}[u]_{j+1}^n - [u]_j^n + \dfrac{1}{2}[u]_{j-1}^n \right)}{\left( \Delta x \right)^2} \\
&- \dfrac{1}{2} \left[ f_j^{n+1} + f_j^n + \dfrac{1}{6} \left( f_{j+1}^n - 2 f_j^n + f_{j-1}^n \right) \right]
\end{align}
$$

$$
\begin{align}
[u]_{j+k}^{n+1} =& [u]_{j+k}^n + \left[D_t u \right]_{j+k}^n \Delta t + \dfrac{1}{2!} \left[D_t^2 u \right]_{j+k}^n \left( \Delta t \right)^2 + O \left( \left( \Delta t \right)^3 \right) ,k=-1,0,1\\
[u]_{j \pm 1}^{n+1} =& [u]_j^{n+1} \pm \left[D_x u \right]_j^{n+1} \Delta x + \dfrac{1}{2!} \left[D_x^2 u \right]_j^{n+1} \left( \Delta x \right)^2 \pm \dfrac{1}{3!} \left[D_x^3 u \right]_j^{n+1} \left( \Delta x \right)^3 \\
&+ \dfrac{1}{4!} \left[D_x^4 u \right]_j^{n+1} \left( \Delta x \right)^4 \pm \dfrac{1}{5!} \left[D_x^5 u \right]_j^{n+1} \left( \Delta x \right)^5 + O \left( \left( \Delta x \right)^6 \right) \\
[u]_{j \pm 1}^n =& [u]_j^n \pm \left[D_x u \right]_j^n \Delta x + \dfrac{1}{2!} \left[D_x^2 u \right]_j^n \left( \Delta x \right)^2 \pm \dfrac{1}{3!} \left[D_x^3 u \right]_j^n \left( \Delta x \right)^3 \\
&+ \dfrac{1}{4!} \left[D_x^4 u \right]_j^n \left( \Delta x \right)^4 \pm \dfrac{1}{5!} \left[D_x^5 u \right]_j^n \left( \Delta x \right)^5 + O \left( \left( \Delta x \right)^6 \right) \\
f_{j \pm 1}^n =& f_j^n \pm D_x f_j^n \Delta x + \dfrac{1}{2!} D_x^2 f_j^n \left( \Delta x \right)^2 \pm \dfrac{1}{3!} D_x^3 f_j^n \left( \Delta x \right)^3 + O \left( \left( \Delta x \right)^4 \right)
\end{align}
$$

$$
\begin{align}
\tau_j^n =& \dfrac{\dfrac{1}{12} \left( \left[D_t u \right]_{j+1}^n \Delta t + \dfrac{1}{2!} \left[D_t^2 u \right]_{j+1}^n \left( \Delta t \right)^2 + O \left( \left( \Delta t \right)^3 \right) \right) + \dfrac{5}{6}\left( \left[D_t u \right]_j^n \Delta t + \dfrac{1}{2!} \left[D_t^2 u \right]_j^n \left( \Delta t \right)^2 + O \left(\left( \Delta t \right)^3 \right) \right) + \dfrac{1}{12} \left( \left[D_t u \right]_{j-1}^n \Delta t + \dfrac{1}{2!} \left[D_t^2 u \right]_{j-1}^n \left( \Delta t \right)^2 + O \left( \left( \Delta t \right)^3 \right) \right) }{\Delta t} \\
&- \dfrac{\left[ \dfrac{1}{2!} \left[D_x^2 u\right]_j^{n+1} \left( \Delta x \right)^2 + \dfrac{1}{4!} \left[ D_x^4 u \right]_j^{n+1} \left( \Delta x \right)^4 + O \left( \left( \Delta x \right)^6 \right) \right] + \left[ \dfrac{1}{2!} \left[D_x^2 u\right]_j^n \left( \Delta x \right)^2 + \dfrac{1}{4!} \left[ D_x^4 u \right]_j^n \left( \Delta x \right)^4 + O \left( \left( \Delta x \right)^6 \right) \right]}{\left( \Delta x \right)^2} \\
&- \dfrac{1}{2} \left[ f_j^{n+1} + f_j^n + \dfrac{1}{6} \dfrac{2}{2!} D_x^2 f_j^n \left( \Delta x \right)^2 + O \left( \left( \Delta x \right)^4 \right) \right] \\
=& \left[ \dfrac{1}{12} \left( \left[D_t u \right]_{j+1}^n + \dfrac{1}{2} \left[D_t^2 u \right]_{j+1}^n \Delta t \right) + \dfrac{5}{6}\left( \left[D_t u \right]_j^n + \dfrac{1}{2} \left[D_t^2 u \right]_j^n \Delta t \right) + \dfrac{1}{12} \left( \left[D_t u \right]_{j-1}^n + \dfrac{1}{2} \left[D_t^2 u \right]_{j-1}^n \Delta t \right) + O \left( \left( \Delta t \right)^2 \right) \right] \\
&- \left[ \dfrac{1}{2} \left[D_x^2 u\right]_j^{n+1} + \dfrac{1}{24} \left[ D_x^4 u \right]_j^{n+1} \left( \Delta x \right)^2 + \dfrac{1}{2} \left[D_x^2 u\right]_j^n + \dfrac{1}{24} \left[ D_x^4 u \right]_j^n \left( \Delta x \right)^2 + O \left( \left( \Delta x \right)^4 \right) \right] \\
&- \dfrac{1}{2} \left[ f_j^{n+1} + f_j^n + \dfrac{1}{6} D_x^2 f_j^n \left( \Delta x \right)^2 + O \left( \left( \Delta x \right)^4 \right) \right] \\
=& \left( \dfrac{1}{12} \left[D_t u \right]_{j+1}^n + \dfrac{5}{6} \left[D_t u \right]_j^n + \dfrac{1}{12} \left[D_t u \right]_{j-1}^n + \dfrac{1}{24} \left[D_t^2 u \right]_{j+1}^n \Delta t  + \dfrac{5}{12} \left[D_t^2 u \right]_j^n \Delta t + \dfrac{1}{24} \left[D_t^2 u \right]_{j-1}^n \Delta t \right) \\
&- \left( \dfrac{1}{2} \left[ D_x^2 u\right]_j^{n+1} + \dfrac{1}{2} \left[D_x^2 u\right]_j^n + \dfrac{1}{24} \left[ D_x^4 u \right]_j^{n+1} \left( \Delta x \right)^2 + \dfrac{1}{24} \left[ D_x^4 u \right]_j^n \left( \Delta x \right)^2 \right) \\
&- \left( \dfrac{1}{2} f_j^{n+1} + \dfrac{1}{2} f_j^n + \dfrac{1}{12} D_x^2 f_j^n \left( \Delta x \right)^2 \right) + O \left( \left( \Delta t \right)^2 \right) + O \left( \left( \Delta x \right)^4 \right) \\
\end{align}
$$

$$
\begin{align}
\left[ D_t u \right]_{j \pm 1}^n =& \left[ D_t u \right]_j^n \pm \left[ D_x D_t u \right]_j^n \Delta x + \dfrac{1}{2!} \left[ D_x^2 D_t u \right]_j^n \left( \Delta x \right)^2 \pm \dfrac{1}{3!} \left[ D_x^3 D_t u \right]_j^n \left( \Delta x \right)^3 + O \left( \left( \Delta x \right)^4 \right) \\
\left[ D_t^2 u \right]_{j \pm 1}^n =& \left[ D_t^2 u \right]_j^n \pm \left[ D_x D_t^2 u \right]_j^n \Delta x  + O \left( \left( \Delta x \right)^2 \right) \\
\left[ D_x^2 u \right]_j^{n+1} =& \left[ D_x^2 u \right]_j^n + \left[ D_t D_t^2 u \right]_j^n \Delta t  + O \left( \left( \Delta t \right)^2 \right) \\
\left[ D_x^4 u \right]_j^{n+1} =& \left[ D_x^4 u \right]_j^n + O \left( \Delta t \right) \\
f_j^{n+1} =& f_j^n + D_t f_j^n \Delta t + O \left( \left( \Delta t \right)^2 \right)
\end{align}
$$

$$
\begin{align}
\tau_j^n =& \left[ \left[D_t u \right]_j^n + \dfrac{1}{12} \left( \dfrac{2}{2!} \left[ D_x^2 D_t u \right]_j^n \left( \Delta x \right)^2 + O\left( \left( \Delta x \right)^4 \right) \right) + \dfrac{5}{12} \left[D_t^2 u \right]_j^n \Delta t + \dfrac{1}{24} \left( 2 \left[ D_t^2 u\right]_j^n \right) \Delta t \right] \\
&- \left[ \dfrac{1}{2} \left( \left[ D_x^2 u \right]_j^n + \left[ D_t D_x^2 u \right]_j^n \Delta t  + O \left( \left( \Delta t \right)^2 \right) \right) + \dfrac{1}{2} \left[D_x^2 u\right]_j^n + \dfrac{1}{24} \left( \left[ D_x^4 u \right]_j^n + O \left( \Delta t \right) \right) \left( \Delta x \right)^2 + \dfrac{1}{24} \left[ D_x^4 u \right]_j^n \left( \Delta x \right)^2 \right] \\
&- \left[ \dfrac{1}{2} \left( f_j^n + D_t f_j^n \Delta t + O \left( \left( \Delta t \right)^2 \right) \right) + \dfrac{1}{2} f_j^n + \dfrac{1}{12} D_x^2 f_j^n \left( \Delta x \right)^2 \right] + O \left( \left( \Delta t \right)^2 \right) + O \left( \left( \Delta x \right)^4 \right) \\
=& \left( \left[D_t u \right]_j^n + \dfrac{1}{12} \left[ D_x^2 D_t u \right]_j^n \left( \Delta x \right)^2 + \dfrac{1}{2} \left[ D_t^2 u\right]_j^n \Delta t \right) \\
&- \left( \left[ D_x^2 u \right]_j^n + \dfrac{1}{2} \left[ D_t D_x^2 u \right]_j^n \Delta t + \dfrac{1}{12} \left[ D_x^4 u \right]_j^n \left( \Delta x \right)^2 \right) \\
&- \left( f_j^n + \dfrac{1}{2} D_t f_j^n \Delta t + \dfrac{1}{12} D_x^2 f_j^n \left( \Delta x \right)^2 \right) + O \left( \left( \Delta t \right)^2 \right) + O \left( \left( \Delta x \right)^4 \right) \\
=& O \left( \left( \Delta t \right)^2 \right) + O \left( \left( \Delta x \right)^4 \right) \\
=&  O \left( \left( \Delta x \right)^4 \right) \qquad \mathbf{Q.E.D.}
\end{align}
$$

### 2.2 稳定性

差分格式的稳定性与定解问题的真解无关，是自身的固有性质，用来刻画数值解关于定解数据的连续依赖性。根据线性叠加原理，稳定性可拆分为初值稳定性和右端项稳定性。考虑线性差分格式
$$
u^{n+1}=Bu^n+\Delta tH^n
$$

并给定离散范数$\left\|\cdot\right\|=\left\|\cdot\right\|_{\Delta x}$​​。

#### 2.2.1 稳定性分解刻画

>定义：
>
>- 对于对应的齐次线性差分格式
> $$
>  u^{n+1}=Bu^n,\quad n=0:N-1
> $$
>  当$(\Delta x,\Delta t)\rightarrow0$时，若其数值解满足
> $$
>  \left\| u^n \right\| \le K \left\| u^0 \right\|,\quad n=0:N
> $$
>  其中常数$K>0$与$\Delta x,\Delta t,u^0$皆无关，则称对应的差分格式按$\left\|\cdot\right\|$​模具有初值稳定性。
>
>- 设差分格式满足$u^0\equiv0$，当$(\Delta x,\Delta t)\rightarrow0$时，若其数值解满足
> $$
>  \left\| u^n \right\| \le M \sum_{m=0}^{n-1} \left\| H^m \right\| \Delta t,\quad n=1:N
> $$
>  其中常数$M>0$与$\Delta x,\Delta t$皆无关，则称对应的差分格式按$\left\|\cdot\right\|$模具有右端项稳定性。

由于定解条件的舍入误差，若差分格式具有初值稳定性，当$\{u^n\}_{n=0}^N$和$\{w^n\}_{n=0}^N$均满足差分格式时，
$$
\left\|u^n-w^n\right\|\le K\left\|u^0-w^0\right\|,\quad n=0:N
$$
即每一时间层的误差均可被初值扰动控制。

> 齐次线性差分格式的初值稳定性蕴含对应非齐次线性差分格式的右端项稳定性。

证明：Duhamel原理的离散化表述

常见稳定性分析方法：

1. 离散最大模原理
2. Fourier方法
3. 能量方法
4. 修正方程方法
5. 直接矩阵方法
6. 分离变量方法
7. GKS方法

#### 2.2.2 离散最大模原理

设$f \equiv 0$，此时希望数值格式继承连续问题的最大模原理：数值解的离散最大模不增。

> - 离散最大模原理是最大模稳定的充分条件（$K=1$），不一定是必要条件。
> - 任意网格点的差分方程均满足凸组合系数结构可以导出离散最大模原理（参考全隐格式推导）。

> 当$f\equiv0$时，利用离散最大模原理，证明(HP)，(HI)，(HD)全显格式具有最大模初值稳定性当且仅当$\mu a\le\dfrac{1}{2}$。

证明：对于(HP)，任一网格点处均满足
$$
u_j^{n+1}=\mu au_{j-1}^n+(1-2\mu a)u_j^n+\mu au_{j+1}^n
$$
当$\mu a\le\dfrac{1}{2}$时，上式表明$u_j^{n+1}$可表为$u^n_{j-1},u^n_j,u^n_{j+1}$的凸组合，从而
$$
\left|u_j^{n+1}\right|\le\max\left\{\left|u_{j-1}^n\right|,\left|u_j^n\right|,\left|u_{j+1}^n\right|\right\}
$$
故$\left\|u^{n+1}\right\|_{\infty}\le\left\|u^n\right\|_{\infty}$，全显格式具有最大模初值稳定性。

另一方面，只需证明$\mu a>\dfrac{1}{2}$不能保证最大模初值稳定性。取初值$u_j^0=(-1)^j,j=0:J$，
$$
\begin{align}
u_j^1&=\mu a[(-1)^{j+1}+(-1)^{j-1}]+(1-2\mu a)(-1)^j \\
&=(1-4\mu a)(-1)^j
\end{align}
$$
归纳可得$u^n_j=(1-4\mu a)^n(-1)^j$。由于$\mu a>\dfrac{1}{2}$，$\left|1-4\mu a\right|>1$，$\left|u_j^n\right|\rightarrow+\infty\quad(n\rightarrow\infty)$，故是不稳定的。 $\qquad \mathbf{Q.E.D}$

> 当$f\equiv0$时，利用离散最大模原理，证明(HP)，(HI)，(HD)的全隐格式无条件具有最大模初值稳定性。

证明：对于(HP)，任一网格点处均满足
$$
(1+2\mu a)u_j^{n+1}=\mu au_{j-1}^n+u_j^n+\mu au_{j+1}^n
$$
即
$$
u_j^{n+1} = \dfrac{\mu a}{1 + 2 \mu a} u_{j-1}^n + \dfrac{1}{1 + 2 \mu a} u_j^n + \dfrac{\mu a}{1 + 2 \mu a} u_{j+1}^n
$$
（凸组合性质）不妨设$\left| u_{j_0}^{n+1} \right| = \left\| u^{n+1} \right\|_\infty$，那么
$$
\begin{align}
	\left\| u^{n+1} \right\|_\infty = \left| u_{j_0}^{n+1} \right| &= \left| \dfrac{\mu a}{1 + 2 \mu a} u_{j_0 - 1}^n + \dfrac{1}{1 + 2 \mu a} u_{j_0}^n + \dfrac{\mu a}{1 + 2 \mu a} u_{j_0 + 1}^n \right| \\
	& \le \dfrac{\mu a}{1 + 2 \mu a} \left\| u^n \right\|_\infty + \dfrac{1}{1 + 2 \mu a} \left\| u^n \right\|_\infty + \dfrac{\mu a}{1 + 2 \mu a} \left\| u^n \right\|_\infty \\
	& \le \left\| u^n \right\|_\infty
\end{align}
$$
故全显格式具有最大模初值稳定性。$\qquad \mathbf{Q.E.D}$

> 当$f \equiv 0$时，热传导方程$u_t=u_{xx}$的Hermite格式
> $$
> \left( 1 + \dfrac{1}{12} \delta_x^2 \right) ( u_j^{n+1} - u_j^n ) = \dfrac{1}{2} \mu \delta_x^2 ( u^{n+1}_j + u^n_j )
> $$
> 若$\mu \in \left[ \dfrac{1}{6}, \dfrac{5}{6} \right]$，那么格式保持离散最大模原理。

证明：计算得到
$$
\left( \dfrac{5}{6} + \mu \right) u_j^{n+1}  = \left( -\dfrac{1}{12} + \dfrac{1}{2} \mu \right) u_{j+1}^{n+1} + \left( -\dfrac{1}{12} + \dfrac{1}{2} \mu \right) u_{j-1}^{n+1} + \left( \dfrac{1}{12} + \dfrac{1}{2} \mu \right) u_{j+1}^n + \left( \dfrac{5}{6} - \mu \right) u_j^n + \left( \dfrac{1}{12} + \dfrac{1}{2} \mu \right) u_{j-1}^n
$$
当$\mu \in \left[ \dfrac{1}{6}, \dfrac{5}{6} \right]$，右端系数皆正且左端系数不低于右端系数和，从而格式保持离散最大模原理。$\qquad \mathbf{Q.E.D}$

#### 2.2.3 离散$L^2$模的Fourier方法

考虑等距时空网格$\mathcal{T}_{\Delta x,\Delta t}$上的线性常系数双层格式
$$
\sum_{s=-l_1}^{r_1} a_s u_{j+s}^{n+1} = \sum_{s=-l_0}^{r_0} b_s u_{j+s}^n
$$
按时间推进计算，考虑在$\mathbb{R}$上等距网络$\mathcal{T}_{\Delta x}$和其上的网格函数$u=\{u_m\}_{m=-\infty}^{+\infty}$。假设$u$的离散$L^2$模有限，即
$$
\left\| u \right\|_2 = \left\| u \right\|_{2,\Delta x} = \left( \sum_{m=-\infty}^{+\infty} \left| u_m \right|^2 \right)^{\frac{1}{2}} < +\infty
$$
对$u$作逐点常值延拓，得到阶梯函数$\tilde{u}:\mathbb{R} \rightarrow \mathbb{R}$
$$
\tilde{u}(x)=u_m,\quad x \in \left( \left( m - \dfrac{1}{2} \right) \Delta x,\left( m + \dfrac{1}{2} \right) \Delta x \right),\quad \forall m
$$
从而$\tilde{u} \in L^2(\mathbb{R})$，进而存在$L^2(\mathbb{R})$函数$\hat{u}(k):\mathbb{R} \rightarrow \mathbb{C}$​
$$
\hat{u}(k) = \mathcal{F} \tilde{u}(x) = \dfrac{1}{\sqrt{2\pi}} \int_{-\infty}^{+\infty} e^{-ikx} \tilde{u}(x) \mathrm{d}x \\
\tilde{u}(x) = \mathcal{F}^{-1} \hat{u}(k) = \dfrac{1}{\sqrt{2\pi}} \int_{-\infty}^{+\infty} e^{ikx} \hat{u}(k) \mathrm{d}k
$$
据Parseval恒等式
$$
\left\| u \right\|^2_{2,\Delta x} = \left\| \tilde{u} \right\|^2_{L^2(\mathbb{R})} = \left\| \hat{u} \right\|^2_{L^2(\mathbb{R})} = \int_{-\infty}^{+\infty} \left| \hat{u}(k) \right|^2 \mathrm{d}k
$$

##### 2.2.3.1 增长因子

根据离散格式，在$\mathbb{R}$上几乎处处成立
$$
\sum_{s=-l_1}^{r_1} a_s \tilde{u}^{n+1} (x + s \Delta x) = \sum_{s=-l_0}^{r_0} b_s \tilde{u}^n (x + s \Delta x) \\
\sum_{s=-l_1}^{r_1} a_s \mathcal{F}[\tilde{u}^{n+1} (x + s \Delta x)] = \sum_{s=-l_0}^{r_0} b_s \mathcal{F}[\tilde{u}^n (x + s \Delta x)]
$$
由于平移性质$\mathcal{F} [\tilde{u}^n (x + s \Delta x)] = e^{isk \Delta x} \mathcal{F} [\tilde{u}^n (x)] = e^{isk \Delta x} \hat{u}^n(k)$，
$$
\sum_{s=-l_1}^{r_1} a_s e^{isk \Delta x} \hat{u}^{n+1}(k) = \sum_{s=-l_0}^{r_0} b_s e^{isk \Delta x} \hat{u}^n(k)
$$
定义增长因子$\lambda(k)=\dfrac{\hat{u}^{n+1} (k)}{\hat{u}^n (k)}$，其绝对值是波数为$k$的简谐波单步时间推进后的振幅变化率，那么
$$
\lambda(k)=\dfrac{\hat{u}^{n+1} (k)}{\hat{u}^n (k)}=\dfrac{\sum_{s=-l_0}^{r_0} b_s e^{isk \Delta x}}{\sum_{s=-l_1}^{r_1} a_s e^{isk \Delta x}}
$$
实际计算中，可任取波数$k \in \mathbb{R}$，将模态解$u_j^n = \lambda^n e^{ikj \Delta x}$代入差分格式，按上式计算增长因子。

##### 2.2.3.2 von Neumann条件

设$T>0$是给定的终止时刻。

> 定理：下列命题等价。
>
> 1. 双层格式具有$L^2$模稳定性。
>
> 2. 存在常数$K>0$​​使
>    $$
>    \left| \lambda(k) \right|^n \le K,\quad \forall k \in \mathbb{R},\quad \forall n \le \dfrac{T}{\Delta t}
>    $$
>
> 3. 当$\Delta t$充分小时，存在与$k,\Delta t$无关的常数$C \ge 0$​使
>    $$
>    \left| \lambda(k) \right| \le 1 + C \Delta t,\quad \forall k \in \mathbb{R}
>    $$
>
> > 命题2中的$K$可能依赖$T$。
> >
> > 命题3称为von Neumann条件。特别地，当$C=0$时称为严格的von Neumann条件。

证明：双层格式数值解满足
$$
\begin{align}
\left\| u^n \right\|_{2,\Delta x} &= \left\| \hat{u}^n \right\|_{L^2(\mathbb{R})} \le \sup_{k \in \mathbb{R}} \left| \lambda(k) \right| \left\| \hat{u}^{n-1} \right\|_{L^2(\mathbb{R})} \\ 
&\le \cdots \\
&\le \left[ \sup_{k \in \mathbb{R}} \left| \lambda(k) \right| \right]^n \left\| \hat{u}^0 \right\|_{L^2(\mathbb{R})} = \left[ \sup_{k \in \mathbb{R}} \left| \lambda(k) \right| \right]^n \left\| \hat{u}^0 \right\|_{2,\Delta x}
\end{align}
$$




> 当$f\equiv 0$时，利用Fourier方法，证明(HI)全显格式具有$L^2$模初值稳定性当且仅当$\mu a \le \dfrac{1}{2}$。

证明：对于全显格式$\Delta_t u_j^n = \mu a \delta_x^2 u_j^n$，即
$$
u_j^{n+1} = \mu a u_{j+1}^n + ( 1 - 2 \mu a) u_j^n + \mu a u_{j-1}^n \\

\begin{align}
\lambda(k) &= \dfrac{\mu a e^{ik\Delta x} + (1 - 2\mu a) + \mu a e^{-ik\Delta x}}{1} \\
&= (1 - 2\mu a) + 2\mu a \cos (k \Delta x) \\
&= 1 - 4 \mu a \sin^2 \left( \dfrac{k}{2}\Delta x \right)
\end{align}
$$
增长因子中无显式$\Delta t$，对应的von Neumann条件为
$$
-1 \le 1 - 4 \mu a \sin^2 \left( \dfrac{k}{2}\Delta x \right) \le 1
$$
即$\mu a \le \dfrac{1}{2}$。$\qquad \mathbf{Q.E.D}$

> 当$f\equiv 0$时，利用Fourier方法，证明(HI)全隐格式无条件具有$L^2$模初值稳定性。

证明：对于全隐格式$\Delta_t u_j^n = \mu a \delta_x^2 u_j^{n+1}$，即
$$
-\mu a u_{j+1}^{n+1} + (1 + 2 \mu a)u_j^{n+1} - \mu a u_{j-1}^{n+1} = u_j^n \\
\begin{align}
\lambda(k) &= \dfrac{1}{-\mu a e^{ik\Delta x} + (1 + 2\mu a) - \mu a e^{-ik\Delta x}} \\
&= \left[ 1 + 4 \mu a \sin^2 \left( \dfrac{k}{2}\Delta x \right) \right]^{-1}
\end{align}
$$
增长因子中无显式$\Delta t$，对应的von Neumann条件为
$$
\left[ 1 + 4 \mu a \sin^2 \left( \dfrac{k}{2}\Delta x \right) \right]^{-1} \le 1
$$
而这是处处成立的。$\qquad \mathbf{Q.E.D}$

#### 2.2.4 直接矩阵方法

>当$f\equiv0$时，利用直接矩阵方法，证明(HD)全显格式具有$L^2$模初值稳定性当且仅当$\mu a \le \dfrac{1}{2}$。



#### 2.2.5 分离变量方法

> 当$f\equiv0$时，利用分离变量方法，证明(HD)全隐格式无条件具有$L^2$模初值稳定性。



### 2.3 收敛性

> 定义：
>
> 

> （Lax-Richtmyer等价定理）假设线性微分方程定解问题是适定的。若线性差分格式是相容的，则稳定性和收敛性是等价的，且误差阶不低于相容阶。







> 考虑一维反应扩散方程$u_{t}=u_{xx}+\alpha u$的纯初值问题或周期边值问题，其中$\alpha = \pm 1$。
>
> 1. 基于等距时空网络，构造其全显格式和全隐格式。
> 2. 建立$L^2$模稳定性结论。

1. 取等距时空网络$\mathcal{T}_{\Delta x,\Delta t}$。对于全显格式，以$(x_j,t^n)$为离散焦点，取差商离散
   $$
   [u_t]_j^n = \dfrac{\Delta_+ [u]_j^n}{\Delta t} \\
   [u_{xx}]_j^n = \dfrac{\delta^2 [u]_j^n}{\left( \Delta x \right)^2} \\
   $$
   据$[u_t]_j^n = [u_{xx}]_j^n +\alpha [u]_j^n$成立
   $$
   \dfrac{\Delta_+ [u]_j^n}{\Delta t} = \dfrac{\delta^2 [u]_j^n}{\left( \Delta x \right)^2} + \alpha [u]_j^n
   $$
   整理得到全显格式
   $$
   u_j^{n+1} = \mu u_{j-1}^n + (1 - 2\mu + \alpha \Delta t)u_j^n + \mu u_{j+1}^n
   $$
   对于全隐格式，以$(x_j,t^{n+1})$​为离散焦点，取差商离散
   $$
   [u_t]_j^{n+1} = \dfrac{\Delta_- [u]_j^{n+1}}{\Delta t} \\
   [u_{xx}]_j^{n+1} = \dfrac{\delta^2 [u]_j^{n+1}}{\left( \Delta x \right)^2} \\
   $$
   据$[u_t]_j^{n+1} = [u_{xx}]_j^{n+1} +\alpha [u]_j^{n+1}$​成立
   $$
   \dfrac{\Delta_- [u]_j^{n+1}}{\Delta t} = \dfrac{\delta^2 [u]_j^{n+1}}{\left( \Delta x \right)^2} + \alpha [u]_j^{n+1}
   $$
   整理得到全隐格式
   $$
   -\mu u_{j-1}^{n+1} + (1 + 2\mu - \alpha \Delta t) u_j^{n+1} - \mu u_{j+1}^{n+1} = u_j^n
   $$

2. 对于全显格式，其增长因子
   $$
   \begin{align}
   \lambda(k) &= \dfrac{\mu e^{-ik\Delta x} + (1 - 2\mu + \alpha \Delta t) + \mu e^{ik\Delta x}}{1} \\
   &= (1 - 2\mu + \alpha \Delta t) + 2 \mu \cos (k \Delta x) \\
   &= (1 + \alpha \Delta t) - 4 \mu \sin^2 \left( \dfrac{k}{2} \Delta x \right)
   \end{align}
   $$

   $$
   \begin{align}
   \left| \lambda(k) \right| \le \left| 1 - 4 \mu \sin^2 \left( \dfrac{k}{2} \Delta x \right) \right| + \left| \alpha \right| \Delta t
   \end{align}
   $$

   据von Neumann条件，当$\mu \le \dfrac{1}{2}$时，无论$\alpha$取何值，格式均具有$L^2$模稳定性。

   对于全隐格式，其增长因子
   $$
   \begin{align}
   \lambda(k) &= \dfrac{1}{-\mu e^{-ik\Delta x} + (1 + 2\mu - \alpha \Delta t) - \mu e^{ik\Delta x}} \\
   &= \left[ (1 + 2\mu - \alpha \Delta t) - 2 \mu \cos (k \Delta x) \right]^{-1} \\
   &= \left[ (1 - \alpha \Delta t) + 4 \mu \sin^2 \left( \dfrac{k}{2} \Delta x \right) \right]^{-1}
   \end{align}
   $$
   当$\Delta t$充分小时，$\left| \lambda(k) \right| \le 1$无条件成立。据严格von Neumann条件，格式均具有$L^2$模稳定性。



## 3 热传导方程

考虑线性常系数热传导方程
$$
u_t = au_{xx},\quad a>0
$$
取等距时空网格$\mathcal{T}_{\Delta x,\Delta t}$，记网比$\mu=\dfrac{\Delta t}{\left( \Delta x \right)^2}$。

### 3.1 相容性

3.1.1 加权平均格式

设$\theta \in [0,1]$为给定权重，将全显格式和全隐格式加权平均即得加权平均格式
$$
\Delta_t u^n_j = \theta \mu a \delta_x^2 u_j^{n+1} + (1 - \theta) \mu a \delta_x^2 u^n_j 
$$

> 加权平均格式至少具有(2,1)阶局部截断误差。

证明：
$$
\begin{align}
\tau_j^n &= \dfrac{[u]_j^{n+1}-[u]_j^n}{\Delta t} - \dfrac{\theta a \delta_x^2 u_j^{n+1} + (1 - \theta) a \delta_x^2 u_j^n}{\left( \Delta x \right)^2} \\
&= \dfrac{1}{\Delta t} \left( [u]_j^{n+1} - [u]_j^n \right) - \dfrac{a}{\left( \Delta x \right)^2}\left[ \dfrac{1}{2} \delta_x^2 \left( [u]_j^{n+1} + [u]_j^n \right) + \left( \theta - \dfrac{1}{2} \right) \delta_x^2 \left( [u]_j^{n+1} - [u]_j^n \right) \right]
\end{align}
$$
以$(x_j,t^{n+\frac{1}{2}})$为展开中心，计算得到
$$
\begin{align}
\tau_j^n &= -\Delta t \left( \theta - \dfrac{1}{2} \right) [D_t^2 u]_j^{n + \frac{1}{2}} - \dfrac{a}{12}\left( \Delta x \right)^2 [D_x^4 u]_j^{n + \frac{1}{2}} + O \left( \left( \Delta x \right)^4 + \left( \Delta t \right)^2 \right) \\
&= - \Delta t \left[ \dfrac{1}{12 \mu a} + \theta - \dfrac{1}{2} \right][D_t^2 u]_j^{n+\frac{1}{2}} + O \left( \left( \Delta x \right)^4 + \left( \Delta t \right)^2 \right)
\end{align}
$$
这表明格式至少具有(2,1)阶局部截断误差。$\qquad \mathbf{Q.E.D.}$

取$\theta=\dfrac{1}{2}$​，即得Crank-Nicolson格式（CN格式），其无条件具有(2,2)阶局部截断误差。
$$
\Delta_t u_j^n = \dfrac{1}{2} \mu a \left( \delta_x^2 u_j^{n+1} + \delta_x^2 u_j^n \right)
$$


取$\theta = \dfrac{1}{2} - \dfrac{1}{12 \mu a}$，即得Douglas格式，其无条件具有(4,2)阶局部截断误差，达到整体四阶相容。

> （习题3.1）当$\mu a = \dfrac{1}{\sqrt{20}}$时，Douglas格式具有整体六阶的局部截断误差。

证明：



> 讨论加权平均格式的$L_2$模稳定性。

证明：设$k \in \mathbb{R}$是任意波数，将模态解$u_j^n = \lambda^n e^{i k j \Delta x}$代入
$$
\lambda = \lambda(k) = \dfrac{1 - 4 \mu a (1 - \theta)\sin^2(\dfrac{1}{2} k \Delta x)}{1 + 4 \mu a \theta \sin^2(\dfrac{1}{2} k \Delta x)}
$$

由于加权平均格式为双层线性格式，其具有$L^2$模稳定性当且仅当von Neumann条件成立，其成立当且仅当$\mu a (1 - 2 \theta) \le \dfrac{1}{2}$。

> 当$\theta <\dfrac{1}{2}$时，格式称作偏显格式，有条件稳定
>
> 当$\theta \ge \dfrac{1}{2}$时，格式称作偏隐格式，无条件稳定


> 讨论加权平均格式的最大模稳定性。

证明：将格式改写为
$$
(1 + 2 \theta \mu a) u_j^{n+1} = \theta \mu a (u_{j-1}^{n+1} + u_{j+1}^{n+1}) + \left[1 - 2(1 -\theta) \mu a \right] u^n_j + (1 - \theta) \mu a (u_{j-1}^n + u_{j+1}^n)
$$
设$\left| u_{j_0}^{n+1} \right| = \left\| u^{n+1} \right\|_\infty$，那么当$(1 - \theta) \mu a \le \dfrac{1}{2}$（这表明下述结论为最大模稳定的充分条件）时
$$
\begin{align}
	(1 + 2 \theta \mu a)\left\| u^{n+1} \right\|_\infty = (1 + 2 \theta \mu a) \left| u_{j_0}^{n+1} \right| &\le 2 \theta \mu a \left\| u^{n+1} \right\|_\infty + \left[ \left| 1 - 2(1 - \theta) \mu a \right| + 2 \left| (1 - \theta) \mu a \right| \right] \left\| u^n \right\|_\infty \\
	&\le 2 \theta \mu a \left\| u^{n+1} \right\|_\infty + \left\| u^n \right\|_\infty
\end{align}
$$




3.1.2 三层格式





3.2 计算效率





3.3 误差估计及收敛分析

3.3.1 基于强正则性假设



3.3.2 基于弱正则性假设





> 考虑周期边值问题
> $$
> \begin{cases}
> 	u_t = u_{xx} \\
> 	u(x, t) = u(x + 2 \pi, t),\quad (x,t) \in \mathbb{R} \times [0,1] \\
> 	u(x, 0) = u_0(x),\quad x \in \mathbb{R}
> \end{cases}
> $$
> 其中$u_0(x) = u_0^{(i)}(x), i=1, 2$，$u_0^{(1)}(x)$是间断函数
> $$
> u_0^{(1)}(x) =
> \begin{cases}
> 	1,\quad x \in \left[ -\dfrac{\pi}{2}, \dfrac{\pi}{2} \right] \\
> 	0,\quad x \in \left[ -\pi, -\dfrac{\pi}{2} \right) \cup \left( \dfrac{\pi}{2}, \pi \right]
> \end{cases}
> $$
> $u_0^{(2)}(x)$是导数间断的连续函数
> $$
> u_0^{(2)}(x) = \pi - \left| x \right|,\quad x \in \left[ -\pi,\pi \right]
> $$
> 对上述问题，固定网比$\mu = 0.4$，分别利用全显格式，全隐格式和CN格式，进行数值模拟。

$$
u^{(1)}(x,t) = \dfrac{1}{2} + \sum_{n=1}^{\infty}\dfrac{2 \sin \left( \dfrac{n \pi}{2} \right) }{n \pi} e^{-n^2t} \cos(nx)
$$

$$
u^{(2)}(x,t) = \dfrac{\pi}{2} + \sum_{n=1}^{\infty} \dfrac{4 \sin^2 \left( \dfrac{n \pi}{2} \right)}{n^2 \pi} e^{-n^2 t} \cos(nx)
$$

对于每层的网格函数$\{u^n_j\}_{j=-J}^J$，
$$
u_{-J} = u_J \\
B_{n+1} u^{n+1}_{-J+1:J} = B_n u^n_{-J+1:J}
$$
对于全显格式



对于全隐格式



对于CN格式







```matlab
clear
clc

J = 18;
error = [0, 0, 0, 0, 0];
order = [0, 0, 0, 0];

for round = 1: 5
    space_net = linspace(-pi, pi, 2 * J + 1);
    delta_x = space_net(2) - space_net(1);
    time_net = linspace(0, 1, 1 / (0.4 * pi * pi / J / J));
    J = J * 2;

    u0 = arrayfun(@u0_1, space_net);
    %u0 = arrayfun(@u0_2, space_net);
    u_star = zeros(size(time_net, 2), size(space_net, 2));
    for i = 1: size(time_net, 2)
        for j = 1: size(space_net, 2)
            u_star(i, j) = u_1(space_net(j), time_net(i));
            %u_star(i,j) = u_2(space_net(j), time_net(i));
        end
    end
    u = heat_classical_explicit(1, space_net, time_net, u0);
    %u = heat_classical_implicit(1, space_net ,time_net, u0);
    %u = heat_cn(1, space_net, time_net, u0);

    error(round) = sqrt(sum((u(end, :) - u_star(end, :)) .^ 2) * delta_x);
end
for round = 1: 4
    order(round) = (log(error(round)) - log(error(round + 1))) / log(2);
end
```

```matlab
function u = heat_classical_explicit(a, space_net, time_net, u0)
    delta_x = space_net(2) - space_net(1);
    delta_t = time_net(2) - time_net(1);
    mu = delta_t / delta_x / delta_x;
    u = zeros(size(time_net, 2), size(space_net, 2));
    u(1, :) = u0;
    for i = 2: size(time_net, 2)
        u(i, 1) = mu * a * (u(i - 1, end - 1) + u(i - 1, 2)) + (1 - 2 * mu * a) * u(i - 1, 1);
        for j = 2: size(space_net, 2) - 1
            u(i, j) = mu * a * (u(i - 1, j - 1) + u(i - 1, j + 1)) + (1 - 2 * mu * a) * u(i - 1, j);
        end
        u(i, end) = mu * a * (u(i - 1, end - 1) + u(i - 1, 2)) + (1 - 2 * mu * a) * u(i - 1, end);
    end
end
```

```matlab
function u = heat_classical_implicit(a, space_net, time_net, u0)
    delta_x = space_net(2) - space_net(1);
    delta_t = time_net(2) - time_net(1);
    mu = delta_t / delta_x / delta_x;
    u = zeros(size(time_net, 2), size(space_net, 2));
    u(1, :) = u0;
    B = zeros(size(space_net, 2) - 1);
    for i = 1: size(B, 1)
        if(i == 1)
            B(i, end) = -mu * a;
        else
            B(i, i - 1) = -mu * a;
        end
        B(i, i) = 1 + 2 * mu * a;
        if(i == size(B, 1))
            B(i, 1) = -mu * a;
        else
            B(i, i + 1) = -mu * a;
        end
    end
    for i = 2: size(time_net, 2)
        u(i, 2: end) = linear_solve(B, transpose(u(i - 1, 2: end)));
        u(i, 1) = u(i, end);
    end
end
```

```matlab
function u = heat_cn(a, space_net, time_net, u0)
    delta_x = space_net(2) - space_net(1);
    delta_t = time_net(2) - time_net(1);
    mu = delta_t / delta_x / delta_x;
    u = zeros(size(time_net, 2), size(space_net, 2));
    u(1, :) = u0;
    B = zeros(size(space_net, 2) - 1);
    for i = 1: size(B, 1)
        if(i == 1)
            B(i, end) = -mu * a / 2;
        else
            B(i, i - 1) = -mu * a / 2;
        end
        B(i, i) = 1 + mu * a;
        if(i == size(B, 1))
            B(i, 1) = -mu * a / 2;
        else
            B(i, i + 1) = -mu * a / 2;
        end
    end
    for i = 2: size(time_net, 2)
        temp = zeros(1, size(u, 2));
        temp(1) = mu * a * (u(i - 1, end - 1) + u(i - 1, 2)) / 2 + (1 - mu * a) * u(i - 1, 1);
        for j = 2: size(u, 2) - 1
            temp(j) = mu * a * (u(i - 1, j - 1) + u(i - 1, j + 1)) / 2 + (1 - mu * a) * u(i - 1, j);
        end
        temp(end) = mu * a * (u(i - 1,end - 1) + u(i - 1, 2)) / 2 + (1 - mu * a) * u(i - 1, end);
        u(i, 2: end) = linear_solve(B, transpose(temp(2: end)));
        u(i, 1) = u(i, end);
    end
end
```

```matlab
function x = linear_solve(A, b)
    x = A \ b;
end
```

```matlab
function u=u_1(x,t)
    u=0.5;
    n=1;
    while true
        part=2*sin(n*pi/2)/n/pi*exp(-n^2*t)*cos(n*x);
        if abs(part)<1e-8
            break;
        end
        u = u + part;
        n = n + 1;
    end
end
```

```matlab
function u = u_2(x,t)
    u = pi / 2;
    n = 1;
    while true
        part = 4 * (sin(n * pi / 2) ^ 2) / (n ^ 2)/ pi * exp(-n ^ 2 * t)*cos(n * x);
        if abs(part) < 1e-8
            break;
        end
        u = u + part;
        n = n + 1;
    end
end
```

```matlab
function u = u0_1(x)
    if abs(x) <= pi / 2
        u = 1;
    else
        u = 0;
    end
end
```

```matlab
function u = u0_2(x)
    u = pi - abs(x);
end
```

全显格式

| $J$  |  $u_0=u^{(1)}_0$  |                   |    $u_0=u^{(2)}_0$    |                   |
| :--: | :---------------: | :---------------: | :-------------------: | :---------------: |
|      |       误差        |      误差阶       |         误差          |      误差阶       |
|  18  | 0.070691184165142 |                   | 9.516450685149100e-04 |                   |
|  36  | 0.035068156760164 | 1.011368706867869 | 2.209013314810101e-04 | 2.107021482750239 |
|  72  | 0.017471651701085 | 1.005145600761598 | 6.119404981949514e-05 | 1.851938829613794 |
| 144  | 0.008721112611955 | 1.002431895079183 | 3.355653316092938e-05 | 0.866797706059126 |
| 288  | 0.004357139541173 | 1.001130883659627 | 3.113077189045240e-05 | 0.108252324755212 |

全隐格式

| $J$  |  $u_0=u^{(1)}_0$  |                   |    $u_0=u^{(2)}_0$    |                   |
| :--: | :---------------: | :---------------: | :-------------------: | :---------------: |
|      |       误差        |      误差阶       |         误差          |      误差阶       |
|  18  | 0.070497869834098 |                   |   0.009581607562645   |                   |
|  36  | 0.035044814005542 | 1.008378698220213 |   0.002355477958212   | 2.024247893896323 |
|  72  | 0.017468782833143 | 1.004421877879651 | 5.851302131711442e-04 | 2.009190213119203 |
| 144  | 0.008720758215217 | 1.002253610699342 | 1.489742043041568e-04 | 1.973695171771822 |
| 288  | 0.004357096478469 | 1.001086514773134 | 4.786330905838536e-05 | 1.638070494781345 |

CN格式

| $J$  |  $u_0=u^{(1)}_0$  |                   |    $u_0=u^{(2)}_0$    |                   |
| :--: | :---------------: | :---------------: | :-------------------: | :---------------: |
|      |       误差        |      误差阶       |         误差          |      误差阶       |
|  18  | 0.070544987124860 |                   |   0.004331275707373   |                   |
|  36  | 0.035050564558861 | 1.009105887723245 |   0.001069609350932   | 2.017708026898790 |
|  72  | 0.017469490567238 | 1.004600143858659 | 2.675597101065143e-04 | 1.999151190817265 |
| 144  | 0.008720845195000 | 1.002297669985149 | 7.326001953173807e-05 | 1.868762896464711 |
| 288  | 0.004357106769506 | 1.001097496465550 | 3.518821729232728e-05 | 1.057933659648913 |
