# Bitcoin Alpha 信任网络的拓扑数据分析项目说明

> 本文件用于课程项目的选题、实现和汇报准备，不是可以直接提交的最终报告。项目目标是清楚地展示课程中的距离、层次聚类、Vietoris--Rips 复形和同调方法如何用于一个真实数据集。整体范围以能够完成、解释并通过课程项目为准。

## 1 项目概述

本项目研究 Bitcoin Alpha 平台中用户信任关系随时间的变化。

Bitcoin Alpha 是一个使用 Bitcoin 进行交易的平台。由于用户可以保持匿名，平台允许用户对其他用户给出信任或不信任评分。把用户看作节点，把评分看作带方向、带正负号并带时间戳的边，就得到一个随时间变化的信任网络。

本项目的核心做法是：

1. 把完整数据按照时间划分为若干窗口，例如每半年一个窗口；
2. 在每个窗口中，根据用户信任对象的相似程度定义 Jaccard 距离；
3. 使用该距离进行 single-linkage 层次聚类；
4. 建立 Vietoris--Rips filtration，并计算零维和一维同调；
5. 用 barcode 和 Betti curve 描述每个时间窗口的拓扑结构；
6. 比较不同时间窗口，从而观察信任网络是否存在不同的发展阶段。

项目可以用一句话概括为：

> 使用拓扑数据分析研究 Bitcoin Alpha 用户信任行为的多尺度结构如何随时间变化。

一个合适的英文研究问题是：

> How does the topological structure of the Bitcoin Alpha trust network change over time?

## 2 数据集介绍

### 2.1 数据来源

本项目使用 Stanford Network Analysis Project 提供的 Bitcoin Alpha 数据集：

- 数据页面：[Bitcoin Alpha trust weighted signed network](https://snap.stanford.edu/data/soc-sign-bitcoin-alpha.html)
- 数据文件：`soc-sign-bitcoinalpha.csv.gz`
- 推荐引用：S. Kumar, F. Spezzano, V. S. Subrahmanian, and C. Faloutsos, *Edge Weight Prediction in Weighted Signed Networks*, IEEE ICDM, 2016.

### 2.2 数据规模

根据 SNAP 的数据说明：

| 项目 | 数值 |
|---|---:|
| 用户节点数 | 3,783 |
| 评分边数 | 24,186 |
| 评分范围 | -10 到 +10 |
| 正评分比例 | 约 93% |

这是一个具有以下四种性质的网络：

- **有向网络**：用户 $u$ 对用户 $v$ 的评分不等于 $v$ 对 $u$ 的评分；
- **带权网络**：每条边的评分在 -10 到 +10 之间；
- **带符号网络**：正数表示信任，负数表示不信任；
- **时间网络**：每条评分带有时间戳，可以研究网络随时间的变化。

### 2.3 文件字段

CSV 文件的每一行表示一条评分记录：

```text
SOURCE, TARGET, RATING, TIME
```

各字段的含义如下：

| 字段 | 含义 |
|---|---|
| `SOURCE` | 给出评分的用户编号 |
| `TARGET` | 被评分的用户编号 |
| `RATING` | 信任评分，范围为 -10 到 +10 |
| `TIME` | Unix 时间戳，即从 1970 年 1 月 1 日开始计算的秒数 |

在主分析中，把 `RATING > 0` 定义为正向信任关系。负评分暂时只用于描述每个时间窗口中负评分的数量和比例，不单独构造 signed distance。这样可以让项目保持在课程介绍的 pseudometric 框架中。

## 3 与课程内容的对应关系

课程把输入写成

\[
I=\{X_d\mid d\in D\},
\]

其中 $D$ 是希望分析的对象集合，每个 $X_d$ 是代表对象 $d$ 的有限数据集。

在本项目中：

- $D$ 是时间窗口的集合；
- $d\in D$ 是某一个半年或季度；
- $X_d$ 是该窗口中满足活跃条件的用户集合；
- 每个 $X_d$ 上定义一个 Jaccard pseudometric；
- 最终希望根据拓扑特征比较并聚类这些时间窗口。

项目步骤与课程概念的对应关系如下：

| 项目步骤 | 课程概念 |
|---|---|
| 把交易历史划分为时间窗口 | 输入 $I=\{X_d\mid d\in D\}$ |
| 比较用户信任对象的相似性 | metric 的 restriction 与 pseudometric |
| Jaccard 距离 | 有限集合上的 Jaccard metric |
| 多尺度用户分组 | $X/d_t$、single linkage 与 dendrogram |
| 建立拓扑空间 | Vietoris--Rips complex 与 filtration |
| 计算连通分量和环 | $H_0$ 与 $H_1$ |
| 表示持续结构 | barcode、persistent homology 与 Betti curve |
| 比较不同窗口 | piecewise-constant functions 上的距离 |
| 对窗口进行聚类 | 在 $D$ 上得到 partition 或 dendrogram |

### 3.1 课程框架的直观解释

如果此前没有学习过这门课，可以把整个课程框架理解为以下转换过程：

\[
\text{原始数据}
\longrightarrow
\text{距离}
\longrightarrow
\text{多尺度几何结构}
\longrightarrow
\text{同调特征}
\longrightarrow
\text{比较与聚类}.
\]

每一步解决一个不同的问题：

1. **原始数据是什么？**  
   在本项目里，原始数据是一系列带时间戳的用户评分。我们把它整理成多个时间窗口，每个窗口对应一个有限用户集合。

2. **怎样比较两个用户？**  
   原始用户编号本身没有几何意义，因此需要定义距离。Jaccard distance 根据两位用户信任对象的重合程度，给出一个 0 到 1 之间的数。

3. **为什么要观察多个尺度？**  
   如果强行选择一个距离阈值，只能得到一个 partition，而且阈值通常很任意。TDA 让阈值 $t$ 从小到大变化，从而同时观察非常相似的小群体和较宽松的大群体。

4. **为什么要建立 Vietoris--Rips complex？**  
   距离矩阵能说明哪些用户相近，但不容易直接概括整体形状。Vietoris--Rips construction 把距离关系转换成 simplicial complex，使我们能够使用 homology。

5. **Homology 给出什么？**  
   $H_0$ 记录 connected components，也就是用户群体怎样合并；$H_1$ 记录一维环状结构。Persistent homology 进一步记录这些结构在多大的尺度范围内存在。

6. **怎样比较两个时间窗口？**  
   每个窗口最终得到 Betti curves。把两条曲线之间的面积差作为距离，就可以比较窗口，并再次使用 hierarchical clustering 得到时间阶段。

### 3.2 课程讲义中的具体位置

| 课程讲义位置 | 主要概念 | 本项目中的用途 |
|---|---|---|
| Points 2--3 | 输入 $I=\{X_d\mid d\in D\}$ 和分析 pipeline | 把每个时间窗口看成一个待比较对象 |
| Points 14--15 | distance、pseudometric 与 restriction | 说明用户空间上的距离怎样合法构造 |
| Points 50、52、54 | Jaccard metric 与 Steinhaus transform | 为用户信任集合定义 Jaccard distance |
| Point 25 | partitions $X/d_t$ | 在固定尺度下形成用户群体 |
| Points 114、117、122 | dendrogram、hierarchical clustering、single linkage | 描述用户群体随尺度合并的过程 |
| Points 187--193 | Vietoris--Rips complex 与 filtration | 把距离矩阵变成多尺度 simplicial complexes |
| Points 198--208 | homology、$H_0$ 与 connected components | 计算和解释 $H_0$、$H_1$ |
| Points 135、158、160、164 | tame vector spaces、bar decomposition 与 barcode | 表示拓扑结构的出生和死亡 |
| Points 55、57、60、61 | piecewise-constant functions 与曲线距离 | 比较不同时间窗口的 Betti curves |

Presentation 中不需要逐条证明这些结论。比较合适的做法是：说明项目用了哪些课程对象、这些对象在数据中代表什么，以及最终图形如何解释。

## 4 数据表示与 Jaccard 距离

### 4.1 时间窗口

先把 Unix 时间戳转换为日期，然后把数据划分为连续时间窗口。

建议先使用**半年窗口**。半年窗口数量不会太多，而且通常比月度窗口包含更多评分，适合作为课程项目的默认设置。如果某些半年窗口仍然非常稀疏，可以把相邻窗口合并。

为了保证窗口之间可以比较，应当在项目开始时固定窗口宽度，不根据最后图形是否“好看”反复改变。

### 4.2 用户的正向信任集合

对时间窗口 $d$ 中的用户 $u$，定义

\[
T_d(u)=\{v\mid u\text{ 在窗口 }d\text{ 中对 }v\text{ 给出正评分}\}.
\]

也就是说，$T_d(u)$ 是用户 $u$ 在该窗口内主动信任的用户集合。

如果同一用户在同一窗口内多次评价同一个对象，建议只保留该窗口中时间最晚的一次评分。这样每一对 `SOURCE-TARGET` 在一个窗口内只对应一个最终状态。

### 4.3 活跃用户筛选

如果 $T_d(u)=\varnothing$，Jaccard 距离会把多个空集合视为完全相同，产生一个没有实际意义的零距离大类。因此建议只保留在该窗口内至少给出 3 个正评分的用户：

\[
X_d=\{u\mid |T_d(u)|\geq 3\}.
\]

阈值 3 是一个简单、易解释的默认选择。实现后可以检查每个窗口剩余多少用户；如果人数太少，可以改为 2，但需要在报告中说明。

为了控制 $H_1$ 的计算量，可以进一步保留每个窗口最活跃的前 $N$ 个用户，例如 $N=200$。所有窗口尽量使用相同的 $N$，避免结果主要由用户数量差异造成。

### 4.4 Jaccard 距离

对 $u,v\in X_d$，定义

\[
d_d(u,v)
=1-\frac{|T_d(u)\cap T_d(v)|}{|T_d(u)\cup T_d(v)|}.
\]

解释如下：

- 如果两位用户信任的对象完全相同，则距离为 0；
- 如果两人的信任对象有部分重合，则距离介于 0 和 1 之间；
- 如果两人的信任对象完全不重合，则距离为 1。

Jaccard distance 是有限集合空间上的 metric。把每个用户 $u$ 映射到集合 $T_d(u)$，再把 Jaccard metric 限制到用户集合上，就得到 $X_d$ 上的 pseudometric。它可能不是严格的 metric，因为两位不同用户可能拥有完全相同的信任集合，从而距离为 0。

需要注意：原始网络是有向的，但得到的距离是对称的。方向信息仍然通过“用户主动信任了谁”进入 $T_d(u)$，之后比较的是两位用户的信任行为是否相似。

## 5 层次聚类

对每个窗口的距离矩阵进行 single-linkage hierarchical clustering，可以得到一个 dendrogram。

在尺度 $t$ 下，如果两位用户能够通过一系列距离不超过 $t$ 的用户连接起来，就把他们放在同一 cluster 中。随着 $t$ 增大：

- 小尺度表示信任行为非常相似的用户；
- 中等尺度形成较大的行为群体；
- 大尺度最终把多数用户合并在一起。

这正对应课程中的 partitions $X_d/(d_d)_t$。对同一个距离考察所有尺度，而不是任意选择一个 cut height，是使用 TDA 的主要理由之一。

项目中不必同时完整比较 single、complete 和 average linkage。为了控制工作量，可以只把 single linkage 作为主方法。如果时间充足，再增加 average linkage 作为附加图。

## 6 Vietoris--Rips filtration

对每个 pseudometric space $(X_d,d_d)$，在尺度 $t$ 建立 Vietoris--Rips complex：

\[
VR_t(X_d)=\{\sigma\subseteq X_d\mid d_d(u,v)\leq t
\text{ for all }u,v\in\sigma\}.
\]

直观地说：

- 每个用户是一个顶点；
- 当两个用户的 Jaccard 距离不超过 $t$ 时，在它们之间加入一条边；
- 如果一组用户两两之间都有边，就加入相应的高维 simplex。

当 $t$ 从 0 增大到 1 时，得到一系列嵌套的 simplicial complexes：

\[
VR_0\subseteq VR_{t_1}\subseteq VR_{t_2}\subseteq\cdots\subseteq VR_1.
\]

这就是 Vietoris--Rips filtration。

## 7 同调和拓扑特征

### 7.1 零维同调 $H_0$

$H_0$ 描述 connected components。开始时用户通常彼此分离；随着尺度增加，用户逐渐合并成更大的 connected components。

$H_0$ barcode 中：

- 每条 bar 对应一个 connected component；
- bar 的死亡时间表示该 component 在什么尺度与其他 component 合并；
- 长 bar 表示一个在较大尺度上仍保持独立的用户群体。

$H_0$ 与 single-linkage dendrogram 表达的是同一类合并信息，因此它最容易解释，也应当是项目的主要结果。

### 7.2 一维同调 $H_1$

$H_1$ 描述 filtration 中出现的环状结构。这里的环不是原始信任网络中的直接交易环，而是**用户信任行为相似性空间中的环**。

如果一个 $H_1$ class 持续的尺度很短，它可能只是局部噪声。如果持续较长，它表示在一段尺度范围内存在一个不能立即被三角形填满的循环相似性结构。

课程项目不需要对每一个 $H_1$ bar 给出社会学解释。展示不同时间窗口中 $H_1$ 的数量和持续时间是否不同即可。如果所有 $H_1$ bars 都很短，这也是可以报告的结果。

### 7.3 Betti curves

定义

\[
\beta_k^d(t)=\dim H_k(VR_t(X_d)),\qquad k=0,1.
\]

其中：

- $\beta_0^d(t)$ 是尺度 $t$ 下 connected components 的数量；
- $\beta_1^d(t)$ 是尺度 $t$ 下独立一维环的数量。

由于不同窗口的用户数量可能不同，建议使用归一化的零维曲线：

\[
\widetilde{\beta}_0^d(t)=\frac{\beta_0^d(t)}{|X_d|}.
\]

如果所有窗口固定使用相同数量 $N$ 的用户，原始 $\beta_0$ 也可以直接比较。

## 8 时间窗口之间的比较

可以用 $l_1$ 距离比较两个窗口 $d_i,d_j$ 的 Betti curves：

\[
D_k(d_i,d_j)
=\int_0^1
\left|
\widetilde{\beta}_k^{d_i}(t)-
\widetilde{\beta}_k^{d_j}(t)
\right|dt.
\]

对于 $H_1$，可以使用同样的积分；如果不同窗口的用户数固定，也可以直接比较原始 $\beta_1$。

将所有时间窗口两两比较后，可以得到一个时间窗口距离矩阵。再对这个矩阵进行 hierarchical clustering，就得到时间窗口的 dendrogram。

最后可以观察：

- 相邻时间窗口是否通常聚在一起；
- 是否存在早期和后期两个主要阶段；
- 某个时间窗口是否表现为异常窗口；
- 拓扑变化是否与用户数、评分数或正负评分比例的变化同时出现。

## 9 实际实施流程

### 9.1 实验设计总览

这个实验包含两层分析：

- **窗口内部分析**：把用户作为点，研究一个时间窗口内部的用户群体结构；
- **窗口之间分析**：把每个时间窗口的 Betti curves 作为特征，研究不同时间窗口是否形成不同的发展阶段。

完整的数据流是：

\[
\begin{aligned}
\text{评分记录}
&\longrightarrow \text{时间窗口}
\longrightarrow \text{用户的正向信任集合}\\
&\longrightarrow \text{Jaccard 距离矩阵}
\longrightarrow \text{dendrogram 与 VR filtration}\\
&\longrightarrow H_0,H_1
\longrightarrow \text{Betti curves}\\
&\longrightarrow \text{窗口距离矩阵}
\longrightarrow \text{时间窗口 dendrogram}.
\end{aligned}
\]

第一次实验可以先固定以下设置：

| 实验选择 | 默认设置 | 选择理由 |
|---|---|---|
| 时间窗口 | 半年 | 比月度窗口更不容易过于稀疏 |
| 正向信任 | 评分大于 0 | 定义简单，容易解释 |
| 重复评分 | 保留窗口内最后一次 | 用最新评分代表该窗口的最终关系 |
| 活跃用户条件 | 至少给出 3 个正评分 | 避免大量空信任集合 |
| 每个窗口最大用户数 | 200 | 控制 $H_1$ 计算量 |
| 用户距离 | Jaccard distance | 是课程中的合法 metric construction |
| Homology 维数 | $H_0$ 和 $H_1$ | 足以覆盖课程核心，且容易解释 |
| 系数域 | $\mathbb F_2$ | 课程中常用，计算简单 |
| Filtration 范围 | $t\in[0,1]$ | Jaccard distance 位于该范围 |
| 用户聚类方法 | single linkage | 与 $H_0$ 和 $X/d_t$ 直接对应 |
| 窗口曲线距离 | $l_1$ | 表示两条 Betti curves 之间的面积差 |

这些是为了让实验统一并控制工作量的默认设置，不需要把它们描述为唯一正确的选择。

### 9.2 步骤 1 下载和读取数据

1. 从 SNAP 下载 `soc-sign-bitcoinalpha.csv.gz`；
2. 读入四列：`SOURCE`、`TARGET`、`RATING`、`TIME`；
3. 把 `TIME` 转换为日期；
4. 按时间排序；
5. 检查缺失值、重复行和评分范围。

**这一阶段的输入与输出：**

- 输入是一张原始评分记录表；
- 输出是一张清洗后的数据表、数据时间范围和基础统计；
- 这一阶段不进行拓扑计算，只确认数据可以正确使用。

**Presentation 中可以说：**

> We begin with a weighted, signed, directed and temporal network. Each row records who rated whom, the rating value, and the time of the rating.

### 9.3 步骤 2 描述性分析

至少计算并展示：

- 每个时间窗口的评分数量；
- 每个窗口的活跃用户数量；
- 正评分与负评分比例；
- 每个用户给出评分数量的分布。

这一步可以帮助判断半年窗口是否合适，也能为最后解释 TDA 结果提供背景。

**这一步的具体输出：**

- 每个窗口的用户数、评分数和正负评分比例表；
- 评分数量随时间变化图；
- 活跃用户数量随时间变化图；
- 早期、中期和后期三个代表窗口。

这些普通统计量是解释 TDA 结果的背景。例如，如果后期的 $\beta_0$ 下降更快，同时后期的边数也明显增加，就应说明拓扑变化可能与网络增长同时发生。

**Presentation 中可以说：**

> We divide the temporal network into six-month windows. Each window becomes one data object that will later be represented by topological features.

### 9.4 步骤 3 构造每个窗口的数据集

对每个时间窗口：

1. 如果同一 `SOURCE-TARGET` 有多条记录，保留时间最晚的一条；
2. 保留 `RATING > 0` 的边；
3. 为每个用户建立 $T_d(u)$；
4. 删除 $|T_d(u)|<3$ 的用户；
5. 如果用户仍然很多，保留最活跃的前 200 名；
6. 计算两两 Jaccard distance matrix。

如果一个窗口保留了 $n$ 位用户，这一步会产生一个 $n\times n$ 的对称距离矩阵。对角线为 0，其余元素位于 0 到 1 之间。

例如，假设

\[
T(A)=\{D,E,F\},\qquad T(B)=\{D,E,G\}.
\]

两人的交集有 2 个元素，并集有 4 个元素，因此

\[
d(A,B)=1-\frac{2}{4}=0.5.
\]

这表示 A 和 B 的信任对象有一定重合，但并不完全相同。A 和 B 距离较小不一定表示他们互相信任，只表示他们信任的对象相似。

**这一步的具体输出：**

- 每个窗口的活跃用户列表；
- 每个窗口的 Jaccard distance matrix；
- 每个窗口的距离值 histogram。

**Presentation 中可以说：**

> We do not measure whether two users directly trust each other. We measure whether they trust similar sets of people. Jaccard distance turns this behavioural similarity into a pseudometric.

### 9.5 步骤 4 层次聚类

对每个窗口的 Jaccard distance matrix 运行 single-linkage clustering，并保存 dendrogram。

不必在报告中展示所有窗口的 dendrogram。可以选择早期、中期和后期各一个窗口进行展示，其余结果用表格总结。

**怎样读取用户 dendrogram：**

- 两个用户或群体在较低高度合并，说明它们的信任行为比较相似；
- 在较高高度才合并，说明它们差异较大；
- 一个长期独立的大分支表示一种较独特的信任行为模式；
- dendrogram 一定能够计算出来，但是否存在清楚的大分支要由实际结果决定。

**这一步的具体输出：**

- 每个窗口的用户 dendrogram；
- 早期、中期和后期三个代表性 dendrogram；
- 不同尺度下 connected components 的数量。

**Presentation 中可以说：**

> Single linkage shows how behavioural groups merge as the distance threshold increases. It is the dendrogram representation of the multi-scale partitions $X/d_t$.

### 9.6 步骤 5 Persistent homology

对每个窗口计算：

- $H_0$ persistence barcode；
- $H_1$ persistence barcode；
- $\beta_0(t)$ 和 $\beta_1(t)$。

优先使用课程中使用的软件。如果课程代码基于 StableBear，可以沿用 StableBear；如果需要 Python 实现，也可以使用 `ripser`、`gudhi`、`scipy`、`numpy` 和 `pandas`。

计算过程可以理解为：

1. 每个用户成为一个 vertex；
2. 当两个用户的距离不超过 $t$ 时加入 edge；
3. 当一组用户两两之间都有 edge 时加入相应的 simplex；
4. 让 $t$ 从 0 增加到 1，形成 Vietoris--Rips filtration；
5. 在 $\mathbb F_2$ 上计算每个尺度的 $H_0$ 和 $H_1$；
6. 把同调类的出生与死亡记录为 barcode。

**$H_0$ 的解释：**

- 每一条 bar 代表一个 connected component；
- bar 结束表示该 component 与其他 component 合并；
- 很多较长的 $H_0$ bars 表示用户行为比较分散；
- $\beta_0(t)$ 下降较快表示用户容易在较小尺度形成大群体。

**$H_1$ 的解释：**

- 每一条 bar 代表一个在某段尺度存在的一维环；
- bar 越长，说明该结构对尺度变化越稳定；
- 这里的环属于信任行为相似性空间，不是原始网络中的直接交易环；
- 如果 $H_1$ 很弱，也可以报告该表示主要反映 connected components。

**Presentation 中可以说：**

> Persistent homology records which structures survive across many distance thresholds. We mainly use $H_0$ for group merging and $H_1$ for loop-like structure in behavioural similarity.

### 9.7 步骤 6 比较时间窗口

1. 在统一的尺度网格 $t\in[0,1]$ 上计算 Betti curves；
2. 计算窗口之间的 $l_1$ 距离；
3. 得到时间窗口距离矩阵；
4. 对时间窗口运行 hierarchical clustering；
5. 结合原始网络的用户数、边数和正负评分比例解释结果。

**这一步的具体输出：**

- 时间窗口之间的距离矩阵；
- 窗口距离 heatmap；
- 最终的时间窗口 dendrogram；
- 对可能的早期、中期和后期阶段的总结。

**怎样读取最终结果：**

- 两个窗口的 Betti curves 越接近，$l_1$ 距离越小，在 dendrogram 中越早合并；
- 如果早期窗口聚在一起、后期窗口聚在一起，可以描述为不同时间阶段；
- 如果某个窗口单独形成长分支，可以说其拓扑结构与其他窗口差异较大；
- 如果所有窗口都很相似，可以说在当前数据表示下，网络拓扑随时间相对稳定。

**Presentation 中可以说：**

> Finally, each time window is represented by its Betti curves. We compare these curves and cluster the windows, so the final dendrogram describes temporal regimes rather than individual users.

## 10 最少需要生成的结果

为了形成一份完整但不过度复杂的课程报告，建议至少准备以下结果：

1. **数据随时间的分布图**：每个窗口的评分数和活跃用户数；
2. **数据统计表**：每个窗口的用户数、边数、正评分比例；
3. **三个用户 dendrogram**：早期、中期和后期各一个；
4. **三个 $H_0/H_1$ barcode 图**：与上述三个窗口对应；
5. **Betti curve 比较图**：把若干窗口的 $\beta_0$ 和 $\beta_1$ 画在一起；
6. **时间窗口距离热图**；
7. **时间窗口 dendrogram**；
8. **结果解释表**：指出哪些窗口相似、哪些窗口不同，以及可能原因。

这些结果已经足以构成一个完整的课程项目。不需要为了增加篇幅额外加入复杂机器学习模型。

## 11 结果可以怎样解释

项目开始前不应假定一定会得到明显的阶段划分。可能出现的结果包括：

- 早期窗口的 $\beta_0$ 下降较慢，表示用户的信任行为较分散；
- 后期窗口的 $\beta_0$ 下降较快，表示用户更容易形成较大的相似行为群体；
- 某些窗口出现较多或较长的 $H_1$ bars，表示信任行为相似性结构较复杂；
- 所有窗口的 $H_1$ 都很弱，说明这个数据表示主要体现 cluster structure，而没有明显的一维循环结构；
- 时间窗口 dendrogram 主要按早期和后期分组；
- 或者没有清楚的时间分组，说明这种 Jaccard 表示下的拓扑结构相对稳定。

“没有明显的 $H_1$”或者“没有清楚的时间阶段”并不意味着项目失败。只要方法正确、结果完整，并能够解释限制，就仍然是有效的数据分析。

## 12 方法限制

报告中可以简要说明以下限制：

1. Jaccard 距离只使用评分是否为正，没有使用 +1 到 +10 的强度差异；
2. 主分析没有把负评分纳入距离；
3. 用户筛选阈值和时间窗口宽度会影响结果；
4. 不同时间窗口中的用户不一定相同；
5. $H_1$ 的计算可能需要限制用户数量；
6. Jaccard 空间中的环不能直接解释为原始信任网络中的交易环；
7. 这是探索性课程项目，不试图预测 Bitcoin 价格，也不对用户是否欺诈作因果结论。

## 13 建议的报告结构

### 1 Introduction

- Bitcoin Alpha 平台为什么需要信任评分；
- 为什么研究信任网络的时间变化；
- 项目的研究问题和主要步骤。

### 2 Mathematical Background

- distance、metric 与 pseudometric；
- Jaccard distance；
- hierarchical clustering 与 dendrogram；
- Vietoris--Rips filtration；
- $H_0$、$H_1$、barcode 与 Betti curve。

### 3 Dataset and Preprocessing

- 数据来源、规模和字段；
- 时间窗口；
- 正评分定义；
- 活跃用户筛选；
- Jaccard distance matrix 的构造。

### 4 Methods

- single-linkage clustering；
- persistent homology；
- Betti curve；
- 时间窗口之间的 $l_1$ 距离和聚类。

### 5 Results

- 描述性统计；
- 代表性窗口的 dendrogram 和 barcode；
- Betti curves；
- 时间窗口距离热图和 dendrogram。

### 6 Discussion and Limitations

- 主要观察；
- 与原始网络统计的关系；
- 方法限制；
- 哪些结果不能被过度解释。

### 7 Conclusion

- 简要回答研究问题；
- 总结 TDA 是否显示出时间结构差异。

## 14 口头答辩需要能够解释的内容

为了通过答辩，组员至少应能够用自己的话回答：

1. 为什么原始数据是有向、带权、带符号、带时间的网络？
2. $T_d(u)$ 表示什么？
3. Jaccard distance 为什么适合比较两个用户？
4. 为什么得到的是 pseudometric，而不一定是 metric？
5. dendrogram 中的合并高度表示什么？
6. Vietoris--Rips complex 如何从距离矩阵构造？
7. $H_0$ 和 $H_1$ 分别表示什么？
8. barcode 中 bar 的长度是什么意思？
9. 为什么要比较不同时间窗口的 Betti curves？
10. 项目最重要的限制是什么？

## 15 最小通过版本检查表

- [ ] 成功读取并描述 Bitcoin Alpha 数据；
- [ ] 完成合理的时间窗口划分；
- [ ] 明确定义 $X_d$、$T_d(u)$ 和 Jaccard distance；
- [ ] 每个窗口得到一个距离矩阵；
- [ ] 至少展示若干 single-linkage dendrogram；
- [ ] 计算并展示 $H_0$ barcode；
- [ ] 尝试计算 $H_1$，即使结果较弱也如实报告；
- [ ] 画出并比较 Betti curves；
- [ ] 得到时间窗口距离矩阵和 dendrogram；
- [ ] 对结果进行不过度的解释；
- [ ] 报告中包含课程数学背景、数据、方法、结果和限制；
- [ ] 每位组员能够解释核心数学概念。

## 16 可选扩展

以下内容不是通过项目所必需的，只有在主分析已经完成后再考虑：

- 比较半年窗口和季度窗口；
- 比较 single linkage 和 average linkage；
- 使用评分强度构造第二种距离；
- 单独研究负评分；
- 对活跃用户阈值做简单敏感性分析；
- 将 TDA 结果与网络 density、平均 degree 等普通网络统计量进行比较。

项目应优先保证主线完整，而不是同时实现很多扩展。

## 17 方法与课程知识最终总结

### 17.1 项目到底需要哪些方法

完成主分析只需要以下七类方法：

1. **普通数据处理**：读取 CSV、转换时间、删除重复记录和划分时间窗口；
2. **描述性统计**：统计每个窗口的用户数、评分数和正负评分比例；
3. **集合表示**：把每位用户表示为其正向信任对象的集合；
4. **Jaccard pseudometric**：把用户行为相似性转换成距离矩阵；
5. **Single-linkage clustering**：观察用户群体怎样随尺度合并；
6. **Persistent homology**：计算 $H_0$、$H_1$、barcodes 和 Betti curves；
7. **曲线距离与第二次聚类**：比较 Betti curves，并对时间窗口进行聚类。

其中真正属于课程核心的部分是第 4--7 项。第 1--3 项是把原始网络整理成课程能够分析的输入。

### 17.2 方法、课程内容和结果之间的对应

| 实验方法 | 课程知识 | 产生的结果 | 最简单的解释 |
|---|---|---|---|
| 时间窗口划分 | Points 2--3 的输入 pipeline | 一组 $X_d$ | 每个窗口是一个待分析对象 |
| Jaccard distance | Points 14--15、50、52、54 | 用户距离矩阵 | 信任对象越相似，用户距离越小 |
| Single linkage | Points 25、117、122 | 用户 dendrogram | 观察用户群体随尺度合并 |
| Vietoris--Rips filtration | Points 187--193 | 一系列 simplicial complexes | 把距离关系变成可计算同调的几何对象 |
| $H_0$ | Points 198--208 | $H_0$ barcode、$\beta_0$ | 描述 connected components 和群体合并 |
| $H_1$ | Homology 部分 | $H_1$ barcode、$\beta_1$ | 描述相似性空间中的一维环 |
| Persistent bar decomposition | Points 135、158、160、164 | Barcodes | 记录结构在哪个尺度出生和死亡 |
| Betti curve 的 $l_1$ 距离 | Points 55、57、60、61 | 窗口距离矩阵 | 曲线面积差越小，窗口越相似 |
| 对窗口做 hierarchical clustering | Point 117 | 时间窗口 dendrogram | 得到可能的时间阶段 |

### 17.3 整个项目的核心逻辑

Presentation 时只要始终围绕下面四句话，就不容易讲乱：

1. 我们先把每个时间窗口中的用户变成一个 pseudometric space；
2. 距离衡量的不是用户是否直接互相信任，而是他们是否信任相似的人；
3. Persistent homology 描述这个用户空间在多个距离尺度下的 connected components 和环；
4. 比较不同窗口的 Betti curves，就能判断网络的拓扑结构是否随时间变化。

也可以把最终结论写成以下形式：

> The Bitcoin Alpha network was divided into temporal windows. For each window, users were represented by their positive-trust neighbourhoods and compared using Jaccard distance. Persistent homology summarised the multi-scale structure of each window, and distances between Betti curves were used to identify temporal regimes.

## 18 Presentation 准备提纲

### 18.1 建议的 10--15 分钟结构

| 时间 | 内容 | 需要回答的问题 |
|---:|---|---|
| 1 分钟 | 项目动机和研究问题 | 我们为什么研究 Bitcoin Alpha 信任网络？ |
| 1 分钟 | 数据集 | 节点、边、评分和时间戳分别是什么？ |
| 2 分钟 | 数据表示 | 时间窗口、$T_d(u)$ 和 Jaccard distance 怎样定义？ |
| 2 分钟 | 课程背景 | pseudometric、dendrogram 和 VR filtration 是什么？ |
| 3 分钟 | 主要结果 | $H_0/H_1$ barcode 和 Betti curves 显示什么？ |
| 2 分钟 | 时间窗口比较 | heatmap 和最终 dendrogram 怎样解释？ |
| 1 分钟 | 限制与结论 | 方法忽略了什么，最终回答是什么？ |

### 18.2 推荐的讲述顺序

可以按下面的逻辑连续讲：

> Bitcoin Alpha users rate one another from distrust to trust. Because the ratings have timestamps, we can study how the structure changes over time. We divide the data into six-month windows. Within each window, every user is represented by the set of people they rate positively. Jaccard distance measures whether two users trust similar people. This gives a pseudometric space for each window. We then apply single-linkage clustering and construct a Vietoris--Rips filtration. The resulting $H_0$ barcode describes how user groups merge, while $H_1$ describes possible loop-like structures in behavioural similarity. Finally, we compare the Betti curves of the windows and cluster the windows themselves to identify possible temporal regimes.

### 18.3 展示结果时应说什么

展示 $\beta_0$ 曲线时：

> A faster decrease means that users merge into larger behavioural groups at smaller distance thresholds. A slower decrease means the trust behaviour is more fragmented.

展示 $H_1$ 时：

> These loops are not direct trust cycles. They are loops in the space defined by similarity of users' positive-trust neighbourhoods.

展示最终时间窗口 dendrogram 时：

> Windows that merge at a low height have similar Betti curves and therefore similar multi-scale topology. Large branches can be interpreted as different temporal regimes.

### 18.4 常见问题的简短回答

**为什么使用 Jaccard distance？**  
因为每位用户被表示为一个有限集合，而 Jaccard distance 能直接衡量两个集合的重合程度，并且它满足课程需要的 metric 条件。

**为什么它只是 pseudometric？**  
因为两位不同用户可能信任完全相同的一组人，于是他们之间的距离为 0。

**为什么主要看 $H_0$？**  
因为 $H_0$ 直接描述 connected components，与用户群体和 single-linkage dendrogram 很容易对应。

**$H_1$ 中的环代表互相信任吗？**  
不一定。它表示行为相似性空间中的环，而不是原始有向网络中的直接信任环。

**为什么要比较 Betti curves？**  
因为不同窗口的用户集合可能不同，但每个窗口都可以产生定义在相同尺度区间上的曲线，因此曲线提供统一的比较对象。

**最终 dendrogram 中聚类的是什么？**  
最终聚类的是时间窗口，不是用户。用户 dendrogram 描述窗口内部结构，时间窗口 dendrogram 描述不同时期之间的相似性。

## 19 参考资料

1. Stanford Network Analysis Project. [Bitcoin Alpha trust weighted signed network](https://snap.stanford.edu/data/soc-sign-bitcoin-alpha.html).
2. S. Kumar, F. Spezzano, V. S. Subrahmanian, and C. Faloutsos. *Edge Weight Prediction in Weighted Signed Networks*. IEEE International Conference on Data Mining, 2016.
3. [SF2956 Topological Data Analysis course notes](../NotesMC.pdf).

