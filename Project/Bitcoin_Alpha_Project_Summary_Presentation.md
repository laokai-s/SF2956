# Bitcoin Alpha TDA Project Summary Presentation

> 本文档用于准备课程最后一周的 15--20 分钟 project summary。此时项目仍处于 proposed project 阶段，因此应当介绍数据集、研究问题、与课程的联系、可能采用的方法、预期输出和可行性，不应把尚未计算的内容说成真实结果。

## 1 这次 summary 需要达到什么目的

根据课程要求，这次 presentation 的核心不是汇报一个已经完成的 final project，而是让老师确认：

1. 所选数据集是真实、可获得并且适合分析的；
2. 研究问题足够清楚；
3. 项目确实使用课程介绍的工具；
4. 方法路线可以在课程项目的时间和计算条件内完成；
5. 小组成员理解准备使用的数学对象；
6. 即使最终结果不符合预期，项目仍然能形成完整的数据分析。

这次 summary 最重要的结论应当是：

> Bitcoin Alpha 是一个带权、带符号、有向并且带时间戳的信任网络。我们计划把网络划分为时间窗口，用 Jaccard pseudometric 表示用户信任行为的相似性，再使用 single-linkage clustering、Vietoris--Rips filtration 和 persistent homology 比较不同时间窗口的拓扑结构。

## 2 项目的一页摘要

### Proposed title

**Topological Analysis of Temporal Trust Structure in the Bitcoin Alpha Network**

### Motivation

Bitcoin Alpha 平台中的用户是匿名的，因此平台通过用户之间的正负评分记录信任与不信任。由于评分带有时间戳，我们可以研究用户信任行为的群体结构是否随时间发生变化。

### Research question

> How does the multi-scale topological structure of user trust behaviour in the Bitcoin Alpha network change over time?

### Dataset

- 数据来源：Stanford Network Analysis Project；
- 节点：3,783 位用户；
- 边：24,186 条评分；
- 评分范围：-10 到 +10；
- 正评分比例：约 93%；
- 每条记录包含 SOURCE、TARGET、RATING 和 TIME；
- 网络具有 directed、weighted、signed 和 temporal 四种性质。

### Proposed representation

把数据划分为半年时间窗口。对窗口 $d$ 中的用户 $u$，定义

\[
T_d(u)=\{v\mid u\text{ 在窗口 }d\text{ 中对 }v\text{ 给出正评分}\}.
\]

用 Jaccard distance 比较两位用户：

\[
d_d(u,v)=1-\frac{|T_d(u)\cap T_d(v)|}{|T_d(u)\cup T_d(v)|}.
\]

### Proposed TDA pipeline

\[
\text{时间窗口}
\longrightarrow
\text{Jaccard 距离}
\longrightarrow
\text{single-linkage dendrogram}
\longrightarrow
\text{VR filtration}
\longrightarrow
H_0,H_1
\longrightarrow
\text{Betti curves}
\longrightarrow
\text{窗口聚类}.
\]

### Expected outputs

- 每个时间窗口的描述性统计；
- 代表性窗口的用户 dendrogram；
- $H_0$ 和 $H_1$ barcodes；
- $\beta_0$ 和 $\beta_1$ curves；
- 时间窗口之间的距离 heatmap；
- 时间窗口 dendrogram；
- 对可能的早期、中期和后期阶段的解释。

## 3 建议的 presentation 时间安排

| 时间 | 页面 | 内容 |
|---:|---:|---|
| 0:00--1:00 | 1 | 标题和一句话项目概述 |
| 1:00--2:30 | 2 | 背景与动机 |
| 2:30--4:30 | 3 | Bitcoin Alpha 数据集 |
| 4:30--5:30 | 4 | 研究问题 |
| 5:30--7:00 | 5 | 时间窗口与数据预处理 |
| 7:00--9:00 | 6 | 用户集合表示与 Jaccard distance |
| 9:00--10:30 | 7 | 与课程 pipeline 的对应 |
| 10:30--12:00 | 8 | Single linkage 与 $H_0$ |
| 12:00--13:30 | 9 | VR filtration、$H_1$ 与 barcode |
| 13:30--15:00 | 10 | 不同时间窗口的比较 |
| 15:00--16:30 | 11 | 预期输出和解释方式 |
| 16:30--17:30 | 12 | 范围、可行性和限制 |
| 17:30--18:30 | 13 | 结论与希望获得的反馈 |

如果实际只有 15 分钟，可以缩短第 2、5、9 和 12 页，不要删除数据、Jaccard distance、$H_0$ 和最终窗口比较。

## 4 Slide 1 项目标题

### Slide title

**Topological Analysis of Temporal Trust Structure in the Bitcoin Alpha Network**

### Slide 上放什么

- 课程：SF2956 Topological Data Analysis；
- 小组成员姓名；
- Bitcoin Alpha trust network；
- 一句话目标：

> Detecting temporal changes in user trust behaviour using multi-scale topological features.

### 讲解重点

开场不需要解释数学，只要说明研究对象和目标：

> Our project studies how trust behaviour in the Bitcoin Alpha network changes over time. We plan to represent users through their positive-trust neighbourhoods and compare temporal windows using tools from topological data analysis.

## 5 Slide 2 背景与动机

### Slide title

**Why Study Trust in Bitcoin Alpha?**

### Slide 上放什么

- Bitcoin 平台用户可以保持匿名；
- 匿名交易需要 reputation 或 trust information；
- 用户可以给其他用户 -10 到 +10 的评分；
- 评分形成一个随时间变化的信任网络；
- 普通统计量能说明网络大小，但不一定能描述多尺度群体结构。

### 建议图形

画一个很小的有向网络示意图：

- 节点代表用户；
- 正边代表 trust；
- 负边代表 distrust；
- 边旁边标注评分；
- 不需要在此页展示真实的大型网络。

### 讲解重点

> A simple count of users or ratings does not describe how groups of similar trust behaviour are organised. We therefore want to study structure across multiple distance scales.

不要把项目描述为 Bitcoin 价格预测，也不要声称能够识别真实诈骗者。

## 6 Slide 3 数据集

### Slide title

**Bitcoin Alpha Dataset**

### Slide 上放什么

| Property | Value |
|---|---:|
| Users | 3,783 |
| Ratings | 24,186 |
| Rating range | -10 to +10 |
| Positive ratings | approximately 93% |
| Network type | directed, weighted, signed, temporal |

每条数据记录：

\[
(\text{SOURCE},\text{TARGET},\text{RATING},\text{TIME}).
\]

数据来源：

[Stanford SNAP Bitcoin Alpha dataset](https://snap.stanford.edu/data/soc-sign-bitcoin-alpha.html)

### 讲解重点

解释四个形容词：

- directed：A 对 B 的评分与 B 对 A 的评分不同；
- weighted：评分有强弱；
- signed：评分可以为正或负；
- temporal：每条评分都有时间戳。

### 建议说法

> The temporal information is essential for our project because it allows us to compare the topology of the network at different stages.

## 7 Slide 4 研究问题

### Slide title

**Research Question**

### 主问题

> How does the multi-scale topological structure of user trust behaviour change over time?

### 可以放在下面的子问题

1. Do users form groups with similar positive-trust behaviour?
2. Do these groups merge differently in early and late time windows?
3. Are the $H_0$ and $H_1$ signatures different across time?
4. Can temporal windows be grouped into different regimes?

### 讲解重点

这里应当强调“exploratory analysis”，而不是预先保证一定找到 cluster。

> We are not assuming that clear regimes must exist. Our goal is to use the course pipeline to test whether the multi-scale structure changes over time.

## 8 Slide 5 时间窗口和预处理

### Slide title

**Temporal Windows and Preprocessing**

### Slide 上放什么

主分析计划：

1. 把 Unix timestamps 转换为日期；
2. 把完整数据划分为半年窗口；
3. 同一 SOURCE--TARGET 在一个窗口内有多次评分时，保留最后一次；
4. 评分大于 0 定义为 positive trust；
5. 负评分用于描述性统计，不进入主距离；
6. 只保留至少给出 3 个正评分的用户；
7. 每个窗口最多保留 200 位最活跃用户。

### 为什么这样处理

- 半年窗口比月度窗口更不容易过于稀疏；
- 活跃用户筛选避免大量空 trust sets；
- 固定最大用户数控制 $H_1$ 的计算量；
- 所有窗口使用相同规则，保证比较具有一致性。

### 建议图形

画一条时间轴，并标出：

\[
d_1,d_2,d_3,\ldots,d_m
\]

每个 $d_i$ 表示一个半年窗口。

### 讲解重点

> Each temporal window is treated as one object in the course notation $I=\{X_d\mid d\in D\}$.

## 9 Slide 6 用户表示和 Jaccard distance

### Slide title

**From Trust Behaviour to a Pseudometric**

### Slide 上放什么

对窗口 $d$ 中的用户 $u$：

\[
T_d(u)=\{v\mid u\text{ positively rates }v\}.
\]

对两位用户：

\[
d_d(u,v)
=1-\frac{|T_d(u)\cap T_d(v)|}{|T_d(u)\cup T_d(v)|}.
\]

### 一个小例子

\[
T(A)=\{D,E,F\},\qquad T(B)=\{D,E,G\}.
\]

\[
d(A,B)=1-\frac{2}{4}=0.5.
\]

### 如何解释

- 距离 0：两位用户信任完全相同的一组人；
- 距离接近 0：信任偏好很相似；
- 距离接近 1：信任对象几乎没有重合。

### 必须说清楚的区别

> A small distance does not necessarily mean that A directly trusts B. It means that A and B trust similar people.

### 为什么是 pseudometric

两位不同的用户可能拥有完全相同的正向信任集合，因此它们之间的距离可能为 0。

## 10 Slide 7 与课程内容的对应

### Slide title

**Connection to the Course Pipeline**

### Slide 上放什么

| Project step | Course concept |
|---|---|
| Temporal windows | $I=\{X_d\mid d\in D\}$ |
| Jaccard distance | metric restriction and pseudometric |
| User groups at scale $t$ | partitions $X/d_t$ |
| Multi-scale clustering | single linkage and dendrogram |
| Geometric construction | Vietoris--Rips filtration |
| Topological invariants | $H_0$, $H_1$ and barcodes |
| Window representation | Betti curves |
| Window comparison | $l_1$ distance and hierarchical clustering |

### 可以提到的讲义位置

- Input pipeline：Points 2--3；
- Pseudometrics：Points 14--15；
- Jaccard metric：Points 50、52、54；
- Hierarchical clustering：Points 25、117、122；
- Vietoris--Rips：Points 187--193；
- Homology：Points 198--208；
- Barcodes：Points 135、158、160、164；
- Piecewise-constant functions：Points 55--61。

### 讲解重点

不要逐个解释所有 Point，也不需要证明 theorem。只需要告诉老师，这个项目的每一个主要步骤都能在课程框架中找到对应对象。

## 11 Slide 8 Single linkage 和零维同调

### Slide title

**Multi-Scale Groups and \(H_0\)**

### Slide 上放什么

- 对每个时间窗口计算 single-linkage dendrogram；
- 距离尺度 $t$ 从 0 增加到 1；
- 在小尺度只合并非常相似的用户；
- 在大尺度形成较大的群体；
- $H_0$ 记录 connected components。

### 图形计划

最终项目中展示：

- 一个早期窗口 dendrogram；
- 一个中期窗口 dendrogram；
- 一个后期窗口 dendrogram；
- 对应的 $H_0$ barcodes 或 $\beta_0$ curves。

### 如何解释

- 较低高度合并：信任行为相似；
- 较高高度才合并：行为差异较大；
- $\beta_0$ 下降快：用户较快形成大群体；
- $\beta_0$ 下降慢：用户行为比较分散。

### 建议说法

> We expect $H_0$ to be the most interpretable part of the analysis because it directly describes how behavioural groups merge across distance scales.

## 12 Slide 9 Vietoris--Rips 和一维同调

### Slide title

**Vietoris--Rips Filtration and \(H_1\)**

### Slide 上放什么

在尺度 $t$：

- 用户是 vertices；
- 距离不超过 $t$ 的用户之间加入 edges；
- 两两相连的一组用户形成 simplex；
- 随着 $t$ 增加，得到 filtration；
- $H_1$ 检测一维环状结构。

### Barcode 的解释

- bar 的起点：拓扑结构出现；
- bar 的终点：拓扑结构消失；
- 较长的 bar：在较宽尺度范围内存在；
- 很短的 bar：可能只对应局部或短暂结构。

### 必须避免的错误解释

$H_1$ 中的环不是原始网络中 A 信任 B、B 信任 C、C 又信任 A 的直接环。它是由 Jaccard distance 定义的用户行为相似性空间中的环。

### 建议说法

> \(H_1\) is secondary in our project. If it is weak, that is still an informative result showing that the representation is dominated by component structure rather than persistent loops.

## 13 Slide 10 比较时间窗口

### Slide title

**Comparing Temporal Windows**

### Slide 上放什么

每个窗口得到：

\[
\beta_0^d(t),\qquad \beta_1^d(t).
\]

用 $l_1$ 距离比较窗口：

\[
D_k(d_i,d_j)
=\int_0^1
\left|
\beta_k^{d_i}(t)-\beta_k^{d_j}(t)
\right|dt.
\]

然后：

1. 构造时间窗口距离矩阵；
2. 画出 heatmap；
3. 对时间窗口运行 hierarchical clustering；
4. 得到最终的 temporal-window dendrogram。

### 两种 dendrogram 不要混淆

- 用户 dendrogram：研究一个窗口内部哪些用户相似；
- 时间窗口 dendrogram：研究哪些时期的拓扑结构相似。

### 讲解重点

> The final clustering is performed on temporal windows, not on users. Windows that have similar Betti curves merge earlier.

## 14 Slide 11 预期输出和可能解释

### Slide title

**Expected Outputs**

### 计划产生的图表

1. 每个窗口的评分数和活跃用户数；
2. 正负评分比例；
3. 代表性用户 dendrograms；
4. $H_0$ 和 $H_1$ barcodes；
5. 多个窗口的 Betti curves；
6. 窗口距离 heatmap；
7. 最终的时间窗口 dendrogram。

### 可能出现的模式

- 后期 $\beta_0$ 下降更快：信任行为更容易形成较大的相似群体；
- 早期 $\beta_0$ 下降较慢：信任行为更分散；
- 某个窗口具有更多或更长的 $H_1$ bars：其行为相似性结构更复杂；
- 早期和后期窗口形成不同分支：可能存在不同发展阶段；
- 所有窗口相似：在当前表示下，拓扑结构随时间较稳定。

### 表达方式

全部使用“possible”、“we expect to examine”或“we will test whether”。不要在尚未完成分析时说“we found”。

## 15 Slide 12 项目范围和可行性

### Slide title

**Scope and Feasibility**

### 主分析只包含

- 一种主要数据表示；
- 一种主要距离：Jaccard distance；
- 一种主要 linkage：single linkage；
- 两个同调维数：$H_0$ 和 $H_1$；
- 一种窗口比较：Betti-curve $l_1$ distance。

### 为什么算得动

- 每个窗口最多保留 200 位活跃用户；
- Jaccard 距离值限制在 0 到 1；
- 主要结果集中在 $H_0$；
- $H_1$ 只计算到维数 1；
- 不建立复杂分类模型；
- 不计算 $H_2$ 或更高维同调。

### 主要限制

- 忽略 +1 到 +10 之间的强度区别；
- 负评分不进入主距离；
- 时间窗口和活跃阈值会影响结果；
- 不同窗口中的用户不完全相同；
- Jaccard 相似性不等于用户之间直接互信；
- 项目是探索性分析，不作因果结论。

## 16 Slide 13 结论与希望获得的反馈

### Slide title

**Proposed Project and Questions**

### 总结

> We propose to study the temporal evolution of the Bitcoin Alpha trust network by representing users through positive-trust neighbourhoods, constructing a Jaccard pseudometric, computing persistent homology, and clustering temporal windows using distances between Betti curves.

### 可以向老师确认的问题

1. Is the Jaccard representation appropriate for the course project?
2. Is using six-month windows a reasonable initial choice?
3. Is it sufficient to focus mainly on $H_0$ and use $H_1$ as a secondary analysis?
4. Should the final comparison use raw or normalised Betti curves?
5. Would the instructors prefer StableBear for the persistent-homology computation?

这些问题能显示小组已经形成清楚方案，同时愿意根据课程重点调整。

## 17 可以直接使用的英文开场

> Our proposed project studies the temporal evolution of trust behaviour in the Bitcoin Alpha network. Bitcoin Alpha users rate one another from minus ten, representing strong distrust, to plus ten, representing strong trust. The dataset contains 3,783 users, 24,186 ratings, and a timestamp for every rating. We divide the data into six-month windows. Within each window, a user is represented by the set of people they rate positively. We then compare users using Jaccard distance. This gives a pseudometric space for every window. We plan to apply single-linkage clustering, construct Vietoris--Rips filtrations, and compute persistent homology in dimensions zero and one. Finally, we compare the Betti curves of the temporal windows and cluster the windows themselves. Our goal is to examine whether the multi-scale organisation of trust behaviour changes over time.

## 18 可以直接使用的英文结尾

> To summarise, the project follows the analysis pipeline introduced in the course: data objects, pseudometrics, multi-scale clustering, Vietoris--Rips filtrations, homology, barcodes, and comparison of piecewise-constant functions. The project is deliberately limited to one main distance, \(H_0\) and \(H_1\), and a moderate number of active users per window. We therefore expect it to be computationally feasible and suitable for the scope of the course. At this stage, we would especially appreciate feedback on the Jaccard representation, the choice of temporal windows, and the proposed comparison of Betti curves.

## 19 老师可能提问的问题

### Why do you need TDA instead of ordinary community detection?

Ordinary community detection usually produces a partition at one chosen setting. TDA allows us to observe how connected components and higher-order structures change across all distance thresholds.

### Why use Jaccard distance?

Each user is naturally represented by a finite set of positively rated users. Jaccard distance directly compares the overlap of these sets and is a metric on finite sets.

### Why is the user space only a pseudometric space?

Two different users may have identical positive-trust neighbourhoods, so their distance can be zero.

### What does an \(H_0\) bar mean?

It represents a connected component that exists over a range of distance thresholds. Its death time is the scale at which it merges with another component.

### What does an \(H_1\) bar mean?

It represents a loop in the user-behaviour similarity space. It is not necessarily a direct cycle of trust edges in the original network.

### How do you compare windows with different users?

We do not match individual users across windows. Each window is summarised by Betti curves defined over the same scale interval, and we compare those curves.

### What happens if there are no clear clusters?

That is still a valid result. We can report that the chosen Jaccard representation does not reveal strongly separated behavioural groups or temporal regimes.

### Why ignore negative ratings?

The main pipeline uses a clean pseudometric based on positive-trust sets. Negative ratings are retained in descriptive statistics. A signed-distance analysis would expand the scope beyond what is needed for the course project.

### What is the main expected contribution?

The contribution is an interpretable application of the course TDA pipeline to a temporal signed network, together with a comparison of the network's multi-scale structure across time.

## 20 Presentation 制作检查表

- [ ] 标题页包含项目名和组员；
- [ ] 清楚说明这是 proposed strategy，而不是已完成结果；
- [ ] 数据规模和字段准确；
- [ ] 有一张简单的信任网络示意图；
- [ ] 有一张时间窗口示意图；
- [ ] 用一个小例子解释 Jaccard distance；
- [ ] 展示项目 pipeline；
- [ ] 解释用户 dendrogram 和窗口 dendrogram 的区别；
- [ ] 解释 $H_0$、$H_1$ 和 barcode；
- [ ] 列出计划生成的结果图；
- [ ] 说明范围和计算可行性；
- [ ] 说明主要限制；
- [ ] 准备 3--5 个希望老师反馈的问题；
- [ ] 把总时长控制在 15--20 分钟；
- [ ] 每位组员都能解释自己负责的部分。

## 21 不建议在 summary 中加入的内容

- 不需要证明 Jaccard distance 的 triangle inequality；
- 不需要完整讲解 homology 的代数定义；
- 不需要展示所有课程 theorem；
- 不要同时提出很多替代距离和 linkage；
- 不要加入 Bitcoin 价格预测；
- 不要把 Costanzo 数据集也放进同一次 proposal；
- 不要声称一定会发现明显 cluster；
- 不要展示尚未计算的虚构结果；
- 不要把 $H_1$ 环解释成直接的互信关系；
- 不需要讨论复杂的 null model 或显著性检验。

## 22 准备 presentation 前需要完成的最少工作

即使 summary 主要介绍计划，最好提前完成以下小检查：

1. 确认数据文件能够下载和读取；
2. 确认四列的含义和数据规模；
3. 把时间戳成功转换为日期；
4. 画出评分数量随时间的简单图；
5. 检查半年窗口中大约有多少活跃用户；
6. 手工准备一个 Jaccard distance 的小例子；
7. 准备一张项目 pipeline 图；
8. 明确主分析只使用 Jaccard、single linkage、$H_0$ 和 $H_1$。

## 23 参考资料

1. Stanford Network Analysis Project. [Bitcoin Alpha trust weighted signed network](https://snap.stanford.edu/data/soc-sign-bitcoin-alpha.html).
2. S. Kumar, F. Spezzano, V. S. Subrahmanian, and C. Faloutsos. *Edge Weight Prediction in Weighted Signed Networks*. IEEE International Conference on Data Mining, 2016.
3. [SF2956 Topological Data Analysis course notes](../NotesMC.pdf).
4. [Detailed Bitcoin Alpha project plan](Bitcoin_Alpha_TDA项目说明.md).

