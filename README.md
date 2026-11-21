# Pain Radar OS · GitHub Actions 云端多维数据采集器

本目录是一个**零外部依赖、开箱即用**的 GitHub Actions 独立仓库模板，专门解决本机 ISP 住宅代理限制 `.gov` 域名（导致 `ClinicalTrials.gov` 与 `openFDA` 无法直连）的问题。

它在 GitHub 官方美国云端机房（Ubuntu Runner）每天自动定时采集 **6 大维度一手公开医药数据**，生成统一的 `data/pain_radar_bundle.json.gz`（同时保留明文 `.json` 与运行摘要 `sync_summary.md`），供本机 `PAIN RADAR OS` 一键拉取或拖入入库。

---

## 采集覆盖的 6 大官方免费源（零商业 API 采购）

1. **ClinicalTrials.gov API v2**：6 个镇痛核心切片（NaV1.8/钠通道、外周KOR/偏向MOR、术后长效局麻/改良NSAID、CGRP偏头痛、含中国中心镇痛临床、提前终止 TERMINATED 排雷），自动过滤针灸/按摩/纯器械噪音。
2. **openFDA Drugs@FDA**：25 个镇痛标杆活性成分的美国 FDA 批准记录（NDA/BLA/ANDA、剂型、规格、持有人、首次/最近批准日）。
3. **openFDA FAERS**：14 个经典镇痛药真实世界 Top 12 不良反应与报告占比（**2.2 类改良型新药立项痛点金矿**）。
4. **Open Targets Platform GraphQL + RCSB PDB**：14 个核心镇痛靶点基因（`SCN10A`/`SCN9A`/`OPRK1`/`OPRM1`/`CALCRL`/`TRPV1`/`NTRK1` 等）的人类遗传学镇痛评分、小分子/抗体成药性、已知临床管线数、安全警示与实验三维晶体结构数量。
5. **ChEMBL REST API (EMBL-EBI)**：镇痛核心分子及靶点关联管线的理化性质（MW / ALogP / Ro5）、最高临床分期、首次获批年份与作用机制（MoA）。
6. **Europe PMC REST API**：2022–2027 镇痛顶刊（NEJM / Lancet / JAMA / Anesthesiology / Pain 等）文献元数据、被引次数与自动提取的关联 `NCT` 试验号。

---

## 5 分钟上线指南（公开或私有仓库均可）

### 第 1 步：在 GitHub 新建一个空仓库
1. 登录 GitHub，点击右上角 **New repository**；
2. 仓库名填例如 `pain-radar-cloud-sync`，可见性可选 **Private（私有，更私密）** 或 **Public（公开，拉取免 Token）**；
3. 不要勾选 README/gitignore，点击 **Create repository**。

### 第 2 步：把本目录推送到 GitHub
在本目录（`D:\AI\Antigravity\Projects\Pharma\cloud_sync`）打开 PowerShell，执行：
```powershell
git init -b main
git add .
git commit -m "init: Pain Radar OS multi-dimensional cloud sync"
git remote add origin https://github.com/<你的GitHub用户名>/pain-radar-cloud-sync.git
git push -u origin main
```

### 第 3 步：手动触发首次全量采集
1. 打开你的 GitHub 仓库 → 点击顶部 **Actions** 选项卡；
2. 左侧选中 **Pain Radar OS · Multi-Dimensional Cloud Sync** → 点击右侧 **Run workflow**；
3. `mode` 选 **`full`**，点击绿色按钮运行（约 2–4 分钟完成）；
4. 跑完后，点进该次运行可直接看到 **Job Summary** 各维度入库数量，且仓库 `data/` 目录下会自动生成 `pain_radar_bundle.json.gz`（也可在运行页下方 Artifacts 直接下载 zip）。

### 第 4 步：在本机 PAIN RADAR OS 一键同步入库
打开本地雷达（`http://127.0.0.1:3722`）→ 进入 **「云端同步与导入」** 或 **「定时与配置」**：
- **若为 Public 仓库**：填入 Raw 地址即可一键拉取（无需 Token）：
  ```
  https://raw.githubusercontent.com/<你的用户名>/pain-radar-cloud-sync/main/data/pain_radar_bundle.json.gz
  ```
- **若为 Private 仓库**：填入同一个 Raw 地址 + 一个只读 Fine-grained PAT Token（仅勾选该仓库 `Contents: Read-only` 权限）即可一键拉取；或者直接把 GitHub Actions 页面下载的 `pain_radar_bundle.json.gz` 拖进导入框！

> 此后每天北京时间 **06:30** 云端自动增量更新（每周日自动全量回扫），你早上打开本地雷达点一下「一键拉取云端包」即可完成 6 大维度全量刷新！
