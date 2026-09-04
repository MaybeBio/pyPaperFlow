# pyPaperFlow - 文献阅读自动化 🔬

<div align="center" markdown="1">

<img src="assets/logo.png" alt="pyPaperFlow Logo" width="180" />

<p><strong>面向科研工作者的自动化文献处理与知识发现平台</strong></p>

<p>批量检索、批量获取、批量解析、批量结构化，把文献真正变成可计算、可复用、可追踪的研究资产。</p>

<p>从文献检索到知识内化，把重复劳动交给流程，把关键判断留给你。</p>

![pyPaperFlow](assets/main.png)

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![PR's Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat)](http://makeapullrequest.com)
[![Workflow](https://img.shields.io/badge/Workflow-7%20Stages-0366d6)](https://github.com/MaybeBio/pyPaperFlow/blob/main/docs/Design.md)
[![PyPI version](https://img.shields.io/pypi/v/pyPaperFlow.svg?logo=pypi&logoColor=white)](https://pypi.org/project/pyPaperFlow/)
[![Python Versions](https://img.shields.io/pypi/pyversions/pyPaperFlow.svg?logo=python&logoColor=white)](https://pypi.org/project/pyPaperFlow/)
[![Downloads](https://static.pepy.tech/badge/pyPaperFlow)](https://pepy.tech/project/pyPaperFlow)

<p>
  <a href="https://github.com/MaybeBio/pyPaperFlow">GitHub 仓库</a> |
  <a href="https://github.com/MaybeBio/pyPaperFlow/blob/main/README.md">English README</a> |
  <a href="https://github.com/MaybeBio/pyPaperFlow/blob/main/README_zh.md">中文 README</a>
</p>

</div>

> **如果该项目对你有帮助, 请麻烦点一个 Star ⭐, 谢谢!**

---

## 快速开始

```bash
pip install pyPaperFlow
paperflow --help
```

详细安装与各模块用法见 [安装](installation.md) 与 [使用方法](usage/index.md)。

## 文档导航

- [简介](overview.md)
- [功能特性](features.md)
- [架构设计哲学](architecture.md)
- [安装](installation.md)
- 使用方法
    - [模块概述](usage/index.md)
    - [1. 研究起点](usage/research-start.md)
    - [2. 文献检索（及元数据抓取）](usage/search.md)
    - [3. 文献获取（及全文下载）](usage/fetch.md)
    - [4. 文献内容提取与结构化处理](usage/extraction.md)
    - [5. 其他文献数据平台的处理](usage/other-databases.md)
    - [6. 批判性阅读与知识图谱分析](usage/downstream.md)
    - [7. Reading与Coding的交点](usage/reading-coding.md)
- [全维度文献章节模块组合利用方案](combination.md)
- [测试示例](test-cases.md)
- [一个完整的文献调研示例](full-survey.md)
- [后续维护待办](todo.md)
- 参考（在 GitHub 仓库中查看）
    - [设计文档](https://github.com/MaybeBio/pyPaperFlow/blob/main/docs/Design.md)
    - [测试用例](https://github.com/MaybeBio/pyPaperFlow/blob/main/docs/Cases.md)
    - [Skills](https://github.com/MaybeBio/pyPaperFlow/blob/main/docs/Skills.md)
    - [MinerU 解析](https://github.com/MaybeBio/pyPaperFlow/blob/main/docs/mineru_parse.md)
    - [undetected 回退](https://github.com/MaybeBio/pyPaperFlow/blob/main/docs/undetected_fallback.md)
    - [PaperDB 各数据库抓取笔记](https://github.com/MaybeBio/pyPaperFlow/tree/main/docs/PaperDB)
