# -*- coding: utf-8 -*-
"""
glycan_etl — 糖类数据库 PDF→结构化入库 解析引擎包
=================================================
为「AI 辅助推断多糖结构平台」提供"已知结构-谱图"真值参照数据的 ETL 引擎。

子模块:
  - core        核心解析入库引擎（原 glycan_etl_v3.py，段落式正文解析）
  - table_parser 表格式文献解析器（分子量/归属表/正文键连，补齐 core 缺口）
  - embeddings  bge-m3 向量化（为 nmr_shifts_1d 生成 1024 维 embedding）
"""
__version__ = "3.1.0"

__all__ = ["core", "table_parser", "embeddings"]
