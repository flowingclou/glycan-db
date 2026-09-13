-- ============================================================================
-- 糖类数据库 · V4 迁移（向量检索能力 + 多糖性质溯源 + 检索索引）
-- 前置: 需先执行 v001_mvp.sql → v002_extend.sql → v003_qc.sql
-- 数据库: PostgreSQL 16+ + pgvector
-- ============================================================================
-- 背景: 仓库早期版本把 pgvector 能力只写进了 README 与 embeddings.py,
--       却没有任何 DDL 创建 nmr_shifts_1d.embedding 列 —— 于是
--       `python3 glycan_etl/embeddings.py` 必然报 "column embedding does not
--       exist", README §5 的近邻检索 SQL 也无法执行。本迁移补齐该能力。
-- ============================================================================

-- 0. 前置扩展（embedding 列需要 vector；glycoct_hash 需要 pgcrypto）
CREATE EXTENSION IF NOT EXISTS vector;
CREATE EXTENSION IF NOT EXISTS pgcrypto;

-- ----------------------------------------------------------------------------
-- 1. 一维位移记录的向量列（BGE-M3 = 1024 维）
--    若换用其他 embedding 模型，维度需与这里保持一致（或另开一列）。
-- ----------------------------------------------------------------------------
ALTER TABLE nmr_shifts_1d
    ADD COLUMN IF NOT EXISTS embedding vector(1024);

-- 余弦近邻索引（HNSW 在中小规模真值库上召回与延迟都优于 IVFFlat）
CREATE INDEX IF NOT EXISTS idx_nmr_shifts_1d_embedding_hnsw
    ON nmr_shifts_1d USING hnsw (embedding vector_cosine_ops);

-- ----------------------------------------------------------------------------
-- 2. polysaccharide_props 溯源列
--    原表只有 sugar_id: 同一结构被多篇文献报道时无法区分数据来源,
--    重复导入也无法按 (结构, 文献) 幂等覆盖。
-- ----------------------------------------------------------------------------
ALTER TABLE polysaccharide_props
    ADD COLUMN IF NOT EXISTS source_id BIGINT REFERENCES literature(source_id);

CREATE INDEX IF NOT EXISTS idx_poly_props_sugar ON polysaccharide_props(sugar_id);

-- ----------------------------------------------------------------------------
-- 3. 常用检索索引（AI 平台回查「候选结构 vs 真值」时的主要过滤条件）
-- ----------------------------------------------------------------------------
CREATE INDEX IF NOT EXISTS idx_sugars_iupac        ON sugars(iupac_short);
CREATE INDEX IF NOT EXISTS idx_exp_nmr_type        ON nmr_experiments(nmr_type);
CREATE INDEX IF NOT EXISTS idx_exp_solvent         ON nmr_experiments(solvent);
CREATE INDEX IF NOT EXISTS idx_shift_nucleus_ppm   ON nmr_shifts_1d(nucleus, shift_ppm);
CREATE INDEX IF NOT EXISTS idx_shift_anomeric      ON nmr_shifts_1d(is_anomeric)
    WHERE is_anomeric;
CREATE INDEX IF NOT EXISTS idx_corr_2d_linkage     ON nmr_correlations_2d(linkage_evidence)
    WHERE linkage_evidence;

-- ----------------------------------------------------------------------------
-- 4. 幂等约束：同一 (结构, 文献, 谱类型, 溶剂) 只应有一条实验头
--    对应 glycan_etl/core.py insert_record() 的「先删后写」策略;
--    该唯一约束同时防止并发/中断重跑产生重复实验头。
-- ----------------------------------------------------------------------------
CREATE UNIQUE INDEX IF NOT EXISTS uq_exp_sugar_source_type_solvent
    ON nmr_experiments(sugar_id, source_id, nmr_type, solvent);

-- ----------------------------------------------------------------------------
-- 5. 验证
-- ----------------------------------------------------------------------------
-- 5.1 确认向量列存在:
-- SELECT column_name, data_type FROM information_schema.columns
-- WHERE table_name = 'nmr_shifts_1d' AND column_name = 'embedding';
--
-- 5.2 待向量化条数（embeddings.py 的处理对象）:
-- SELECT count(*) FROM nmr_shifts_1d WHERE embedding IS NULL;
