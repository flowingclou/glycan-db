-- ============================================================================
-- 糖类数据库 · V5 迁移（结构表达等级 + 残基组成式）
-- 前置: 需先执行 v001 → v002 → v003 → v004
-- 数据库: PostgreSQL 16+
-- ============================================================================
-- 背景（P0-1）:
--   GlycoCT 只能表达**确定结构**。文献里的多糖常只给出甲基化/组成信息，
--   残基之间的连接顺序未知；此时硬生成 GlycoCT 等于伪造结构。
--   因此显式记录"这条记录的结构表达到什么程度"，供下游决定能否做结构比对：
--     complete         —— 单糖/寡糖，连接明确，glycoct 为完整结构编码
--     repeat_unit      —— 单一重复单元多糖，glycoct 为一个重复单元
--     composition_only —— 仅有残基组成，glycoct 为空，看 composition
-- ============================================================================

ALTER TABLE sugars
    ADD COLUMN IF NOT EXISTS structure_level TEXT
    CHECK (structure_level IN ('complete', 'repeat_unit', 'composition_only'));

-- 残基组成式，如 "GalA6,Gal3,Ara3,Rha2,GlcA"（按数量降序）
ALTER TABLE sugars
    ADD COLUMN IF NOT EXISTS composition TEXT;

CREATE INDEX IF NOT EXISTS idx_sugars_structure_level ON sugars(structure_level);
CREATE INDEX IF NOT EXISTS idx_sugars_composition     ON sugars(composition);

-- ----------------------------------------------------------------------------
-- 回填提示: 本迁移只加列，不自动改写既有数据。
-- 既有记录的 glycoct 若是旧版伪编码（如 "RES 1b:a-lglcp-1:5|2:x"），
-- 请运行: python3 pipelines/backfill_glycoct.py --config <你的配置>
-- ----------------------------------------------------------------------------

-- 验证:
-- SELECT structure_level, count(*) FROM sugars GROUP BY 1 ORDER BY 2 DESC;
-- SELECT iupac_short, composition FROM sugars WHERE glycoct IS NULL LIMIT 10;
