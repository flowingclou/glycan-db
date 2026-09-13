-- ============================================================================
-- 糖类数据库 · V6 迁移（域级结构表达等级 domain_only + 域架构字段）
-- 前置: 需先执行 v001 → v002 → v003 → v004 → v005
-- 数据库: PostgreSQL 16+
-- ============================================================================
-- 背景:
--   V5 的三级表达等级（complete / repeat_unit / composition_only）无法表达
--   "域级骨架"这一档。果胶类多糖文献（如山楂多糖 HP）常只在正文给一句：
--
--     "HP is mainly composed of a large number of HG domains and a small
--      number of RG-I domains with side chains."
--
--   这类句子**不含连接式箭头**，抓不到逐残基序列；但它明确给出了"主链是什么域、
--   侧链挂在哪个域上"，信息量远高于纯组成百分比，硬塞进 composition_only 会
--   丢掉这一层，硬生成 GlycoCT 又是伪造数据。因此新增一档 domain_only。
-- ============================================================================

-- ----------------------------------------------------------------------------
-- 1. 放宽 structure_level 的 CHECK 约束，加入 'domain_only'
--    旧约束是 ALTER TABLE ... ADD COLUMN ... CHECK (...) 自动命名的，
--    名字在不同库上可能不同，因此按 pg_constraint 动态查找后重建（幂等）。
-- ----------------------------------------------------------------------------
DO $$
DECLARE
    cname TEXT;
BEGIN
    SELECT con.conname INTO cname
    FROM pg_constraint con
    JOIN pg_class rel ON rel.oid = con.conrelid
    JOIN pg_attribute att ON att.attrelid = rel.oid AND att.attnum = ANY (con.conkey)
    WHERE rel.relname = 'sugars'
      AND con.contype = 'c'
      AND att.attname = 'structure_level'
    LIMIT 1;

    IF cname IS NOT NULL THEN
        EXECUTE format('ALTER TABLE sugars DROP CONSTRAINT %I', cname);
    END IF;

    ALTER TABLE sugars
        ADD CONSTRAINT sugars_structure_level_check
        CHECK (structure_level IN ('complete', 'repeat_unit', 'domain_only',
                                   'composition_only'));
END $$;

-- ----------------------------------------------------------------------------
-- 2. 域架构字段：{"domains": [{"name": "HG", "quantity": "large",
--                             "side_chains": false}, ...],
--               "evidence": "<原文结论句>", "source": "text_conclusion"}
-- ----------------------------------------------------------------------------
ALTER TABLE sugars
    ADD COLUMN IF NOT EXISTS domain_architecture JSONB;

CREATE INDEX IF NOT EXISTS idx_sugars_domain_arch
    ON sugars USING gin (domain_architecture);

-- ----------------------------------------------------------------------------
-- 3. 验证
-- ----------------------------------------------------------------------------
-- 3.1 约束已包含新等级:
-- SELECT pg_get_constraintdef(oid) FROM pg_constraint
-- WHERE conname = 'sugars_structure_level_check';
--
-- 3.2 各等级分布:
-- SELECT coalesce(structure_level, '(null)'), count(*) FROM sugars GROUP BY 1 ORDER BY 2 DESC;
--
-- 3.3 域级记录:
-- SELECT iupac_short, domain_architecture FROM sugars WHERE structure_level = 'domain_only';
