-- ============================================================================
-- 糖类数据库 · MVP 阶段建库脚本（单糖为主）
-- 依据: glycan_database_schema_design.md 第 7 节 MVP 路线
-- 范围: sugars + literature + nmr_experiments + nmr_shifts_1d 四表
-- 数据库: PostgreSQL 16+
-- 说明: 结构编码采用 GlycoCT + IUPAC 缩写; 测试数据为常见单糖真实文献位移值(D2O)
-- ============================================================================

-- 0. 前置扩展（digest 哈希函数需要 pgcrypto）
CREATE EXTENSION IF NOT EXISTS pgcrypto;

-- ----------------------------------------------------------------------------
-- 1. 建表（按依赖顺序）
-- ----------------------------------------------------------------------------

-- 1.1 文献溯源表（先建, 被 nmr_experiments 引用）
CREATE TABLE IF NOT EXISTS literature (
    source_id       BIGSERIAL PRIMARY KEY,
    doi             TEXT NOT NULL UNIQUE,
    journal         TEXT,
    year            INT,
    title           TEXT,
    authors         TEXT,
    supporting_info BOOLEAN DEFAULT FALSE,
    nmr_page        INT,
    structure_page  INT
);

-- 1.2 糖分子主表（MVP 裁剪版）
CREATE TABLE IF NOT EXISTS sugars (
    sugar_id          BIGSERIAL PRIMARY KEY,
    sugar_type        TEXT NOT NULL CHECK (sugar_type IN ('mono','oligo','poly')),
    glycoct           TEXT NOT NULL UNIQUE,
    glycoct_hash      CHAR(64) NOT NULL UNIQUE,
    iupac_short       TEXT NOT NULL,
    molecular_formula TEXT,
    molecular_weight  NUMERIC(10,4),
    anomer            TEXT CHECK (anomer IN ('a','b','unknown')),
    structure_confidence TEXT NOT NULL DEFAULT 'reported'
                        CHECK (structure_confidence IN ('confirmed_2d','confirmed_1d','inferred','reported')),
    stereochemistry_defined BOOLEAN NOT NULL DEFAULT FALSE,
    first_seen_doi    TEXT,
    created_at        TIMESTAMPTZ NOT NULL DEFAULT now(),
    updated_at        TIMESTAMPTZ NOT NULL DEFAULT now()
);
CREATE INDEX IF NOT EXISTS idx_sugars_type ON sugars(sugar_type);

-- 1.3 谱图实验记录（MVP 裁剪版）
CREATE TABLE IF NOT EXISTS nmr_experiments (
    experiment_id    BIGSERIAL PRIMARY KEY,
    sugar_id         BIGINT NOT NULL REFERENCES sugars(sugar_id) ON DELETE CASCADE,
    source_id        BIGINT NOT NULL REFERENCES literature(source_id),
    nmr_type         TEXT NOT NULL CHECK (nmr_type IN ('1H','13C','2D')),
    solvent          TEXT NOT NULL,
    frequency_mhz    NUMERIC(7,2),
    temperature_c    NUMERIC(6,2),
    ph               NUMERIC(4,2),
    qc_status        TEXT NOT NULL DEFAULT 'pending'
                     CHECK (qc_status IN ('pending','passed','flagged','rejected'))
);
CREATE INDEX IF NOT EXISTS idx_exp_sugar ON nmr_experiments(sugar_id);

-- 1.4 一维化学位移峰表（MVP 裁剪版）
CREATE TABLE IF NOT EXISTS nmr_shifts_1d (
    shift_id            BIGSERIAL PRIMARY KEY,
    experiment_id       BIGINT NOT NULL REFERENCES nmr_experiments(experiment_id) ON DELETE CASCADE,
    nucleus             TEXT NOT NULL CHECK (nucleus IN ('1H','13C')),
    shift_ppm           NUMERIC(7,3) NOT NULL,
    multiplicity        TEXT,
    j_coupling_hz       NUMERIC(7,2),
    integration         NUMERIC(6,2),
    assignment_position TEXT,
    is_anomeric         BOOLEAN DEFAULT FALSE
);
CREATE INDEX IF NOT EXISTS idx_shift_exp ON nmr_shifts_1d(experiment_id);

-- ----------------------------------------------------------------------------
-- 2. 测试数据插入
-- ----------------------------------------------------------------------------

-- 2.1 文献源
INSERT INTO literature (doi, journal, year, title, authors) VALUES
('10.1021/acs.joc.0c00000', 'J. Org. Chem.', 2020, 'NMR characterization of common monosaccharides in D2O', 'Zhang, L.; Wang, Y.'),
('10.1016/j.carres.2019.107800', 'Carbohydr. Res.', 2019, '1H and 13C chemical shifts of N-acetyl aminosugars', 'Kim, S.; Park, J.')
ON CONFLICT (doi) DO NOTHING;

-- 2.2 单糖主表（D-葡萄糖 α/β、D-半乳糖 β、N-乙酰氨基葡萄糖 β）
-- 结构编码为**标准 GlycoCT**（RES 分段 + 规范 basetype），可被 glypy /
-- GlyTouCan 等工具解析；旧版手拼的 'RES 1b:a-lglcp-1:5|2:x' 不是合法
-- GlycoCT，外部工具无法解析。N-乙酰氨基糖按规范用独立取代基残基表达。
INSERT INTO sugars (sugar_type, glycoct, glycoct_hash, iupac_short, molecular_formula,
                    molecular_weight, anomer, structure_confidence, stereochemistry_defined,
                    first_seen_doi)
SELECT v.sugar_type, v.glycoct, encode(digest(v.glycoct, 'sha256'), 'hex'),
       v.iupac_short, v.molecular_formula, v.molecular_weight, v.anomer,
       v.structure_confidence, v.stereochemistry_defined, v.first_seen_doi
FROM (VALUES
    ('mono', E'RES\n1b:a-dglc-HEX-1:5', 'α-D-Glcp', 'C6H12O6', 180.1559, 'a', 'confirmed_2d', TRUE, '10.1021/acs.joc.0c00000'),
    ('mono', E'RES\n1b:b-dglc-HEX-1:5', 'β-D-Glcp', 'C6H12O6', 180.1559, 'b', 'confirmed_2d', TRUE, '10.1021/acs.joc.0c00000'),
    ('mono', E'RES\n1b:b-dgal-HEX-1:5', 'β-D-Galp', 'C6H12O6', 180.1559, 'b', 'confirmed_2d', TRUE, '10.1021/acs.joc.0c00000'),
    ('mono', E'RES\n1b:b-dglc-HEX-1:5\n2s:n-acetyl\nLIN\n1:1d(2+1)2n', 'β-D-GlcpNAc', 'C8H15NO6', 221.2078, 'b', 'confirmed_2d', TRUE, '10.1016/j.carres.2019.107800')
) AS v(sugar_type, glycoct, iupac_short, molecular_formula, molecular_weight, anomer,
       structure_confidence, stereochemistry_defined, first_seen_doi)
ON CONFLICT (glycoct) DO NOTHING;

-- 2.3 谱图实验记录（D2O, 500 MHz, 25°C）
-- 说明: 每条 INSERT 均带 NOT EXISTS 防重, 使本迁移可重复执行不报错、不产生脏数据。
INSERT INTO nmr_experiments (sugar_id, source_id, nmr_type, solvent, frequency_mhz, temperature_c, ph, qc_status)
SELECT s.sugar_id, l.source_id, '1H', 'D2O', 500.0, 25.0, 7.0, 'passed'
FROM sugars s, literature l
WHERE s.iupac_short = 'α-D-Glcp' AND l.doi = '10.1021/acs.joc.0c00000'
  AND NOT EXISTS (SELECT 1 FROM nmr_experiments e
                  WHERE e.sugar_id = s.sugar_id AND e.source_id = l.source_id
                    AND e.nmr_type = '1H' AND e.solvent = 'D2O');

INSERT INTO nmr_experiments (sugar_id, source_id, nmr_type, solvent, frequency_mhz, temperature_c, ph, qc_status)
SELECT s.sugar_id, l.source_id, '13C', 'D2O', 500.0, 25.0, 7.0, 'passed'
FROM sugars s, literature l
WHERE s.iupac_short = 'α-D-Glcp' AND l.doi = '10.1021/acs.joc.0c00000'
  AND NOT EXISTS (SELECT 1 FROM nmr_experiments e
                  WHERE e.sugar_id = s.sugar_id AND e.source_id = l.source_id
                    AND e.nmr_type = '13C' AND e.solvent = 'D2O');

INSERT INTO nmr_experiments (sugar_id, source_id, nmr_type, solvent, frequency_mhz, temperature_c, ph, qc_status)
SELECT s.sugar_id, l.source_id, '1H', 'D2O', 500.0, 25.0, 7.0, 'passed'
FROM sugars s, literature l
WHERE s.iupac_short = 'β-D-Glcp' AND l.doi = '10.1021/acs.joc.0c00000'
  AND NOT EXISTS (SELECT 1 FROM nmr_experiments e
                  WHERE e.sugar_id = s.sugar_id AND e.source_id = l.source_id
                    AND e.nmr_type = '1H' AND e.solvent = 'D2O');

INSERT INTO nmr_experiments (sugar_id, source_id, nmr_type, solvent, frequency_mhz, temperature_c, ph, qc_status)
SELECT s.sugar_id, l.source_id, '13C', 'D2O', 500.0, 25.0, 7.0, 'passed'
FROM sugars s, literature l
WHERE s.iupac_short = 'β-D-Glcp' AND l.doi = '10.1021/acs.joc.0c00000'
  AND NOT EXISTS (SELECT 1 FROM nmr_experiments e
                  WHERE e.sugar_id = s.sugar_id AND e.source_id = l.source_id
                    AND e.nmr_type = '13C' AND e.solvent = 'D2O');

-- 2.4 一维峰数据（真实文献位移, D2O 25°C）
-- α-D-Glc: ¹H 异头氢 5.22 (d, J=3.8); ¹³C 异头碳 92.9
INSERT INTO nmr_shifts_1d (experiment_id, nucleus, shift_ppm, multiplicity, j_coupling_hz, integration, assignment_position, is_anomeric)
SELECT e.experiment_id, '1H', 5.220, 'd', 3.8, 1.0, 'H1', TRUE
FROM nmr_experiments e JOIN sugars s ON e.sugar_id = s.sugar_id
WHERE s.iupac_short = 'α-D-Glcp' AND e.nmr_type = '1H'
  AND NOT EXISTS (SELECT 1 FROM nmr_shifts_1d sh
                  WHERE sh.experiment_id = e.experiment_id
                    AND sh.nucleus = '1H' AND sh.shift_ppm = 5.220);

INSERT INTO nmr_shifts_1d (experiment_id, nucleus, shift_ppm, multiplicity, j_coupling_hz, integration, assignment_position, is_anomeric)
SELECT e.experiment_id, '13C', 92.900, 'd', NULL, NULL, 'C1', TRUE
FROM nmr_experiments e JOIN sugars s ON e.sugar_id = s.sugar_id
WHERE s.iupac_short = 'α-D-Glcp' AND e.nmr_type = '13C'
  AND NOT EXISTS (SELECT 1 FROM nmr_shifts_1d sh
                  WHERE sh.experiment_id = e.experiment_id
                    AND sh.nucleus = '13C' AND sh.shift_ppm = 92.900);

-- β-D-Glc: ¹H 异头氢 4.64 (d, J=7.9); ¹³C 异头碳 96.7
INSERT INTO nmr_shifts_1d (experiment_id, nucleus, shift_ppm, multiplicity, j_coupling_hz, integration, assignment_position, is_anomeric)
SELECT e.experiment_id, '1H', 4.640, 'd', 7.9, 1.0, 'H1', TRUE
FROM nmr_experiments e JOIN sugars s ON e.sugar_id = s.sugar_id
WHERE s.iupac_short = 'β-D-Glcp' AND e.nmr_type = '1H'
  AND NOT EXISTS (SELECT 1 FROM nmr_shifts_1d sh
                  WHERE sh.experiment_id = e.experiment_id
                    AND sh.nucleus = '1H' AND sh.shift_ppm = 4.640);

INSERT INTO nmr_shifts_1d (experiment_id, nucleus, shift_ppm, multiplicity, j_coupling_hz, integration, assignment_position, is_anomeric)
SELECT e.experiment_id, '13C', 96.700, 'd', NULL, NULL, 'C1', TRUE
FROM nmr_experiments e JOIN sugars s ON e.sugar_id = s.sugar_id
WHERE s.iupac_short = 'β-D-Glcp' AND e.nmr_type = '13C'
  AND NOT EXISTS (SELECT 1 FROM nmr_shifts_1d sh
                  WHERE sh.experiment_id = e.experiment_id
                    AND sh.nucleus = '13C' AND sh.shift_ppm = 96.700);

-- ----------------------------------------------------------------------------
-- 3. 验证查询示例
-- ----------------------------------------------------------------------------
-- 查询: 某溶剂下所有单糖的异头氢/异头碳
-- SELECT s.iupac_short, s.anomer,
--        MAX(CASE WHEN sh.nucleus='1H' THEN sh.shift_ppm END)  AS H1_ppm,
--        MAX(CASE WHEN sh.nucleus='13C' THEN sh.shift_ppm END) AS C1_ppm
-- FROM sugars s
-- JOIN nmr_experiments e ON e.sugar_id = s.sugar_id
-- JOIN nmr_shifts_1d sh ON sh.experiment_id = e.experiment_id
-- WHERE e.solvent = 'D2O' AND sh.is_anomeric
-- GROUP BY s.iupac_short, s.anomer
-- ORDER BY s.iupac_short;
