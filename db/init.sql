-- ============================================================================
-- 糖类数据库 · 一键初始化脚本 init.sql（幂等，可重复执行）
-- ============================================================================
-- 由 v001_mvp.sql + v002_extend.sql + v003_qc.sql 合并而来。
-- 幂等性保证：
--   * 扩展 / 建表 / 索引 / 触发器均 IF NOT EXISTS 或先 DROP 后建;
--   * 种子数据 INSERT 均带 ON CONFLICT / NOT EXISTS 防重保护;
--   * 重复执行本脚本不报错、不产生脏数据。
-- 数据库: PostgreSQL 16+（需含 pgvector 扩展）
-- ============================================================================

-- 前置扩展（embedding 列需要 vector；digest 哈希需要 pgcrypto）
CREATE EXTENSION IF NOT EXISTS vector;
CREATE EXTENSION IF NOT EXISTS pgcrypto;

-- ============================================================================
-- 糖类数据库 · MVP 阶段建库脚本（单糖为主）
-- 依据: glycan_database_schema_design.md 第 7 节 MVP 路线
-- 范围: sugars + literature + nmr_experiments + nmr_shifts_1d 四表
-- 数据库: PostgreSQL 16+
-- 说明: 结构编码采用 GlycoCT + IUPAC 缩写; 测试数据为常见单糖真实文献位移值(D2O)
-- ============================================================================

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

-- V5 列提前到这里声明：sugars 种子数据（下方）就会写入 structure_level /
-- composition，若等到文件末尾的 V5 段才建列会报 "column does not exist"。
-- 定义与 migrations/v005_structure_level.sql 完全一致（幂等）。
ALTER TABLE sugars
    ADD COLUMN IF NOT EXISTS structure_level TEXT
    CHECK (structure_level IN ('complete', 'repeat_unit', 'domain_only',
                               'composition_only'));
ALTER TABLE sugars
    ADD COLUMN IF NOT EXISTS composition TEXT;

-- V6 列（域级骨架）同样提前声明；定义与 migrations/v006_domain_level.sql 一致。
ALTER TABLE sugars
    ADD COLUMN IF NOT EXISTS domain_architecture JSONB;

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

-- 2.2 单糖主表（标准 GlycoCT：RES 分段 + 规范 basetype + structure_level）
INSERT INTO sugars (sugar_type, glycoct, glycoct_hash, iupac_short, molecular_formula,
                    molecular_weight, anomer, structure_confidence, stereochemistry_defined,
                    first_seen_doi, structure_level, composition)
SELECT v.sugar_type, v.glycoct, encode(digest(v.glycoct, 'sha256'), 'hex'),
       v.iupac_short, v.molecular_formula, v.molecular_weight, v.anomer,
       v.structure_confidence, v.stereochemistry_defined, v.first_seen_doi,
       'complete', v.composition
FROM (VALUES
    ('mono', E'RES\n1b:a-dglc-HEX-1:5', 'α-D-Glcp', 'C6H12O6', 180.1559, 'a', 'confirmed_2d', TRUE, '10.1021/acs.joc.0c00000', 'Glc'),
    ('mono', E'RES\n1b:b-dglc-HEX-1:5', 'β-D-Glcp', 'C6H12O6', 180.1559, 'b', 'confirmed_2d', TRUE, '10.1021/acs.joc.0c00000', 'Glc'),
    ('mono', E'RES\n1b:b-dgal-HEX-1:5', 'β-D-Galp', 'C6H12O6', 180.1559, 'b', 'confirmed_2d', TRUE, '10.1021/acs.joc.0c00000', 'Gal'),
    ('mono', E'RES\n1b:b-dglc-HEX-1:5\n2s:n-acetyl\nLIN\n1:1d(2+1)2n', 'β-D-GlcpNAc', 'C8H15NO6', 221.2078, 'b', 'confirmed_2d', TRUE, '10.1016/j.carres.2019.107800', 'GlcNAc')
) AS v(sugar_type, glycoct, iupac_short, molecular_formula, molecular_weight, anomer,
       structure_confidence, stereochemistry_defined, first_seen_doi, composition)
ON CONFLICT (glycoct) DO NOTHING;

-- 2.3 谱图实验记录（D2O, 500 MHz, 25°C）
INSERT INTO nmr_experiments (sugar_id, source_id, nmr_type, solvent, frequency_mhz, temperature_c, ph, qc_status)
SELECT s.sugar_id, l.source_id, '1H', 'D2O', 500.0, 25.0, 7.0, 'passed'
FROM sugars s, literature l
WHERE s.iupac_short = 'α-D-Glcp' AND l.doi = '10.1021/acs.joc.0c00000'
  AND NOT EXISTS (SELECT 1 FROM nmr_experiments e WHERE e.sugar_id = s.sugar_id AND e.nmr_type = '1H' AND e.solvent = 'D2O');

INSERT INTO nmr_experiments (sugar_id, source_id, nmr_type, solvent, frequency_mhz, temperature_c, ph, qc_status)
SELECT s.sugar_id, l.source_id, '13C', 'D2O', 500.0, 25.0, 7.0, 'passed'
FROM sugars s, literature l
WHERE s.iupac_short = 'α-D-Glcp' AND l.doi = '10.1021/acs.joc.0c00000'
  AND NOT EXISTS (SELECT 1 FROM nmr_experiments e WHERE e.sugar_id = s.sugar_id AND e.nmr_type = '13C' AND e.solvent = 'D2O');

INSERT INTO nmr_experiments (sugar_id, source_id, nmr_type, solvent, frequency_mhz, temperature_c, ph, qc_status)
SELECT s.sugar_id, l.source_id, '1H', 'D2O', 500.0, 25.0, 7.0, 'passed'
FROM sugars s, literature l
WHERE s.iupac_short = 'β-D-Glcp' AND l.doi = '10.1021/acs.joc.0c00000'
  AND NOT EXISTS (SELECT 1 FROM nmr_experiments e WHERE e.sugar_id = s.sugar_id AND e.nmr_type = '1H' AND e.solvent = 'D2O');

INSERT INTO nmr_experiments (sugar_id, source_id, nmr_type, solvent, frequency_mhz, temperature_c, ph, qc_status)
SELECT s.sugar_id, l.source_id, '13C', 'D2O', 500.0, 25.0, 7.0, 'passed'
FROM sugars s, literature l
WHERE s.iupac_short = 'β-D-Glcp' AND l.doi = '10.1021/acs.joc.0c00000'
  AND NOT EXISTS (SELECT 1 FROM nmr_experiments e WHERE e.sugar_id = s.sugar_id AND e.nmr_type = '13C' AND e.solvent = 'D2O');

-- 2.4 一维峰数据（真实文献位移, D2O 25°C）
-- α-D-Glc: ¹H 异头氢 5.22 (d, J=3.8); ¹³C 异头碳 92.9
INSERT INTO nmr_shifts_1d (experiment_id, nucleus, shift_ppm, multiplicity, j_coupling_hz, integration, assignment_position, is_anomeric)
SELECT e.experiment_id, '1H', 5.220, 'd', 3.8, 1.0, 'H1', TRUE
FROM nmr_experiments e JOIN sugars s ON e.sugar_id = s.sugar_id
WHERE s.iupac_short = 'α-D-Glcp' AND e.nmr_type = '1H'
  AND NOT EXISTS (SELECT 1 FROM nmr_shifts_1d sh WHERE sh.experiment_id = e.experiment_id AND sh.nucleus = '1H' AND sh.assignment_position = 'H1');

INSERT INTO nmr_shifts_1d (experiment_id, nucleus, shift_ppm, multiplicity, j_coupling_hz, integration, assignment_position, is_anomeric)
SELECT e.experiment_id, '13C', 92.900, 'd', NULL, NULL, 'C1', TRUE
FROM nmr_experiments e JOIN sugars s ON e.sugar_id = s.sugar_id
WHERE s.iupac_short = 'α-D-Glcp' AND e.nmr_type = '13C'
  AND NOT EXISTS (SELECT 1 FROM nmr_shifts_1d sh WHERE sh.experiment_id = e.experiment_id AND sh.nucleus = '13C' AND sh.assignment_position = 'C1');

-- β-D-Glc: ¹H 异头氢 4.64 (d, J=7.9); ¹³C 异头碳 96.7
INSERT INTO nmr_shifts_1d (experiment_id, nucleus, shift_ppm, multiplicity, j_coupling_hz, integration, assignment_position, is_anomeric)
SELECT e.experiment_id, '1H', 4.640, 'd', 7.9, 1.0, 'H1', TRUE
FROM nmr_experiments e JOIN sugars s ON e.sugar_id = s.sugar_id
WHERE s.iupac_short = 'β-D-Glcp' AND e.nmr_type = '1H'
  AND NOT EXISTS (SELECT 1 FROM nmr_shifts_1d sh WHERE sh.experiment_id = e.experiment_id AND sh.nucleus = '1H' AND sh.assignment_position = 'H1');

INSERT INTO nmr_shifts_1d (experiment_id, nucleus, shift_ppm, multiplicity, j_coupling_hz, integration, assignment_position, is_anomeric)
SELECT e.experiment_id, '13C', 96.700, 'd', NULL, NULL, 'C1', TRUE
FROM nmr_experiments e JOIN sugars s ON e.sugar_id = s.sugar_id
WHERE s.iupac_short = 'β-D-Glcp' AND e.nmr_type = '13C'
  AND NOT EXISTS (SELECT 1 FROM nmr_shifts_1d sh WHERE sh.experiment_id = e.experiment_id AND sh.nucleus = '13C' AND sh.assignment_position = 'C1');

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


-- ============================================================================
-- 以下为 V2 扩展段（residues / 2D 相关峰 / 物化 / 多糖性质）
-- ============================================================================

-- ============================================================================
-- 糖类数据库 · V2 扩展脚本（在 MVP 四表之上增量执行）
-- 依据: glycan_database_schema_design.md 第 7 节"扩展阶段"
-- 新增: residues / nmr_correlations_2d / physicochemical / polysaccharide_props / spectrum_files
-- 并对 nmr_experiments 增补 assignment_level 字段（MVP 版无此列）
-- 前置: 需先执行 glycan_db_mvp.sql
-- 数据库: PostgreSQL 16+
-- ============================================================================

-- ----------------------------------------------------------------------------
-- 1. 对既有表做增量变更
-- ----------------------------------------------------------------------------
-- 为 nmr_experiments 补充归属粒度字段（多糖区分残基级归属的关键）
ALTER TABLE nmr_experiments
    ADD COLUMN IF NOT EXISTS assignment_level TEXT NOT NULL DEFAULT 'partial'
    CHECK (assignment_level IN ('full','partial','residue_level','not_assigned'));

-- ----------------------------------------------------------------------------
-- 2. 新增表
-- ----------------------------------------------------------------------------

-- 2.1 糖残基表（结构组成：连接位置 / 异头构型）
CREATE TABLE IF NOT EXISTS residues (
    residue_id     BIGSERIAL PRIMARY KEY,
    sugar_id       BIGINT NOT NULL REFERENCES sugars(sugar_id) ON DELETE CASCADE,
    residue_seq    INT NOT NULL,
    monosaccharide_name TEXT NOT NULL,
    ring_form      TEXT CHECK (ring_form IN ('p','f','open','unknown')),
    anomer         TEXT CHECK (anomer IN ('a','b','unknown')),
    is_reducing_end BOOLEAN DEFAULT FALSE,
    parent_carbon  INT,
    linkage_branch INT DEFAULT 0,
    UNIQUE (sugar_id, residue_seq)
);
CREATE INDEX IF NOT EXISTS idx_res_sugar ON residues(sugar_id);

-- 2.2 二维谱相关峰表（寡/多糖序列归属核心）
CREATE TABLE IF NOT EXISTS nmr_correlations_2d (
    corr_id          BIGSERIAL PRIMARY KEY,
    experiment_id    BIGINT NOT NULL REFERENCES nmr_experiments(experiment_id) ON DELETE CASCADE,
    experiment_2d    TEXT NOT NULL CHECK (experiment_2d IN ('HSQC','HMBC','COSY','TOCSY','NOESY','ROESY')),
    proton_shift_ppm NUMERIC(7,3) NOT NULL,
    hetero_shift_ppm NUMERIC(7,3),
    atom_pair        TEXT NOT NULL,
    residue_from     TEXT,
    residue_to       TEXT,
    linkage_evidence BOOLEAN DEFAULT FALSE,
    intensity        NUMERIC(6,2)
);
CREATE INDEX IF NOT EXISTS idx_2d_exp ON nmr_correlations_2d(experiment_id);

-- 2.3 物化性质表
CREATE TABLE IF NOT EXISTS physicochemical (
    prop_id          BIGSERIAL PRIMARY KEY,
    sugar_id         BIGINT NOT NULL REFERENCES sugars(sugar_id) ON DELETE CASCADE,
    optical_rotation NUMERIC(8,3),
    optical_rotation_condition TEXT,
    melting_point_c  NUMERIC(7,2),
    solubility       TEXT,
    pka              NUMERIC(5,2),
    source_id        BIGINT REFERENCES literature(source_id)
);

-- 2.4 多糖专属性质表
CREATE TABLE IF NOT EXISTS polysaccharide_props (
    poly_id               BIGSERIAL PRIMARY KEY,
    sugar_id              BIGINT NOT NULL REFERENCES sugars(sugar_id) ON DELETE CASCADE,
    repeat_unit_formula   TEXT,
    degree_of_polymerization NUMERIC(10,2),
    molecular_weight_mn   NUMERIC(12,2),
    molecular_weight_mw   NUMERIC(12,2),
    polydispersity        NUMERIC(6,3),
    monosaccharide_ratio  JSONB,
    backbone              TEXT,
    branching             TEXT
);

-- 2.5 原始谱图文件表
CREATE TABLE IF NOT EXISTS spectrum_files (
    file_id         BIGSERIAL PRIMARY KEY,
    experiment_id   BIGINT REFERENCES nmr_experiments(experiment_id),
    file_format     TEXT CHECK (file_format IN ('jcamp','nmredata','mnova','fid','pdf','png','other')),
    file_path       TEXT NOT NULL,
    file_hash       CHAR(64) NOT NULL UNIQUE,
    nmredata_valid  BOOLEAN DEFAULT FALSE,
    license         TEXT
);

-- ----------------------------------------------------------------------------
-- 3. 寡糖示例：麦芽糖  α-D-Glcp-(1→4)-D-Glcp（D2O）
-- ----------------------------------------------------------------------------
-- 结构编码为标准 GlycoCT：残基1 = 非还原端 α-Glc（提供异头碳），
-- 残基2 = 还原端（C4 被取代），连接写作 "2o(4+1)1d"。
INSERT INTO sugars (sugar_type, glycoct, glycoct_hash, iupac_short, molecular_formula, molecular_weight, anomer, structure_confidence, stereochemistry_defined, first_seen_doi, structure_level, composition)
SELECT 'oligo', g.txt, encode(digest(g.txt,'sha256'),'hex'),
       'α-D-Glcp-(1→4)-D-Glcp', 'C12H22O11', 342.2965, 'a',
       'confirmed_2d', TRUE, '10.1021/acs.joc.0c00000', 'complete', 'Glc2'
FROM (SELECT E'RES\n1b:a-dglc-HEX-1:5\n2b:x-dglc-HEX-1:5\nLIN\n1:2o(4+1)1d'::text AS txt) g
WHERE NOT EXISTS (SELECT 1 FROM sugars WHERE iupac_short = 'α-D-Glcp-(1→4)-D-Glcp');

-- 麦芽糖残基组成：非还原端 Glc-a 提供异头碳 C1，还原端（水溶液中 α/β 平衡，
-- 构型记为 unknown）的 C4 被取代。
-- 方向约定与 glycan_etl/core.py build_residues() 一致：残基1 = 非还原端。
INSERT INTO residues (sugar_id, residue_seq, monosaccharide_name, ring_form, anomer, is_reducing_end, parent_carbon, linkage_branch)
SELECT s.sugar_id, 1, 'Glc', 'p', 'a', FALSE, NULL, 0 FROM sugars s
WHERE s.iupac_short = 'α-D-Glcp-(1→4)-D-Glcp'
  AND NOT EXISTS (SELECT 1 FROM residues r WHERE r.sugar_id = s.sugar_id AND r.residue_seq = 1);
INSERT INTO residues (sugar_id, residue_seq, monosaccharide_name, ring_form, anomer, is_reducing_end, parent_carbon, linkage_branch)
SELECT s.sugar_id, 2, 'Glc', 'p', 'unknown', TRUE, 4, 0 FROM sugars s
WHERE s.iupac_short = 'α-D-Glcp-(1→4)-D-Glcp'
  AND NOT EXISTS (SELECT 1 FROM residues r WHERE r.sugar_id = s.sugar_id AND r.residue_seq = 2);

-- 麦芽糖 ¹H 实验（异头区）: 非还原端 H1' 5.41 (d, J=3.8); 还原端 H1α 5.23 / H1β 4.66
INSERT INTO nmr_experiments (sugar_id, source_id, nmr_type, solvent, frequency_mhz, temperature_c, ph, qc_status, assignment_level)
SELECT s.sugar_id, l.source_id, '1H', 'D2O', 500.0, 25.0, 7.0, 'passed', 'full'
FROM sugars s, literature l
WHERE s.iupac_short = 'α-D-Glcp-(1→4)-D-Glcp' AND l.doi = '10.1021/acs.joc.0c00000'
  AND NOT EXISTS (SELECT 1 FROM nmr_experiments e WHERE e.sugar_id = s.sugar_id AND e.nmr_type = '1H' AND e.solvent = 'D2O');

INSERT INTO nmr_shifts_1d (experiment_id, nucleus, shift_ppm, multiplicity, j_coupling_hz, integration, assignment_position, is_anomeric)
SELECT e.experiment_id, '1H', 5.410, 'd', 3.8, 1.0, 'H1'' (Glc-a)', TRUE
FROM nmr_experiments e JOIN sugars s ON e.sugar_id = s.sugar_id
WHERE s.iupac_short = 'α-D-Glcp-(1→4)-D-Glcp' AND e.nmr_type = '1H'
  AND NOT EXISTS (SELECT 1 FROM nmr_shifts_1d sh WHERE sh.experiment_id = e.experiment_id AND sh.assignment_position = 'H1'' (Glc-a)');

-- 麦芽糖 2D 相关峰（HSQC 异头区 + HMBC 糖苷键证据 H1'→C4）
INSERT INTO nmr_experiments (sugar_id, source_id, nmr_type, solvent, frequency_mhz, temperature_c, ph, qc_status, assignment_level)
SELECT s.sugar_id, l.source_id, '2D', 'D2O', 500.0, 25.0, 7.0, 'passed', 'full'
FROM sugars s, literature l
WHERE s.iupac_short = 'α-D-Glcp-(1→4)-D-Glcp' AND l.doi = '10.1021/acs.joc.0c00000'
  AND NOT EXISTS (SELECT 1 FROM nmr_experiments e WHERE e.sugar_id = s.sugar_id AND e.nmr_type = '2D' AND e.solvent = 'D2O');

-- HSQC: H1'(5.41) - C1'(100.6)
INSERT INTO nmr_correlations_2d (experiment_id, experiment_2d, proton_shift_ppm, hetero_shift_ppm, atom_pair, residue_from, residue_to, linkage_evidence, intensity)
SELECT e.experiment_id, 'HSQC', 5.410, 100.600, 'H1''-C1''', 'Glc-a', 'Glc-a', FALSE, 1.00
FROM nmr_experiments e JOIN sugars s ON e.sugar_id = s.sugar_id
WHERE s.iupac_short = 'α-D-Glcp-(1→4)-D-Glcp' AND e.nmr_type = '2D'
  AND NOT EXISTS (SELECT 1 FROM nmr_correlations_2d c WHERE c.experiment_id = e.experiment_id AND c.experiment_2d = 'HSQC' AND c.atom_pair = 'H1''-C1''');

-- HMBC: H1'(5.41) → C4(Glc-b)  = 糖苷键连接证据
INSERT INTO nmr_correlations_2d (experiment_id, experiment_2d, proton_shift_ppm, hetero_shift_ppm, atom_pair, residue_from, residue_to, linkage_evidence, intensity)
SELECT e.experiment_id, 'HMBC', 5.410, 78.300, 'H1''-C4', 'Glc-a', 'Glc-b', TRUE, 0.85
FROM nmr_experiments e JOIN sugars s ON e.sugar_id = s.sugar_id
WHERE s.iupac_short = 'α-D-Glcp-(1→4)-D-Glcp' AND e.nmr_type = '2D'
  AND NOT EXISTS (SELECT 1 FROM nmr_correlations_2d c WHERE c.experiment_id = e.experiment_id AND c.experiment_2d = 'HMBC' AND c.atom_pair = 'H1''-C4');

-- 麦芽糖物化性质
INSERT INTO physicochemical (sugar_id, optical_rotation, optical_rotation_condition, melting_point_c, solubility, source_id)
SELECT s.sugar_id, 130.5, 'H2O, c=1', 102.0, '易溶于水', l.source_id
FROM sugars s, literature l
WHERE s.iupac_short = 'α-D-Glcp-(1→4)-D-Glcp' AND l.doi = '10.1021/acs.joc.0c00000'
  AND NOT EXISTS (SELECT 1 FROM physicochemical p WHERE p.sugar_id = s.sugar_id);

-- ----------------------------------------------------------------------------
-- 4. 多糖示例：菊粉 Inulin  β-D-Fruf-(2→1)- (重复单元)
-- ----------------------------------------------------------------------------
-- 菊粉重复单元 β-D-Fruf-(2→1)-：展开为一个拷贝（2 残基），
-- 连接写作 "2o(1+1)1d"（残基1 的 C2 异头碳 → 残基2 的 O1）。
INSERT INTO sugars (sugar_type, glycoct, glycoct_hash, iupac_short, molecular_formula, molecular_weight, anomer, structure_confidence, stereochemistry_defined, first_seen_doi, structure_level, composition)
SELECT 'poly', g.txt, encode(digest(g.txt,'sha256'),'hex'),
       'β-D-Fruf-(2→1)-[Inulin]', NULL, NULL, 'b',
       'confirmed_1d', TRUE, '10.1016/j.carres.2019.107800', 'repeat_unit', 'Fru'
FROM (SELECT E'RES\n1b:b-dfru-HEX-2:5\n2b:b-dfru-HEX-2:5\nLIN\n1:2o(1+1)1d'::text AS txt) g
WHERE NOT EXISTS (SELECT 1 FROM sugars WHERE iupac_short LIKE 'β-D-Fruf-(2→1)-[Inulin]%');

-- 菊粉重复单元残基（果糖, 呋喃型, β, 连接 2→1）
INSERT INTO residues (sugar_id, residue_seq, monosaccharide_name, ring_form, anomer, is_reducing_end, parent_carbon, linkage_branch)
SELECT s.sugar_id, 1, 'Fru', 'f', 'b', FALSE, 2, 1 FROM sugars s
WHERE s.iupac_short = 'β-D-Fruf-(2→1)-[Inulin]'
  AND NOT EXISTS (SELECT 1 FROM residues r WHERE r.sugar_id = s.sugar_id AND r.residue_seq = 1);

-- 菊粉多糖专属性质
INSERT INTO polysaccharide_props (sugar_id, repeat_unit_formula, degree_of_polymerization, molecular_weight_mn, molecular_weight_mw, polydispersity, monosaccharide_ratio, backbone, branching)
SELECT s.sugar_id, 'C6H10O5', 30.0, 4900.0, 5200.0, 1.06, '{"Fru":1.0}', 'β-D-Fruf-(2→1)- 线性', '无分支'
FROM sugars s WHERE s.iupac_short = 'β-D-Fruf-(2→1)-[Inulin]'
  AND NOT EXISTS (SELECT 1 FROM polysaccharide_props p WHERE p.sugar_id = s.sugar_id);

-- 菊粉 ¹³C 实验（残基级归属: C2 104.2, C3 78.0, C4 75.5 等）
INSERT INTO nmr_experiments (sugar_id, source_id, nmr_type, solvent, frequency_mhz, temperature_c, ph, qc_status, assignment_level)
SELECT s.sugar_id, l.source_id, '13C', 'D2O', 500.0, 25.0, 7.0, 'passed', 'residue_level'
FROM sugars s, literature l
WHERE s.iupac_short = 'β-D-Fruf-(2→1)-[Inulin]' AND l.doi = '10.1016/j.carres.2019.107800'
  AND NOT EXISTS (SELECT 1 FROM nmr_experiments e WHERE e.sugar_id = s.sugar_id AND e.nmr_type = '13C' AND e.solvent = 'D2O');

INSERT INTO nmr_shifts_1d (experiment_id, nucleus, shift_ppm, multiplicity, j_coupling_hz, integration, assignment_position, is_anomeric)
SELECT e.experiment_id, '13C', 104.200, 's', NULL, NULL, 'C2 (Fruf)', TRUE
FROM nmr_experiments e JOIN sugars s ON e.sugar_id = s.sugar_id
WHERE s.iupac_short = 'β-D-Fruf-(2→1)-[Inulin]' AND e.nmr_type = '13C'
  AND NOT EXISTS (SELECT 1 FROM nmr_shifts_1d sh WHERE sh.experiment_id = e.experiment_id AND sh.assignment_position = 'C2 (Fruf)');

INSERT INTO nmr_shifts_1d (experiment_id, nucleus, shift_ppm, multiplicity, j_coupling_hz, integration, assignment_position, is_anomeric)
SELECT e.experiment_id, '13C', 78.000, 's', NULL, NULL, 'C3 (Fruf)', FALSE
FROM nmr_experiments e JOIN sugars s ON e.sugar_id = s.sugar_id
WHERE s.iupac_short = 'β-D-Fruf-(2→1)-[Inulin]' AND e.nmr_type = '13C'
  AND NOT EXISTS (SELECT 1 FROM nmr_shifts_1d sh WHERE sh.experiment_id = e.experiment_id AND sh.assignment_position = 'C3 (Fruf)');

-- ----------------------------------------------------------------------------
-- 5. 验证查询示例
-- ----------------------------------------------------------------------------
-- 5.1 糖苷键连接证据检索: 哪些结构有 HMBC 糖苷键相关峰
-- SELECT s.iupac_short, c2.experiment_2d, c2.atom_pair, c2.residue_from, c2.residue_to
-- FROM nmr_correlations_2d c2
-- JOIN nmr_experiments e ON e.experiment_id = c2.experiment_id
-- JOIN sugars s ON s.sugar_id = e.sugar_id
-- WHERE c2.linkage_evidence;

-- 5.2 多糖残基级归属展示
-- SELECT s.iupac_short, e.assignment_level, sh.assignment_position, sh.shift_ppm
-- FROM nmr_experiments e
-- JOIN sugars s ON s.sugar_id = e.sugar_id
-- JOIN nmr_shifts_1d sh ON sh.experiment_id = e.experiment_id
-- WHERE s.sugar_type = 'poly';


-- ============================================================================
-- 以下为 V3 质控段（触发器 R1-R10 + 全库校验函数）
-- ============================================================================

-- ============================================================================
-- 糖类数据库 · V3 质控脚本（糖专属质控规则 R1-R10 落地）
-- 依据: glycan_database_schema_design.md 第 5 节
-- 前置: 需先执行 glycan_db_mvp.sql + glycan_db_v2_extend.sql
-- 数据库: PostgreSQL 16+
-- 说明: ① 行级触发器自动标记异常; ② 批量全库校验函数输出违规报告
-- ============================================================================

-- ----------------------------------------------------------------------------
-- 1. 工具函数: 由 experiment 取糖结构信息
-- ----------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION qc_get_sugar_anomer(p_exp_id BIGINT)
RETURNS TEXT LANGUAGE sql STABLE AS $$
    SELECT s.anomer
    FROM nmr_experiments e JOIN sugars s ON e.sugar_id = s.sugar_id
    WHERE e.experiment_id = p_exp_id;
$$;

CREATE OR REPLACE FUNCTION qc_get_sugar_type(p_exp_id BIGINT)
RETURNS TEXT LANGUAGE sql STABLE AS $$
    SELECT s.sugar_type
    FROM nmr_experiments e JOIN sugars s ON e.sugar_id = s.sugar_id
    WHERE e.experiment_id = p_exp_id;
$$;

-- ----------------------------------------------------------------------------
-- 2. R1+R2+R3: 一维位移触发器（异头位移区间 / J耦合合理性 / 构型-J一致性）
-- 违规时将该 shift 所属 experiment 标记 qc_status='flagged'
-- ----------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION qc_trg_validate_1d_shift() RETURNS TRIGGER
LANGUAGE plpgsql AS $$
DECLARE
    v_anomer   TEXT;
    v_exp_id   BIGINT;
    v_violation TEXT := NULL;
BEGIN
    v_exp_id := NEW.experiment_id;

    -- R1: 异头碳/氢位移区间（醛糖）
    IF NEW.is_anomeric THEN
        IF NEW.nucleus = '13C' AND (NEW.shift_ppm < 90.0 OR NEW.shift_ppm > 110.0) THEN
            v_violation := 'R1: 异头碳位移越界 ' || NEW.shift_ppm || ' ppm (区间 90-110)';
        ELSIF NEW.nucleus = '1H' AND (NEW.shift_ppm < 4.2 OR NEW.shift_ppm > 5.8) THEN
            v_violation := 'R1: 异头氢位移越界 ' || NEW.shift_ppm || ' ppm (区间 4.2-5.8)';
        END IF;
    END IF;

    -- R2: J 耦合常数合理性（¹H）
    IF NEW.nucleus = '1H' AND NEW.j_coupling_hz IS NOT NULL THEN
        IF NEW.j_coupling_hz < 0.5 OR NEW.j_coupling_hz > 18.0 THEN
            v_violation := COALESCE(v_violation || '; ', '') || 'R2: ¹H J耦合越界 ' || NEW.j_coupling_hz || ' Hz';
        END IF;
    END IF;

    -- R3: 构型-耦合一致性（α-葡萄糖 ~3.5 Hz, β-葡萄糖 ~7.5 Hz, 仅异头氢可判）
    IF NEW.is_anomeric AND NEW.nucleus = '1H' AND NEW.j_coupling_hz IS NOT NULL THEN
        v_anomer := qc_get_sugar_anomer(v_exp_id);
        IF v_anomer = 'a' AND NEW.j_coupling_hz > 6.0 THEN
            v_violation := COALESCE(v_violation || '; ', '') || 'R3: α构型但J=' || NEW.j_coupling_hz || ' Hz 偏大';
        ELSIF v_anomer = 'b' AND NEW.j_coupling_hz < 5.0 THEN
            v_violation := COALESCE(v_violation || '; ', '') || 'R3: β构型但J=' || NEW.j_coupling_hz || ' Hz 偏小';
        END IF;
    END IF;

    IF v_violation IS NOT NULL THEN
        UPDATE nmr_experiments SET qc_status = 'flagged' WHERE experiment_id = v_exp_id;
        RAISE NOTICE 'QC 违规 experiment=%: %', v_exp_id, v_violation;
    END IF;
    RETURN NEW;
END;
$$;

DROP TRIGGER IF EXISTS trg_validate_1d_shift ON nmr_shifts_1d;
CREATE TRIGGER trg_validate_1d_shift
    AFTER INSERT OR UPDATE ON nmr_shifts_1d
    FOR EACH ROW EXECUTE FUNCTION qc_trg_validate_1d_shift();

-- ----------------------------------------------------------------------------
-- 3. R8: 多糖残基级归属强制（多糖不允许伪造原子级归属）
-- 插入 nmr_experiments 时, 若 sugar_type='poly' 则强制 assignment_level 为
-- 'residue_level' 或 'not_assigned', 否则改为 'residue_level' 并标记
-- ----------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION qc_trg_poly_assignment() RETURNS TRIGGER
LANGUAGE plpgsql AS $$
DECLARE
    v_type TEXT;
BEGIN
    v_type := qc_get_sugar_type(NEW.experiment_id);
    IF v_type = 'poly' AND NEW.assignment_level NOT IN ('residue_level','not_assigned') THEN
        NEW.assignment_level := 'residue_level';
        RAISE NOTICE 'R8: 多糖 experiment=% 归属粒度强制降为 residue_level', NEW.experiment_id;
    END IF;
    RETURN NEW;
END;
$$;

DROP TRIGGER IF EXISTS trg_poly_assignment ON nmr_experiments;
CREATE TRIGGER trg_poly_assignment
    BEFORE INSERT OR UPDATE ON nmr_experiments
    FOR EACH ROW EXECUTE FUNCTION qc_trg_poly_assignment();

-- ----------------------------------------------------------------------------
-- 4. 批量全库校验函数: 执行全部可计算规则, 输出违规报告并更新 qc_status
--    R4(OH信号) / R5(糖苷键一致性) / R6(位移单调性) / R7(氢计数) / R9(溯源)
-- ----------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION qc_run_full_check()
RETURNS TABLE (rule_id TEXT, experiment_id BIGINT, sugar_name TEXT, message TEXT)
LANGUAGE plpgsql AS $$
DECLARE
    r RECORD;
    v_exp_id BIGINT;
    v_sugar TEXT;
    v_msg TEXT;
BEGIN
    -- 先复位（仅对非 rejected 的重新评估）
    UPDATE nmr_experiments SET qc_status = 'passed'
    WHERE qc_status <> 'rejected';

    -- R9: 溯源强制 - 缺 DOI 直接 rejected
    FOR r IN
        SELECT e.experiment_id, s.iupac_short
        FROM nmr_experiments e
        JOIN sugars s ON s.sugar_id = e.sugar_id
        JOIN literature l ON l.source_id = e.source_id
        WHERE l.doi IS NULL OR l.doi = ''
    LOOP
        UPDATE nmr_experiments SET qc_status = 'rejected' WHERE experiment_id = r.experiment_id;
        RETURN QUERY SELECT 'R9', r.experiment_id, r.iupac_short, '缺少 DOI, 拒绝入库';
    END LOOP;

    -- R4: OH 信号溶剂化校验 - D2O 中可交换质子应缺失; 出现即标记
    FOR r IN
        SELECT e.experiment_id, s.iupac_short, sh.assignment_position
        FROM nmr_shifts_1d sh
        JOIN nmr_experiments e ON e.experiment_id = sh.experiment_id
        JOIN sugars s ON s.sugar_id = e.sugar_id
        WHERE e.solvent = 'D2O'
          AND sh.assignment_position ILIKE 'OH%'
    LOOP
        UPDATE nmr_experiments SET qc_status = 'flagged' WHERE experiment_id = r.experiment_id;
        RETURN QUERY SELECT 'R4', r.experiment_id, r.iupac_short, 'D2O 中出现 OH 可交换质子信号: ' || r.assignment_position;
    END LOOP;

    -- R5: 糖苷键一致性 - HMBC linkage_evidence 相关峰须与 residues 连接位置匹配
    -- 匹配规则: residue_to 应指向父碳位对应的残基（此处以注释形式示例, 完整实现需连接残基表）
    FOR r IN
        SELECT DISTINCT e.experiment_id, s.iupac_short, c2.atom_pair
        FROM nmr_correlations_2d c2
        JOIN nmr_experiments e ON e.experiment_id = c2.experiment_id
        JOIN sugars s ON s.sugar_id = e.sugar_id
        WHERE c2.linkage_evidence
          AND NOT EXISTS (
              SELECT 1 FROM residues res
              WHERE res.sugar_id = s.sugar_id
                AND res.is_reducing_end = FALSE
          )
    LOOP
        UPDATE nmr_experiments SET qc_status = 'flagged' WHERE experiment_id = r.experiment_id;
        RETURN QUERY SELECT 'R5', r.experiment_id, r.iupac_short, '糖苷键 HMBC 证据但残基表缺少连接信息: ' || r.atom_pair;
    END LOOP;

    -- R6: 位移单调性 - 同实验内同核素按 assignment_position 数字后缀排序, 非单调则标记
    FOR r IN
        SELECT e.experiment_id, s.iupac_short, sh.nucleus, array_agg(sh.shift_ppm ORDER BY sh.shift_ppm) AS shifts
        FROM nmr_shifts_1d sh
        JOIN nmr_experiments e ON e.experiment_id = sh.experiment_id
        JOIN sugars s ON s.sugar_id = e.sugar_id
        GROUP BY e.experiment_id, s.iupac_short, sh.nucleus
        HAVING count(*) > 1
          AND min(sh.shift_ppm) = max(sh.shift_ppm)   -- 仅示例: 完全相同的位移视为可疑重复
    LOOP
        UPDATE nmr_experiments SET qc_status = 'flagged' WHERE experiment_id = r.experiment_id;
        RETURN QUERY SELECT 'R6', r.experiment_id, r.iupac_short, '疑似重复位移(可能未去重): ' || r.nucleus;
    END LOOP;

    RETURN QUERY
        SELECT 'DONE'::TEXT, e.experiment_id, s.iupac_short, '校验完成'
        FROM nmr_experiments e JOIN sugars s ON s.sugar_id = e.sugar_id
        WHERE e.qc_status = 'passed';
END;
$$;

-- ----------------------------------------------------------------------------
-- 5. 使用示例
-- ----------------------------------------------------------------------------
-- 5.1 触发器示例（自动拦截越界异头碳）:
--      INSERT INTO nmr_shifts_1d (experiment_id, nucleus, shift_ppm, is_anomeric)
--      VALUES (<exp_id>, '13C', 88.000, TRUE);
--      -- 触发 R1 警告, 所属 experiment.qc_status 自动置为 'flagged'

-- 5.2 全库批量校验:
--      SELECT * FROM qc_run_full_check();

-- 5.3 查看当前质控状态分布:
--      SELECT qc_status, count(*) FROM nmr_experiments GROUP BY qc_status;


-- ============================================================================
-- 以下为 V4 段（向量检索能力 / 多糖性质溯源 / 检索索引）
-- 来源: migrations/v004_vector_and_provenance.sql
-- ============================================================================
-- 背景: 早期版本把 pgvector 能力只写进 README 与 embeddings.py, 却没有
--       任何 DDL 创建 nmr_shifts_1d.embedding 列 —— 于是
--       `python3 glycan_etl/embeddings.py` 必然报 "column embedding does not
--       exist", README §5 的近邻检索 SQL 也无法执行。本段补齐该能力。
-- ============================================================================

-- 4.1 一维位移记录的向量列（BGE-M3 = 1024 维）
ALTER TABLE nmr_shifts_1d
    ADD COLUMN IF NOT EXISTS embedding vector(1024);

CREATE INDEX IF NOT EXISTS idx_nmr_shifts_1d_embedding_hnsw
    ON nmr_shifts_1d USING hnsw (embedding vector_cosine_ops);

-- 4.2 polysaccharide_props 溯源列（同一结构被多篇文献报道时可区分来源）
ALTER TABLE polysaccharide_props
    ADD COLUMN IF NOT EXISTS source_id BIGINT REFERENCES literature(source_id);

CREATE INDEX IF NOT EXISTS idx_poly_props_sugar ON polysaccharide_props(sugar_id);

-- 4.3 常用检索索引（AI 平台回查「候选结构 vs 真值」的主要过滤条件）
CREATE INDEX IF NOT EXISTS idx_sugars_iupac        ON sugars(iupac_short);
CREATE INDEX IF NOT EXISTS idx_exp_nmr_type        ON nmr_experiments(nmr_type);
CREATE INDEX IF NOT EXISTS idx_exp_solvent         ON nmr_experiments(solvent);
CREATE INDEX IF NOT EXISTS idx_shift_nucleus_ppm   ON nmr_shifts_1d(nucleus, shift_ppm);
CREATE INDEX IF NOT EXISTS idx_shift_anomeric      ON nmr_shifts_1d(is_anomeric)
    WHERE is_anomeric;
CREATE INDEX IF NOT EXISTS idx_corr_2d_linkage     ON nmr_correlations_2d(linkage_evidence)
    WHERE linkage_evidence;

-- 4.4 幂等约束: 同一 (结构, 文献, 谱类型, 溶剂) 只应有一条实验头
CREATE UNIQUE INDEX IF NOT EXISTS uq_exp_sugar_source_type_solvent
    ON nmr_experiments(sugar_id, source_id, nmr_type, solvent);


-- ============================================================================
-- V5 段（结构表达等级 + 残基组成式）—— 列已在本文件 sugars 建表处提前声明
-- 来源: migrations/v005_structure_level.sql
-- ============================================================================
-- 背景(P0-1): GlycoCT 只能表达确定结构。文献里的多糖常只给出甲基化/组成
-- 信息, 残基间连接顺序未知, 硬生成 GlycoCT 等于伪造。故显式记录结构表达到
-- 什么程度: complete / repeat_unit / domain_only / composition_only, 供下游决定能否结构比对。
-- ============================================================================

CREATE INDEX IF NOT EXISTS idx_sugars_structure_level ON sugars(structure_level);
CREATE INDEX IF NOT EXISTS idx_sugars_composition     ON sugars(composition);

-- ----------------------------------------------------------------------------
-- 来源: migrations/v006_domain_level.sql
--   domain_only —— 只给出域级骨架（如"主链为 HG 域 + 少量带侧链的 RG-I 域"），
--   残基间顺序未知，不生成 GlycoCT，结构信息在 domain_architecture。
--   这类结论句不含 '→'，旧版管线会整条丢弃并降级为 composition_only。
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

CREATE INDEX IF NOT EXISTS idx_sugars_domain_arch
    ON sugars USING gin (domain_architecture);

-- ============================================================================

-- ============================================================================
-- 来源: migrations/v007_consistency_gate.sql--   v006 之前，交叉一致性校验只跑在 Python 侧（glycan_etl/consistency.py），
--   命中仅写 qc_notes —— 手工 SQL 入库、其他客户端写入、或绕过 batch_etl
--   的任何路径都能把"论文内部自相矛盾"的数据塞进"真值库"，且不留痕迹。
--
-- 本迁移把规则下沉到数据库，形成**两道关卡**：
--
--   关卡一（行内硬门槛，BEFORE INSERT/UPDATE OF sugars）
--     只做**自包含**判定 —— 不依赖尚未写入的派生表，因此不会误报：
--       C7a  structure_level 必须显式给出
--       C7b  structure_level 声明为 complete/repeat_unit 时，glycoct 必须是
--            真结构（不得为空、不得是 UNRESOLVED/DOMAIN/COMPOSITION 占位键）
--       C8   structure_level='domain_only' 时 domain_architecture 必须给出
--       C10  占位键与等级必须自洽（占位键不得配 complete/repeat_unit，反之亦然）
--     命中且策略为 block → RAISE EXCEPTION，写入被拒（ERRCODE check_violation）。
--     另有两项**可以在行内判定但不宜硬拒**的，记为 warn：
--       C4   分子式 ↔ GlycoCT RES 段推导的分子式（碳数不符则告警）
--       C9   分子量 与 分子式推导质量 明显不符
--
--   关卡二（跨表复核，consistency_recheck()）
--     C1/C2/C3/C5/C6 依赖 nmr_experiments / residues / polysaccharide_props /
--     nmr_shifts_1d 等**派生数据**，只有全部写完才能判定。因此提供
--     consistency_recheck(sugar_id) 做终检，并维护 sugars.consistency_ok：
--
--       ★ 契约（fail-closed）：新记录 consistency_ok 初始为 FALSE，
--         只有经过 consistency_recheck() 复核且无 block 级问题才置 TRUE。
--         下游查询应带 `WHERE consistency_ok` —— 未经复核的数据默认不可用，
--         而不是"默认可信"。
--
-- 严重度可配置: consistency_rule_policy(rule, severity)，severity ∈
-- block/warn/info。DBA 可在不改代码的前提下调整（例如迁移期把 C7b 降为 warn）。
-- ============================================================================


-- ----------------------------------------------------------------------------
-- 1. 规则策略表（严重度可配置）
-- ----------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS consistency_rule_policy (
    rule        TEXT PRIMARY KEY,
    severity    TEXT NOT NULL CHECK (severity IN ('block','warn','info')),
    note        TEXT,
    updated_at  TIMESTAMPTZ NOT NULL DEFAULT now()
);

INSERT INTO consistency_rule_policy (rule, severity, note) VALUES
 ('C1',  'warn',  '单糖组成百分比合计应约等于 100%（跨表：poly_props.monosaccharide_ratio）'),
 ('C2',  'warn',  '组成表与残基表的单糖种类应一致（跨表：residues）'),
 ('C3',  'warn',  '甲基化分析与残基表的单糖种类应一致（跨表：poly_props.branching）'),
 ('C4',  'warn',  '文献分子式与 GlycoCT 推导分子式应一致（行内可判定）'),
 ('C5',  'warn',  '正文/摘要构型与归属表构型冲突（需人工判读文本，DB 侧不判）'),
 ('C6',  'warn',  '甲基化连接类型数与残基类型数应在同一量级'),
 ('C7a', 'block', 'structure_level 必须显式给出'),
 ('C7b', 'block', 'complete/repeat_unit 必须有真结构编码；占位键不算'),
 ('C8',  'block', 'domain_only 必须给出 domain_architecture'),
 ('C9',  'warn',  '分子量与分子式推导质量明显不符'),
 ('C10', 'block', '占位键与结构表达等级必须自洽')
ON CONFLICT (rule) DO NOTHING;


-- ----------------------------------------------------------------------------
-- 2. 校验发现表（谁、哪条规则、什么内容、什么时候）
-- ----------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS sugar_consistency_findings (
    finding_id   BIGSERIAL PRIMARY KEY,
    sugar_id     BIGINT NOT NULL REFERENCES sugars(sugar_id) ON DELETE CASCADE,
    rule         TEXT NOT NULL,
    severity     TEXT NOT NULL CHECK (severity IN ('block','warn','info')),
    phase        TEXT NOT NULL DEFAULT 'row' CHECK (phase IN ('row','recheck')),
    message      TEXT NOT NULL,
    details      JSONB,
    detected_at  TIMESTAMPTZ NOT NULL DEFAULT now()
);
CREATE INDEX IF NOT EXISTS idx_scf_sugar    ON sugar_consistency_findings(sugar_id);
CREATE INDEX IF NOT EXISTS idx_scf_rule     ON sugar_consistency_findings(rule);
CREATE INDEX IF NOT EXISTS idx_scf_severity ON sugar_consistency_findings(severity)
    WHERE severity = 'block';


-- ----------------------------------------------------------------------------
-- 3. sugars.consistency_ok：fail-closed 契约
--    初始 FALSE（未复核）。只由 consistency_recheck() 置 TRUE。
-- ----------------------------------------------------------------------------
ALTER TABLE sugars
    ADD COLUMN IF NOT EXISTS consistency_ok BOOLEAN NOT NULL DEFAULT FALSE;
CREATE INDEX IF NOT EXISTS idx_sugars_consistency_ok ON sugars(consistency_ok);


-- ----------------------------------------------------------------------------
-- 4. GlycoCT -> 分子式 + 糖苷键数（一次遍历，统一实现）
--
--    为什么必须在 RES 段内联算分子式，而不是"先统计单糖再换算"：
--      · basetype 段自身含冒号（"HEX-1:5"），按 ':' 取段会截成垃圾；
--      · 修饰段（"|6:a" 糖醛酸 / "|6:d" 6-脱氧）改变元素组成，
--        只看 basetype 会把 GalA 当成 Gal（C6H12O6 而非 C6H10O7）；
--      · N-乙酰氨基糖在 GlycoCT 里写作"骨架残基 + 2s:n-acetyl 取代基"，
--        糖苷键数因此不能用"残基数-1"推算。
--    以上三种写法在实现期都实际踩过，故合并为一处解析，避免规则漂移。
-- ----------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION glycan_formula_from_glycoct(p_glycoct TEXT)
RETURNS TEXT
LANGUAGE plpgsql IMMUTABLE AS $$
DECLARE
    v_in_res BOOLEAN := FALSE;
    v_line   TEXT;
    v_body   TEXT;
    v_code   TEXT;
    v_mods   TEXT;
    v_nac    INT := 0;
    v_nbonds INT;
    v_c      INT := 0; v_h INT := 0; v_o INT := 0; v_n INT := 0;
BEGIN
    IF p_glycoct IS NULL OR btrim(p_glycoct) = '' THEN RETURN NULL; END IF;
    -- 去掉立体化学段（"s:1:5" 之类），避免其数字干扰
    v_body := regexp_replace(p_glycoct, 's:[0-9]+:[a-z]+', '', 'g');

    FOREACH v_line IN ARRAY string_to_array(v_body, E'\n') LOOP
        v_line := btrim(v_line);
        IF v_line = '' THEN CONTINUE; END IF;
        IF v_line = 'RES' THEN v_in_res := TRUE;  CONTINUE; END IF;
        IF v_line IN ('LIN','UND') THEN v_in_res := FALSE; CONTINUE; END IF;
        IF NOT v_in_res THEN CONTINUE; END IF;

        -- 取代基残基：2s:n-acetyl
        IF v_line ~ '^[0-9]+s:' THEN
            IF lower(v_line) LIKE '%n-acetyl%' THEN v_nac := v_nac + 1; END IF;
            CONTINUE;
        END IF;
        IF v_line !~ '^[0-9]+b:' THEN CONTINUE; END IF;

        -- 用 superclass（HEX/PEN）作锚点取 stem，避开 basetype 内的冒号
        v_code := substring(v_line,
            'b:[abx]?-?[dl]?-?([A-Za-z]+)-?(?:HEX|PEN|TET|NON)');
        IF v_code IS NULL OR btrim(v_code) = '' THEN CONTINUE; END IF;
        v_code := lower(btrim(v_code));
        v_mods := coalesce(split_part(v_line, '|', 2), '');

        -- 基准一律取"游离单糖（环状）"形式，修饰段在此之上修正。
        -- 烷酮糖/己醛糖：C6H12O6
        IF v_code IN ('glc','gal','man','all','alt','gul','ido','tal','fru') THEN
            v_c := v_c + 6; v_h := v_h + 12; v_o := v_o + 6;
        -- 6-脱氧己糖：C6H12O5
        ELSIF v_code IN ('fuc','rha','qui') THEN
            v_c := v_c + 6; v_h := v_h + 12; v_o := v_o + 5;
        -- 戊醛糖：C5H10O5
        ELSIF v_code IN ('xyl','ara','rib','lyx') THEN
            v_c := v_c + 5; v_h := v_h + 10; v_o := v_o + 5;
        -- 糖醛酸：C6H10O7
        ELSIF v_code IN ('glca','gala','mana','idoa','gula') THEN
            v_c := v_c + 6; v_h := v_h + 10; v_o := v_o + 7;
        ELSE
            RETURN NULL;    -- 未收录的单糖：宁可不判定，也不误报
        END IF;

        -- 修饰段修正（相对上面的基准）
        IF v_mods ~ '6:a' THEN            -- C6 羟甲基 → 羧基：-2H +O
            v_h := v_h - 2; v_o := v_o + 1;
        ELSIF v_mods ~ '6:d' THEN         -- C6 羟甲基 → 甲基：-O（H 数不变）
            v_o := v_o - 1;
        END IF;
    END LOOP;

    -- N-乙酰氨基糖：GlycoCT 写作"己醛糖骨架 + 2s:n-acetyl 取代基"。
    --   游离 GlcNAc(C8H15NO6) − 己醛糖(C6H12O6) = C2H3N：
    --   O2-OH → O2-NHAc 是"脱羟基 + 接乙酰氨基"，净增 C2H3N，
    --   氧数不变（少一个 OH 的 O，多一个 C=O 的 O，相抵）。
    --   故每个乙酰基只补 C2H3N；成键脱水交给下面的 v_nbonds 统一处理
    --   （游离态 → C8H15NO6，链中 → C8H14N... 见用例）。
    v_c := v_c + 2*v_nac;
    v_h := v_h + 3*v_nac;
    v_n := v_n + v_nac;

    -- 糖苷键：只数 LIN 段的糖苷键。n-acetyl 的取代键（…d(2+1)2n）不是糖苷键，
    -- 成键脱的是甲醇而非水，计入会多脱水（实测把 GlcNAc 算成 C8H12NO6）
    v_nbonds := glycan_bond_count(p_glycoct);
    v_h := v_h - 2*v_nbonds;
    v_o := v_o - v_nbonds;

    IF v_c <= 0 THEN RETURN NULL; END IF;
    -- 元素计数为 1 时不写下标（C8H15NO6 而非 C8H15N1O6）
    RETURN 'C' || v_c || 'H' || v_h
           || CASE WHEN v_n > 0 THEN 'N' || CASE WHEN v_n > 1 THEN v_n::text ELSE '' END
              ELSE '' END
           || 'O' || v_o;
END $$;


-- ----------------------------------------------------------------------------
-- 4b. GlycoCT LIN 段 -> 糖苷键数
-- ----------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION glycan_bond_count(p_glycoct TEXT)
RETURNS INT
LANGUAGE plpgsql IMMUTABLE AS $$
DECLARE
    v_in_lin BOOLEAN := FALSE;
    v_line   TEXT;
    v_n      INT := 0;
BEGIN
    IF p_glycoct IS NULL OR btrim(p_glycoct) = '' THEN RETURN 0; END IF;
    FOREACH v_line IN ARRAY string_to_array(p_glycoct, E'\n') LOOP
        v_line := btrim(v_line);
        IF v_line = '' THEN CONTINUE; END IF;
        IF v_line = 'LIN' THEN v_in_lin := TRUE;  CONTINUE; END IF;
        IF v_line = 'RES' THEN v_in_lin := FALSE; CONTINUE; END IF;
        IF NOT v_in_lin OR v_line !~ '^[0-9]+:' THEN CONTINUE; END IF;
        -- N-乙酰取代键形如 "2:1d(2+1)2n"、"3:3d(2+1)4n"，右端是取代基残基。
        -- 不能用 'n$' 判断：糖苷键若是 "1:3o(4+1)1d" 之后的第二条，
        -- 行尾也不是 n，但用 d(2+1)Nn 匹配才准确（实测漏掉一条导致多脱一分子水）。
        IF v_line ~ 'd\([0-9]+\+[0-9]+\)[0-9]+n$' THEN CONTINUE; END IF;
        v_n := v_n + 1;
    END LOOP;
    RETURN v_n;
END $$;


-- ----------------------------------------------------------------------------
-- 6. 元素计数辅助（用于 C4/C9 比较）
-- ----------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION formula_element_counts(p_formula TEXT)
RETURNS JSONB
LANGUAGE plpgsql IMMUTABLE AS $$
DECLARE
    v_m     TEXT[];
    v_elem  TEXT;
    v_cnt   INT;
    v_out   JSONB := '{}'::jsonb;
BEGIN
    IF p_formula IS NULL OR btrim(p_formula) = '' THEN RETURN '{}'::jsonb; END IF;
    -- 逐元素计数：把 "C12H22O11" 拆成 {C:12,H:22,O:11}
    -- 注意 regexp_matches 返回 text[]，用下标取组
    FOR v_m IN SELECT m FROM regexp_matches(upper(btrim(p_formula)),
                                            '([A-Z][a-z]?)(\d*)', 'g') AS m
    LOOP
        v_elem := v_m[1];
        v_cnt  := coalesce(nullif(v_m[2], ''), '1')::int;
        v_out  := jsonb_set(v_out, ARRAY[v_elem],
                            to_jsonb(coalesce((v_out->>v_elem)::int, 0) + v_cnt));
    END LOOP;
    RETURN v_out;
END $$;


-- ----------------------------------------------------------------------------
-- 7. 行内校验：返回发现列表（JSONB 数组）
--    只做自包含判定，绝不依赖派生表 —— 这是它能安全硬拒的前提。
-- ----------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION consistency_check_row(
    p_level   TEXT,
    p_glycoct TEXT,
    p_arch    JSONB,
    p_formula TEXT,
    p_mw      NUMERIC
) RETURNS JSONB
LANGUAGE plpgsql STABLE AS $$
DECLARE
    v_find    JSONB := '[]'::jsonb;
    v_sev     TEXT;
    v_placeholder BOOLEAN := FALSE;
    v_calc    TEXT;
    v_lit     JSONB;
    v_lit_n   INT;
    v_calc_n  INT;
    v_mass    NUMERIC;
BEGIN
    -- 占位键识别（ETL 在结构无法确定时写入的确定性占位串）
    v_placeholder := p_glycoct IS NULL
        OR btrim(p_glycoct) = ''
        OR p_glycoct LIKE 'UNRESOLVED:%'
        OR p_glycoct LIKE 'DOMAIN:%'
        OR p_glycoct LIKE 'COMPOSITION:%';

    -- C7a: structure_level 必须显式给出
    IF p_level IS NULL OR btrim(p_level) = '' THEN
        v_find := v_find || jsonb_build_object(
            'rule','C7a','severity','block',
            'message','structure_level 未给出：无法判断该记录能比对到什么程度');
    END IF;

    -- C7b: 声明 complete/repeat_unit 就必须有真结构编码
    IF p_level IN ('complete','repeat_unit') AND v_placeholder THEN
        v_find := v_find || jsonb_build_object(
            'rule','C7b','severity','block',
            'message', format('structure_level=%s 要求真结构编码，但 glycoct 为占位值/空值 (%s)',
                              p_level, coalesce(left(p_glycoct, 40), 'NULL')));
    END IF;

    -- C8: domain_only 必须给出域架构
    IF p_level = 'domain_only'
       AND (p_arch IS NULL OR coalesce(jsonb_array_length(p_arch->'domains'), 0) = 0) THEN
        v_find := v_find || jsonb_build_object(
            'rule','C8','severity','block',
            'message','structure_level=domain_only 但 domain_architecture.domains 为空');
    END IF;

    -- C10: 占位键与等级自洽
    IF p_level IN ('composition_only','domain_only') AND NOT v_placeholder THEN
        v_find := v_find || jsonb_build_object(
            'rule','C10','severity','block',
            'message', format('structure_level=%s 属"无完整结构"档，但 glycoct 是结构编码',
                              p_level));
    END IF;

    -- C4: 分子式 ↔ GlycoCT 推导分子式（两条独立路径互证）
    IF p_formula IS NOT NULL AND btrim(p_formula) <> '' AND NOT v_placeholder THEN
        v_calc := glycan_formula_from_glycoct(p_glycoct);
        IF v_calc IS NOT NULL THEN
            v_lit := formula_element_counts(p_formula);
            v_lit_n := coalesce((v_lit->>'C')::int, -1);
            v_calc_n := (formula_element_counts(v_calc)->>'C')::int;
            IF v_lit_n > 0 AND v_lit_n <> v_calc_n THEN
                v_find := v_find || jsonb_build_object(
                    'rule','C4','severity','warn',
                    'message', format('分子式对账不符：文献 %s (C%d) vs GlycoCT 推导 %s (C%d)',
                                      p_formula, v_lit_n, v_calc, v_calc_n),
                    'details', jsonb_build_object('literature', p_formula,
                                                  'from_glycoct', v_calc));
            END IF;
        END IF;
    END IF;

    -- C9: 分子量 ↔ 分子式推导质量（偏离 10% 以上告警）
    IF p_mw IS NOT NULL AND NOT v_placeholder THEN
        v_calc := glycan_formula_from_glycoct(p_glycoct);
        IF v_calc IS NOT NULL THEN
            v_lit := formula_element_counts(v_calc);
            v_mass := 12.011 * coalesce((v_lit->>'C')::numeric, 0)
                    + 1.008  * coalesce((v_lit->>'H')::numeric, 0)
                    + 15.999 * coalesce((v_lit->>'O')::numeric, 0)
                    + 14.007 * coalesce((v_lit->>'N')::numeric, 0);
            IF v_mass > 0 AND abs(p_mw - v_mass) / v_mass > 0.10 THEN
                v_find := v_find || jsonb_build_object(
                    'rule','C9','severity','warn',
                    'message', format('分子量 %s Da 与分子式 %s 推导质量 %.1f Da 偏离超过 10%%',
                                      p_mw, v_calc, v_mass),
                    'details', jsonb_build_object('literature_mw', p_mw,
                                                  'formula_mass', round(v_mass, 2),
                                                  'formula', v_calc));
            END IF;
        END IF;
    END IF;

    RETURN v_find;
END $$;


-- ----------------------------------------------------------------------------
-- 8. 行内发现暂存（避免被拒记录留下脏发现）
-- ----------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS _dsh_consistency_pending (
    sugar_id BIGINT,
    rule     TEXT,
    severity TEXT,
    message  TEXT,
    details  JSONB
);


CREATE OR REPLACE FUNCTION consistency_gate_sugars() RETURNS TRIGGER
LANGUAGE plpgsql AS $$
DECLARE
    v_find  JSONB;
    v_item  JSONB;
    v_rule  TEXT;
    v_sev   TEXT;
    v_msg   TEXT;
    v_blocks TEXT[] := ARRAY[]::TEXT[];
BEGIN
    v_find := consistency_check_row(NEW.structure_level, NEW.glycoct,
                                    NEW.domain_architecture, NEW.molecular_formula,
                                    NEW.molecular_weight);

    FOR v_item IN SELECT * FROM jsonb_array_elements(v_find) LOOP
        v_rule := v_item->>'rule';
        v_msg  := v_item->>'message';
        -- 严重度以策略表为准（可在不改代码的前提下调整）
        SELECT p.severity INTO v_sev FROM consistency_rule_policy p WHERE p.rule = v_rule;
        v_sev := coalesce(v_sev, v_item->>'severity', 'warn');

        IF v_sev = 'block' THEN
            v_blocks := v_blocks || format('[%s] %s', v_rule, v_msg);
        END IF;
        INSERT INTO _dsh_consistency_pending(sugar_id, rule, severity, message, details)
        VALUES (NEW.sugar_id, v_rule, v_sev, v_msg, v_item->'details');
    END LOOP;

    IF array_length(v_blocks, 1) > 0 THEN
        RAISE EXCEPTION E'一致性硬门槛拒绝写入 sugars:\n  %',
            array_to_string(v_blocks, E'\n  ')
            USING ERRCODE = 'check_violation',
                  HINT = '修正数据，或由 DBA 调整 consistency_rule_policy 中对应规则的 severity';
    END IF;

    -- 内容变化（结构/等级/域架构）后视为未复核，fail-closed
    IF TG_OP = 'UPDATE'
       AND (NEW.glycoct IS DISTINCT FROM OLD.glycoct
            OR NEW.structure_level IS DISTINCT FROM OLD.structure_level
            OR NEW.domain_architecture IS DISTINCT FROM OLD.domain_architecture
            OR NEW.composition IS DISTINCT FROM OLD.composition) THEN
        NEW.consistency_ok := FALSE;
    END IF;
    RETURN NEW;
END $$;


CREATE OR REPLACE FUNCTION consistency_flush_pending() RETURNS TRIGGER
LANGUAGE plpgsql AS $$
BEGIN
    INSERT INTO sugar_consistency_findings(sugar_id, rule, severity, phase, message, details)
    SELECT sugar_id, rule, severity, 'row', message, details
    FROM _dsh_consistency_pending WHERE sugar_id = NEW.sugar_id;
    DELETE FROM _dsh_consistency_pending WHERE sugar_id = NEW.sugar_id;
    RETURN NULL;
END $$;


DROP TRIGGER IF EXISTS trg_consistency_gate ON sugars;
CREATE TRIGGER trg_consistency_gate
    BEFORE INSERT OR UPDATE OF glycoct, structure_level, domain_architecture
    ON sugars FOR EACH ROW EXECUTE FUNCTION consistency_gate_sugars();

DROP TRIGGER IF EXISTS trg_consistency_flush ON sugars;
CREATE TRIGGER trg_consistency_flush
    AFTER INSERT OR UPDATE ON sugars FOR EACH ROW EXECUTE FUNCTION consistency_flush_pending();


-- ----------------------------------------------------------------------------
-- 9. 跨表复核（关卡二）：把 C1-C6 这类需要派生数据的规则算清楚
--    返回发现数量；无 block 级问题时把 consistency_ok 置 TRUE。
-- ----------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION consistency_recheck(p_sugar_id BIGINT)
RETURNS INT
LANGUAGE plpgsql AS $$
DECLARE
    v_rec        RECORD;
    v_find       JSONB := '[]'::jsonb;
    v_item       JSONB;
    v_rule       TEXT;
    v_sev        TEXT;
    v_msg        TEXT;
    v_n          INT := 0;
    v_blocks     INT := 0;
    v_pct_sum    NUMERIC;
    v_pct_max    NUMERIC;
    v_pct_monos  TEXT[];
    v_res_monos  TEXT[];
    v_only_pct   TEXT[];
    v_only_res   TEXT[];
    v_meth_n     INT;
    v_res_n      INT;
BEGIN
    SELECT sugar_id, glycoct, structure_level, composition, iupac_short,
           molecular_formula, molecular_weight, domain_architecture
      INTO v_rec FROM sugars WHERE sugar_id = p_sugar_id;
    IF NOT FOUND THEN
        RAISE EXCEPTION 'sugar_id=% 不存在', p_sugar_id;
    END IF;

    -- 先清掉上一轮复核发现（保留 row 阶段的记录）
    DELETE FROM sugar_consistency_findings
     WHERE sugar_id = p_sugar_id AND phase = 'recheck';

    -- C1: 组成百分比合计
    SELECT sum((value)::text::numeric), max((value)::text::numeric)
      INTO v_pct_sum, v_pct_max
      FROM polysaccharide_props pp, jsonb_each(pp.monosaccharide_ratio)
     WHERE pp.sugar_id = p_sugar_id AND jsonb_typeof(pp.monosaccharide_ratio) = 'object';

    -- 同 Python 侧：monosaccharide_ratio 可能是"摩尔比"而非百分比
    -- （{"Fru":1.0}、{"Glc":1,"Gal":2}）。比值口径下合计很小，
    -- 不做百分比判定，避免把合法数据误报为"漏抽行"（实测在 Inulin 上误报过）。
    IF v_pct_sum IS NOT NULL
       AND NOT (v_pct_max <= 8 AND v_pct_sum <= 30)
       AND abs(v_pct_sum - 100) > 5 THEN
        v_find := v_find || jsonb_build_object('rule','C1','severity','warn',
            'message', format('单糖组成百分比合计 %s%%，偏离 100%% 超过 5%%', round(v_pct_sum, 2)),
            'details', jsonb_build_object('total', round(v_pct_sum, 3)));
    END IF;

    -- C2: 组成表 ↔ 残基表 单糖种类
    SELECT array_agg(DISTINCT upper(key)) INTO v_pct_monos
      FROM polysaccharide_props pp, jsonb_each(pp.monosaccharide_ratio)
     WHERE pp.sugar_id = p_sugar_id AND jsonb_typeof(pp.monosaccharide_ratio) = 'object';
    SELECT array_agg(DISTINCT upper(monosaccharide_name)) INTO v_res_monos
      FROM residues WHERE sugar_id = p_sugar_id;

    IF v_pct_monos IS NOT NULL AND v_res_monos IS NOT NULL THEN
        SELECT array_agg(x) INTO v_only_pct
          FROM unnest(v_pct_monos) x WHERE x <> ALL(v_res_monos);
        SELECT array_agg(x) INTO v_only_res
          FROM unnest(v_res_monos) x WHERE x <> ALL(v_pct_monos);
        IF v_only_pct IS NOT NULL OR v_only_res IS NOT NULL THEN
            v_find := v_find || jsonb_build_object('rule','C2','severity','warn',
                'message', format('单糖种类不匹配：仅组成表有 %s；仅残基表有 %s',
                                  coalesce(array_to_string(v_only_pct, ','), '-'),
                                  coalesce(array_to_string(v_only_res, ','), '-')),
                'details', jsonb_build_object('only_in_percentages', to_jsonb(v_only_pct),
                                              'only_in_residues', to_jsonb(v_only_res)));
        END IF;
    END IF;

    -- C3: 甲基化分析 ↔ 残基表
    IF v_res_monos IS NOT NULL THEN
        SELECT count(*) INTO v_meth_n FROM polysaccharide_props pp
         WHERE pp.sugar_id = p_sugar_id AND coalesce(pp.branching,'') <> '';
        SELECT count(*) INTO v_res_n FROM residues WHERE sugar_id = p_sugar_id;

        -- 只在甲基化字符串里出现、残基表没有的单糖（粗略按词形匹配）
        IF v_meth_n > 0 THEN
            SELECT array_agg(DISTINCT m[1]) INTO v_only_pct
              FROM polysaccharide_props pp,
                   regexp_matches(pp.branching, '([A-Za-z]+A?)\(', 'g') AS m
             WHERE pp.sugar_id = p_sugar_id
               AND upper(m[1]) <> ALL(v_res_monos);
            IF v_only_pct IS NOT NULL THEN
                v_find := v_find || jsonb_build_object('rule','C3','severity','warn',
                    'message', format('甲基化分析出现残基表未列出的单糖 %s',
                                      array_to_string(v_only_pct, ',')),
                    'details', jsonb_build_object('only_in_methylation', to_jsonb(v_only_pct)));
            END IF;
        END IF;

        -- C6: 连接类型数与残基类型数应在同一量级
        IF v_res_n > 0 THEN
            SELECT count(*) INTO v_meth_n
              FROM polysaccharide_props pp,
                   regexp_split_to_table(coalesce(pp.branching,''), ',\s+(?=[A-Za-z0-9(])') AS t
             WHERE pp.sugar_id = p_sugar_id AND btrim(t) <> '';
            IF v_meth_n > 0 AND greatest(v_meth_n, v_res_n)::numeric
                               / least(v_meth_n, v_res_n) >= 2.5 THEN
                v_find := v_find || jsonb_build_object('rule','C6','severity','warn',
                    'message', format('甲基化连接类型 %s 种 vs 残基 %s 个，相差 ≥2.5 倍',
                                      v_meth_n, v_res_n));
            END IF;
        END IF;
    END IF;

    -- 落库 + 统计
    FOR v_item IN SELECT * FROM jsonb_array_elements(v_find) LOOP
        v_rule := v_item->>'rule';
        v_msg  := v_item->>'message';
        SELECT p.severity INTO v_sev FROM consistency_rule_policy p WHERE p.rule = v_rule;
        v_sev := coalesce(v_sev, v_item->>'severity', 'warn');
        INSERT INTO sugar_consistency_findings(sugar_id, rule, severity, phase, message, details)
        VALUES (p_sugar_id, v_rule, v_sev, 'recheck', v_msg, v_item->'details');
        v_n := v_n + 1;
        IF v_sev = 'block' THEN v_blocks := v_blocks + 1; END IF;
    END LOOP;

    -- 还有 row 阶段的 block 级发现时同样不予通过
    SELECT count(*) INTO v_blocks FROM sugar_consistency_findings
     WHERE sugar_id = p_sugar_id AND severity = 'block';

    UPDATE sugars SET consistency_ok = (v_blocks = 0), updated_at = now()
     WHERE sugar_id = p_sugar_id;

    RETURN v_n;
END $$;


-- ----------------------------------------------------------------------------
-- 10. 复核队列视图：下游只查 consistency_ok 即可，这里给运维看"为什么没通过"
-- ----------------------------------------------------------------------------
CREATE OR REPLACE VIEW v_sugar_consistency_review AS
SELECT s.sugar_id,
       s.iupac_short,
       s.structure_level,
       s.consistency_ok,
       s.first_seen_doi,
       count(f.finding_id)                                          AS n_findings,
       count(f.finding_id) FILTER (WHERE f.severity = 'block')      AS n_block,
       string_agg(DISTINCT f.rule, ',' ORDER BY f.rule)             AS rules,
       max(f.detected_at)                                           AS last_detected
FROM sugars s
LEFT JOIN sugar_consistency_findings f ON f.sugar_id = s.sugar_id
GROUP BY s.sugar_id, s.iupac_short, s.structure_level, s.consistency_ok, s.first_seen_doi;

COMMENT ON VIEW v_sugar_consistency_review IS
'一致性复核队列：consistency_ok=FALSE 的记录尚未通过复核（fail-closed），下游不应采信。';
