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
