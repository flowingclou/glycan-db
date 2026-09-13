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
INSERT INTO sugars (sugar_type, glycoct, glycoct_hash, iupac_short, molecular_formula, molecular_weight, anomer, structure_confidence, stereochemistry_defined, first_seen_doi)
SELECT 'oligo',
       'RES 1b:a-dglcp-1:5(1:4)|2:x,1a:b-dglcp-1:5|2:x',
       encode(digest('RES 1b:a-dglcp-1:5(1:4)|2:x,1a:b-dglcp-1:5|2:x','sha256'),'hex'),
       'α-D-Glcp-(1→4)-D-Glcp', 'C12H22O11', 342.2965, 'a',
       'confirmed_2d', TRUE, '10.1021/acs.joc.0c00000'
WHERE NOT EXISTS (SELECT 1 FROM sugars WHERE iupac_short = 'α-D-Glcp-(1→4)-D-Glcp');

-- 麦芽糖残基组成（还原端 Glc-b, 非还原端 Glc-a, 连接 1→4）
INSERT INTO residues (sugar_id, residue_seq, monosaccharide_name, ring_form, anomer, is_reducing_end, parent_carbon, linkage_branch)
SELECT s.sugar_id, 1, 'Glc', 'p', 'b', TRUE,  NULL, 0 FROM sugars s WHERE s.iupac_short = 'α-D-Glcp-(1→4)-D-Glcp';
INSERT INTO residues (sugar_id, residue_seq, monosaccharide_name, ring_form, anomer, is_reducing_end, parent_carbon, linkage_branch)
SELECT s.sugar_id, 2, 'Glc', 'p', 'a', FALSE, 4, 0 FROM sugars s WHERE s.iupac_short = 'α-D-Glcp-(1→4)-D-Glcp';

-- 麦芽糖 ¹H 实验（异头区）: 非还原端 H1' 5.41 (d, J=3.8); 还原端 H1α 5.23 / H1β 4.66
INSERT INTO nmr_experiments (sugar_id, source_id, nmr_type, solvent, frequency_mhz, temperature_c, ph, qc_status, assignment_level)
SELECT s.sugar_id, l.source_id, '1H', 'D2O', 500.0, 25.0, 7.0, 'passed', 'full'
FROM sugars s, literature l
WHERE s.iupac_short = 'α-D-Glcp-(1→4)-D-Glcp' AND l.doi = '10.1021/acs.joc.0c00000';

INSERT INTO nmr_shifts_1d (experiment_id, nucleus, shift_ppm, multiplicity, j_coupling_hz, integration, assignment_position, is_anomeric)
SELECT e.experiment_id, '1H', 5.410, 'd', 3.8, 1.0, 'H1'' (Glc-a)', TRUE
FROM nmr_experiments e JOIN sugars s ON e.sugar_id = s.sugar_id
WHERE s.iupac_short = 'α-D-Glcp-(1→4)-D-Glcp' AND e.nmr_type = '1H';

-- 麦芽糖 2D 相关峰（HSQC 异头区 + HMBC 糖苷键证据 H1'→C4）
INSERT INTO nmr_experiments (sugar_id, source_id, nmr_type, solvent, frequency_mhz, temperature_c, ph, qc_status, assignment_level)
SELECT s.sugar_id, l.source_id, '2D', 'D2O', 500.0, 25.0, 7.0, 'passed', 'full'
FROM sugars s, literature l
WHERE s.iupac_short = 'α-D-Glcp-(1→4)-D-Glcp' AND l.doi = '10.1021/acs.joc.0c00000';

-- HSQC: H1'(5.41) - C1'(100.6)
INSERT INTO nmr_correlations_2d (experiment_id, experiment_2d, proton_shift_ppm, hetero_shift_ppm, atom_pair, residue_from, residue_to, linkage_evidence, intensity)
SELECT e.experiment_id, 'HSQC', 5.410, 100.600, 'H1''-C1''', 'Glc-a', 'Glc-a', FALSE, 1.00
FROM nmr_experiments e JOIN sugars s ON e.sugar_id = s.sugar_id
WHERE s.iupac_short = 'α-D-Glcp-(1→4)-D-Glcp' AND e.nmr_type = '2D';

-- HMBC: H1'(5.41) → C4(Glc-b)  = 糖苷键连接证据
INSERT INTO nmr_correlations_2d (experiment_id, experiment_2d, proton_shift_ppm, hetero_shift_ppm, atom_pair, residue_from, residue_to, linkage_evidence, intensity)
SELECT e.experiment_id, 'HMBC', 5.410, 78.300, 'H1''-C4', 'Glc-a', 'Glc-b', TRUE, 0.85
FROM nmr_experiments e JOIN sugars s ON e.sugar_id = s.sugar_id
WHERE s.iupac_short = 'α-D-Glcp-(1→4)-D-Glcp' AND e.nmr_type = '2D';

-- 麦芽糖物化性质
INSERT INTO physicochemical (sugar_id, optical_rotation, optical_rotation_condition, melting_point_c, solubility, source_id)
SELECT s.sugar_id, 130.5, 'H2O, c=1', 102.0, '易溶于水', l.source_id
FROM sugars s, literature l
WHERE s.iupac_short = 'α-D-Glcp-(1→4)-D-Glcp' AND l.doi = '10.1021/acs.joc.0c00000';

-- ----------------------------------------------------------------------------
-- 4. 多糖示例：菊粉 Inulin  β-D-Fruf-(2→1)- (重复单元)
-- ----------------------------------------------------------------------------
INSERT INTO sugars (sugar_type, glycoct, glycoct_hash, iupac_short, molecular_formula, molecular_weight, anomer, structure_confidence, stereochemistry_defined, first_seen_doi)
SELECT 'poly',
       'RES 1b:b-dfruf-2:1|2:6(1:2)',
       encode(digest('RES 1b:b-dfruf-2:1|2:6(1:2)','sha256'),'hex'),
       'β-D-Fruf-(2→1)-[Inulin]', NULL, NULL, 'b',
       'confirmed_1d', TRUE, '10.1016/j.carres.2019.107800'
WHERE NOT EXISTS (SELECT 1 FROM sugars WHERE iupac_short LIKE 'β-D-Fruf-(2→1)-[Inulin]%');

-- 菊粉重复单元残基（果糖, 呋喃型, β, 连接 2→1）
INSERT INTO residues (sugar_id, residue_seq, monosaccharide_name, ring_form, anomer, is_reducing_end, parent_carbon, linkage_branch)
SELECT s.sugar_id, 1, 'Fru', 'f', 'b', FALSE, 2, 1 FROM sugars s WHERE s.iupac_short = 'β-D-Fruf-(2→1)-[Inulin]';

-- 菊粉多糖专属性质
INSERT INTO polysaccharide_props (sugar_id, repeat_unit_formula, degree_of_polymerization, molecular_weight_mn, molecular_weight_mw, polydispersity, monosaccharide_ratio, backbone, branching)
SELECT s.sugar_id, 'C6H10O5', 30.0, 4900.0, 5200.0, 1.06, '{"Fru":1.0}', 'β-D-Fruf-(2→1)- 线性', '无分支'
FROM sugars s WHERE s.iupac_short = 'β-D-Fruf-(2→1)-[Inulin]';

-- 菊粉 ¹³C 实验（残基级归属: C2 104.2, C3 78.0, C4 75.5 等）
INSERT INTO nmr_experiments (sugar_id, source_id, nmr_type, solvent, frequency_mhz, temperature_c, ph, qc_status, assignment_level)
SELECT s.sugar_id, l.source_id, '13C', 'D2O', 500.0, 25.0, 7.0, 'passed', 'residue_level'
FROM sugars s, literature l
WHERE s.iupac_short = 'β-D-Fruf-(2→1)-[Inulin]' AND l.doi = '10.1016/j.carres.2019.107800';

INSERT INTO nmr_shifts_1d (experiment_id, nucleus, shift_ppm, multiplicity, j_coupling_hz, integration, assignment_position, is_anomeric)
SELECT e.experiment_id, '13C', 104.200, 's', NULL, NULL, 'C2 (Fruf)', TRUE
FROM nmr_experiments e JOIN sugars s ON e.sugar_id = s.sugar_id
WHERE s.iupac_short = 'β-D-Fruf-(2→1)-[Inulin]' AND e.nmr_type = '13C';

INSERT INTO nmr_shifts_1d (experiment_id, nucleus, shift_ppm, multiplicity, j_coupling_hz, integration, assignment_position, is_anomeric)
SELECT e.experiment_id, '13C', 78.000, 's', NULL, NULL, 'C3 (Fruf)', FALSE
FROM nmr_experiments e JOIN sugars s ON e.sugar_id = s.sugar_id
WHERE s.iupac_short = 'β-D-Fruf-(2→1)-[Inulin]' AND e.nmr_type = '13C';

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
