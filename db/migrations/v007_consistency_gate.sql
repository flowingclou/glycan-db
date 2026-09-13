-- ============================================================================
-- 糖类数据库 · V7 迁移（一致性规则落库为**硬门槛**）
-- 前置: 需先执行 v001 → v002 → v003 → v004 → v005 → v006
-- 数据库: PostgreSQL 16+
-- ============================================================================
-- 背景:
--   v006 之前，交叉一致性校验只跑在 Python 侧（glycan_etl/consistency.py），
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


-- ----------------------------------------------------------------------------
-- 11. 验证
-- ----------------------------------------------------------------------------
-- 11.1 规则策略:
-- SELECT * FROM consistency_rule_policy ORDER BY rule;
--
-- 11.2 复核队列:
-- SELECT * FROM v_sugar_consistency_review ORDER BY consistency_ok, sugar_id;
--
-- 11.3 硬门槛自测（应被拒绝）:
-- INSERT INTO sugars (sugar_type, glycoct, glycoct_hash, iupac_short, structure_level)
-- VALUES ('poly','UNRESOLVED:x:1:abc', encode(digest('UNRESOLVED:x:1:abc','sha256'),'hex'),
--         'bad','complete');          -- C7b: 占位键不得配 complete
--
-- 11.4 全库复核:
-- SELECT consistency_recheck(sugar_id) FROM sugars ORDER BY sugar_id;
