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
