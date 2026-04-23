use std::collections::BTreeMap;

#[derive(Debug, Clone)]
pub struct RawSumstats {
    pub chrom_tokens: Vec<String>,
    pub pos_tokens: Vec<String>,
    pub p_tokens: Vec<String>,
    pub mlog10p_tokens: Vec<String>,
    pub has_p_column: bool,
}

#[derive(Debug, Clone)]
pub struct CleanSumstats {
    pub chr: Vec<i32>,
    pub pos: Vec<f64>,
    pub p: Vec<f64>,
    pub mlog10p: Vec<f64>,
    /// Max POS per chromosome **before** `--skip` filtering (matches Python `df.attrs["chr_max_pos"]`).
    pub chr_max_pos: BTreeMap<i32, f64>,
}

fn normalize_chr_token(value: &str) -> Option<i32> {
    let mut token = value.trim().to_uppercase();
    if token.starts_with("CHR") {
        token = token[3..].to_string();
    }
    match token.as_str() {
        "X" => return Some(23),
        "Y" => return Some(24),
        "MT" | "M" => return Some(25),
        _ => {}
    }
    let iv: i32 = token.parse().ok()?;
    if iv < 1 {
        return None;
    }
    Some(iv)
}

fn parse_f64(s: &str) -> Option<f64> {
    let t = s.trim();
    if t.is_empty() {
        return None;
    }
    t.parse::<f64>().ok()
}

/// Convert CHR/POS/P or CHR/POS/MLOG10P and compute mlog10p; attach chr_max_pos metadata.
pub fn preprocess_sumstats(raw: RawSumstats, skip: f64) -> Result<CleanSumstats, String> {
    let n = raw.chrom_tokens.len();
    if raw.pos_tokens.len() != n {
        return Err("row count mismatch".to_string());
    }

    let mut chr = Vec::new();
    let mut pos = Vec::new();
    let mut p = Vec::new();
    let mut mlog10p = Vec::new();

    if raw.has_p_column {
        if raw.p_tokens.len() != n {
            return Err("row count mismatch".to_string());
        }
        for i in 0..n {
            let c = normalize_chr_token(&raw.chrom_tokens[i]);
            let po = parse_f64(&raw.pos_tokens[i]);
            let pv = parse_f64(&raw.p_tokens[i]);
            let (Some(c), Some(po), Some(pv)) = (c, po, pv) else {
                continue;
            };
            if !(pv > 0.0 && pv <= 1.0) {
                continue;
            }
            chr.push(c);
            pos.push(po);
            p.push(pv);
            mlog10p.push(-pv.log10());
        }
    } else {
        if raw.mlog10p_tokens.len() != n {
            return Err("row count mismatch".to_string());
        }
        for i in 0..n {
            let c = normalize_chr_token(&raw.chrom_tokens[i]);
            let po = parse_f64(&raw.pos_tokens[i]);
            let ml = parse_f64(&raw.mlog10p_tokens[i]);
            let (Some(c), Some(po), Some(ml)) = (c, po, ml) else {
                continue;
            };
            let pv = 10f64.powf(-ml);
            if !(pv > 0.0 && pv <= 1.0) {
                continue;
            }
            chr.push(c);
            pos.push(po);
            p.push(pv);
            mlog10p.push(ml);
        }
    }

    if chr.is_empty() {
        return Ok(CleanSumstats {
            chr: vec![],
            pos: vec![],
            p: vec![],
            mlog10p: vec![],
            chr_max_pos: BTreeMap::new(),
        });
    }

    let mut chr_max_pos: BTreeMap<i32, f64> = BTreeMap::new();
    for i in 0..chr.len() {
        let c = chr[i];
        let po = pos[i];
        chr_max_pos
            .entry(c)
            .and_modify(|m| *m = (*m).max(po))
            .or_insert(po);
    }

    let mut keep: Vec<bool> = vec![true; chr.len()];
    if skip > 0.0 {
        for i in 0..chr.len() {
            if mlog10p[i] < skip {
                keep[i] = false;
            }
        }
    }

    let mut idx: Vec<usize> = (0..chr.len()).collect();
    idx.sort_by(|&a, &b| {
        chr[a]
            .cmp(&chr[b])
            .then_with(|| pos[a].partial_cmp(&pos[b]).unwrap_or(std::cmp::Ordering::Equal))
    });

    let mut chr2 = Vec::new();
    let mut pos2 = Vec::new();
    let mut p2 = Vec::new();
    let mut m2 = Vec::new();
    for &i in &idx {
        if keep[i] {
            chr2.push(chr[i]);
            pos2.push(pos[i]);
            p2.push(p[i]);
            m2.push(mlog10p[i]);
        }
    }

    Ok(CleanSumstats {
        chr: chr2,
        pos: pos2,
        p: p2,
        mlog10p: m2,
        chr_max_pos,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;
    use crate::io::load_sumstats;

    fn fixture(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/fixtures")
            .join(name)
    }

    #[test]
    fn preprocess_has_required_columns() {
        let raw = load_sumstats(&fixture("sumstats_small.tsv"), "\t", None, None, None, None).unwrap();
        let c = preprocess_sumstats(raw, 3.0).unwrap();
        assert!(!c.chr.is_empty());
        assert_eq!(c.chr.len(), c.pos.len());
    }

    #[test]
    fn chr_order_sorted() {
        let raw = load_sumstats(&fixture("sumstats_small.tsv"), "\t", None, None, None, None).unwrap();
        let c = preprocess_sumstats(raw, 3.0).unwrap();
        let mut prev = 0;
        for &x in &c.chr {
            assert!(x >= prev);
            prev = x;
        }
        // MT row is below skip=3.0 in the small fixture, so max remaining chromosome is X (23).
        assert_eq!(*c.chr.iter().max().unwrap(), 23);
    }
}
