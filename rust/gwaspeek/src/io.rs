use flate2::read::GzDecoder;
use serde::Deserialize;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

use crate::preprocess::RawSumstats;

fn is_gzip_path(path: &Path) -> bool {
    path.to_string_lossy()
        .as_ref()
        .to_ascii_lowercase()
        .ends_with(".gz")
}

const ALIAS_JSON: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../src/gwaspeek/data/column_aliases_from_formatbook.json"
));

pub fn repo_data_dir() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../src/gwaspeek/data")
}

pub fn normalize_sep(sep: &str) -> String {
    match sep {
        r"\t" => "\t".to_string(),
        r"\n" => "\n".to_string(),
        r"\r" => "\r".to_string(),
        _ => sep.to_string(),
    }
}

#[derive(Debug, Deserialize)]
struct AliasFile {
    aliases: HashMap<String, Vec<String>>,
}

fn column_alias_table() -> &'static AliasFile {
    static TABLE: std::sync::OnceLock<AliasFile> = std::sync::OnceLock::new();
    TABLE.get_or_init(|| serde_json::from_str(ALIAS_JSON).expect("valid column alias JSON"))
}

fn header_lower_index(headers: &[String]) -> HashMap<String, String> {
    let mut out = HashMap::new();
    for h in headers {
        let key = h.trim().to_lowercase();
        out.entry(key).or_insert_with(|| h.clone());
    }
    out
}

fn resolve_explicit(
    col: &str,
    headers: &[String],
    lower_index: &HashMap<String, String>,
) -> Result<String, String> {
    if headers.iter().any(|h| h == col) {
        return Ok(col.to_string());
    }
    let key = col.trim().to_lowercase();
    if let Some(found) = lower_index.get(&key) {
        return Ok(found.clone());
    }
    let preview: Vec<String> = headers.iter().take(12).map(|c| format!("{c:?}")).collect();
    Err(format!(
        "Column {col:?} not found in input header (first columns: {}…)",
        preview.join(", ")
    ))
}

fn pick_alias(
    headers_set: &HashSet<String>,
    lower_index: &HashMap<String, String>,
    aliases: &[String],
) -> Option<String> {
    let mut ranked: Vec<&String> = aliases.iter().collect();
    ranked.sort_by(|a, b| {
        let la = a.len();
        let lb = b.len();
        lb.cmp(&la).then_with(|| a.to_lowercase().cmp(&b.to_lowercase()))
    });
    for alias in ranked {
        let a = alias.as_str();
        if headers_set.contains(a) {
            return Some(a.to_string());
        }
        let key = a.trim().to_lowercase();
        if let Some(h) = lower_index.get(&key) {
            return Some(h.clone());
        }
    }
    None
}

/// Returns (chrom_src, pos_src, p_src, mlog10p_src) where exactly one of (p_src, mlog10p_src) is Some.
pub fn detect_sumstat_columns(
    headers: &[String],
    chrom_col: Option<&str>,
    pos_col: Option<&str>,
    p_col: Option<&str>,
    mlog10p_col: Option<&str>,
) -> Result<(String, String, Option<String>, Option<String>), String> {
    let spec = column_alias_table();
    let raw_aliases = &spec.aliases;
    let aliases_chr: Vec<String> = raw_aliases
        .get("CHR")
        .cloned()
        .unwrap_or_else(|| vec!["CHR".to_string()]);
    let aliases_pos: Vec<String> = raw_aliases
        .get("POS")
        .cloned()
        .unwrap_or_else(|| vec!["POS".to_string()]);
    let aliases_p: Vec<String> = raw_aliases
        .get("P")
        .cloned()
        .unwrap_or_else(|| vec!["P".to_string()]);
    let aliases_m: Vec<String> = raw_aliases
        .get("MLOG10P")
        .cloned()
        .unwrap_or_else(|| vec!["MLOG10P".to_string()]);

    let headers_list: Vec<String> = headers.iter().cloned().collect();
    let headers_set: HashSet<String> = headers_list.iter().cloned().collect();
    let lower_index = header_lower_index(&headers_list);

    let chrom_src = if let Some(c) = chrom_col {
        resolve_explicit(c, &headers_list, &lower_index)?
    } else {
        pick_alias(&headers_set, &lower_index, &aliases_chr)
            .ok_or_else(|| "Could not auto-detect chromosome column (CHR). Set --chr.".to_string())?
    };

    let pos_src = if let Some(c) = pos_col {
        resolve_explicit(c, &headers_list, &lower_index)?
    } else {
        pick_alias(&headers_set, &lower_index, &aliases_pos)
            .ok_or_else(|| "Could not auto-detect position column (POS). Set --pos.".to_string())?
    };

    let mut p_src: Option<String> = None;
    let mut mlog_src: Option<String> = None;
    if let Some(c) = p_col {
        p_src = Some(resolve_explicit(c, &headers_list, &lower_index)?);
    }
    if let Some(c) = mlog10p_col {
        mlog_src = Some(resolve_explicit(c, &headers_list, &lower_index)?);
    }

    if p_src.is_none() && mlog_src.is_none() {
        p_src = pick_alias(&headers_set, &lower_index, &aliases_p);
        if p_src.is_none() {
            mlog_src = pick_alias(&headers_set, &lower_index, &aliases_m);
        }
        if p_src.is_none() && mlog_src.is_none() {
            return Err(
                "Could not auto-detect P-value or -log10(P) column. Set --p and/or --mlog10p."
                    .to_string(),
            );
        }
    } else if p_src.is_some() && mlog_src.is_some() {
        return Err("Specify only one of --p or --mlog10p, not both.".to_string());
    }

    Ok((chrom_src, pos_src, p_src, mlog_src))
}

pub fn load_sumstats(
    path: &Path,
    sep: &str,
    chrom_col: Option<&str>,
    pos_col: Option<&str>,
    p_col: Option<&str>,
    mlog10p_col: Option<&str>,
) -> Result<RawSumstats, String> {
    let normalized_sep = normalize_sep(sep);
    let delim = if normalized_sep.len() == 1 {
        normalized_sep.chars().next().unwrap()
    } else {
        return Err("Multi-character separators are not supported.".to_string());
    };

    let f = File::open(path).map_err(|e| e.to_string())?;
    let reader: Box<dyn Read> = if is_gzip_path(path) {
        Box::new(BufReader::new(GzDecoder::new(f)))
    } else {
        Box::new(BufReader::new(f))
    };
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(delim as u8)
        .has_headers(true)
        .flexible(true)
        .from_reader(reader);

    let headers: Vec<String> = rdr
        .headers()
        .map_err(|e| e.to_string())?
        .iter()
        .map(|s| s.to_string())
        .collect();

    let (chrom_src, pos_src, p_src, mlog_src) =
        detect_sumstat_columns(&headers, chrom_col, pos_col, p_col, mlog10p_col)?;

    let chrom_i = headers
        .iter()
        .position(|h| h == &chrom_src)
        .ok_or_else(|| "chrom column missing".to_string())?;
    let pos_i = headers
        .iter()
        .position(|h| h == &pos_src)
        .ok_or_else(|| "pos column missing".to_string())?;
    let stat_i = if let Some(ref p) = p_src {
        headers
            .iter()
            .position(|h| h == p)
            .ok_or_else(|| "p column missing".to_string())?
    } else {
        let m = mlog_src.as_ref().unwrap();
        headers
            .iter()
            .position(|h| h == m)
            .ok_or_else(|| "mlog10p column missing".to_string())?
    };

    let mode_p = p_src.is_some();
    let mut chrom = Vec::new();
    let mut pos = Vec::new();
    let mut p_vals = Vec::new();
    let mut mlog_vals = Vec::new();

    for rec in rdr.records() {
        let rec = rec.map_err(|e| e.to_string())?;
        let c = rec.get(chrom_i).unwrap_or("").to_string();
        let p = rec.get(pos_i).unwrap_or("").to_string();
        let s = rec.get(stat_i).unwrap_or("").to_string();
        chrom.push(c);
        pos.push(p);
        if mode_p {
            p_vals.push(s);
            mlog_vals.push(String::new());
        } else {
            p_vals.push(String::new());
            mlog_vals.push(s);
        }
    }

    Ok(RawSumstats {
        chrom_tokens: chrom,
        pos_tokens: pos,
        p_tokens: p_vals,
        mlog10p_tokens: mlog_vals,
        has_p_column: mode_p,
    })
}

// In-memory variant (used by tests); file-based loader is `load_sumstats`.
pub fn load_sumstats_from_reader<R: std::io::Read>(
    r: R,
    sep: &str,
    chrom_col: Option<&str>,
    pos_col: Option<&str>,
    p_col: Option<&str>,
    mlog10p_col: Option<&str>,
) -> Result<RawSumstats, String> {
    let normalized_sep = normalize_sep(sep);
    let delim = if normalized_sep.len() == 1 {
        normalized_sep.chars().next().unwrap()
    } else {
        return Err("Multi-character separators are not supported.".to_string());
    };

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(delim as u8)
        .has_headers(true)
        .flexible(true)
        .from_reader(r);

    let headers: Vec<String> = rdr
        .headers()
        .map_err(|e| e.to_string())?
        .iter()
        .map(|s| s.to_string())
        .collect();

    let (chrom_src, pos_src, p_src, mlog_src) =
        detect_sumstat_columns(&headers, chrom_col, pos_col, p_col, mlog10p_col)?;

    let chrom_i = headers.iter().position(|h| h == &chrom_src).unwrap();
    let pos_i = headers.iter().position(|h| h == &pos_src).unwrap();
    let stat_i = if let Some(ref p) = p_src {
        headers.iter().position(|h| h == p).unwrap()
    } else {
        headers.iter().position(|h| h == mlog_src.as_ref().unwrap()).unwrap()
    };

    let mode_p = p_src.is_some();
    let mut chrom = Vec::new();
    let mut pos = Vec::new();
    let mut p_vals = Vec::new();
    let mut mlog_vals = Vec::new();

    for rec in rdr.records() {
        let rec = rec.map_err(|e| e.to_string())?;
        chrom.push(rec.get(chrom_i).unwrap_or("").to_string());
        pos.push(rec.get(pos_i).unwrap_or("").to_string());
        let s = rec.get(stat_i).unwrap_or("").to_string();
        if mode_p {
            p_vals.push(s);
            mlog_vals.push(String::new());
        } else {
            p_vals.push(String::new());
            mlog_vals.push(s);
        }
    }

    Ok(RawSumstats {
        chrom_tokens: chrom,
        pos_tokens: pos,
        p_tokens: p_vals,
        mlog10p_tokens: mlog_vals,
        has_p_column: mode_p,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;
    use std::path::PathBuf;

    fn fixture(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/fixtures")
            .join(name)
    }

    #[test]
    fn detect_columns_default_fixture_headers() {
        let p = fixture("sumstats_small.tsv");
        let raw = load_sumstats(&p, "\t", None, None, None, None).unwrap();
        assert!(raw.has_p_column);
        assert_eq!(raw.chrom_tokens.len(), 7);
    }

    #[test]
    fn load_escaped_tab_sep() {
        let p = fixture("sumstats_small.tsv");
        let raw = load_sumstats(&p, r"\t", None, None, None, None).unwrap();
        assert_eq!(raw.chrom_tokens.len(), 7);
    }

    #[test]
    fn load_gzipped_fixture_matches_plain() {
        let p = fixture("sumstats_small.tsv");
        let plain = load_sumstats(&p, "\t", None, None, None, None).unwrap();

        let data = std::fs::read(&p).unwrap();
        let mut gz_path = std::env::temp_dir();
        gz_path.push(format!(
            "gwaspeek_gz_test_{}.tsv.gz",
            std::process::id()
        ));
        {
            let f = File::create(&gz_path).unwrap();
            let mut enc = GzEncoder::new(f, Compression::default());
            enc.write_all(&data).unwrap();
            enc.finish().unwrap();
        }
        let gz = load_sumstats(&gz_path, "\t", None, None, None, None).unwrap();
        let _ = std::fs::remove_file(&gz_path);

        assert_eq!(gz.chrom_tokens, plain.chrom_tokens);
        assert_eq!(gz.pos_tokens, plain.pos_tokens);
        assert_eq!(gz.p_tokens, plain.p_tokens);
        assert_eq!(gz.mlog10p_tokens, plain.mlog10p_tokens);
        assert_eq!(gz.has_p_column, plain.has_p_column);
    }
}
