use std::collections::{BTreeMap, HashMap};
use std::sync::OnceLock;

use crate::preprocess::CleanSumstats;

pub const DEFAULT_BUILD: &str = "37";

fn canonical_chr_lengths(build: &str) -> &'static HashMap<i32, f64> {
    static M37: OnceLock<HashMap<i32, f64>> = OnceLock::new();
    static M38: OnceLock<HashMap<i32, f64>> = OnceLock::new();
    match build {
        "38" => M38.get_or_init(|| {
            let mut m = HashMap::new();
            m.insert(1, 248_956_422.0);
            m.insert(2, 242_193_529.0);
            m.insert(3, 198_295_559.0);
            m.insert(4, 190_214_555.0);
            m.insert(5, 181_538_259.0);
            m.insert(6, 170_805_979.0);
            m.insert(7, 159_345_973.0);
            m.insert(8, 145_138_636.0);
            m.insert(9, 138_394_717.0);
            m.insert(10, 133_797_422.0);
            m.insert(11, 135_086_622.0);
            m.insert(12, 133_275_309.0);
            m.insert(13, 114_364_328.0);
            m.insert(14, 107_043_718.0);
            m.insert(15, 101_991_189.0);
            m.insert(16, 90_338_345.0);
            m.insert(17, 83_257_441.0);
            m.insert(18, 80_373_285.0);
            m.insert(19, 58_617_616.0);
            m.insert(20, 64_444_167.0);
            m.insert(21, 46_709_983.0);
            m.insert(22, 50_818_468.0);
            m.insert(23, 156_040_895.0);
            m.insert(24, 57_227_415.0);
            m.insert(25, 16_569.0);
            m
        }),
        _ => M37.get_or_init(|| {
            let mut m = HashMap::new();
            m.insert(1, 249_250_621.0);
            m.insert(2, 243_199_373.0);
            m.insert(3, 198_022_430.0);
            m.insert(4, 191_154_276.0);
            m.insert(5, 180_915_260.0);
            m.insert(6, 171_115_067.0);
            m.insert(7, 159_138_663.0);
            m.insert(8, 146_364_022.0);
            m.insert(9, 141_213_431.0);
            m.insert(10, 135_534_747.0);
            m.insert(11, 135_006_516.0);
            m.insert(12, 133_851_895.0);
            m.insert(13, 115_169_878.0);
            m.insert(14, 107_349_540.0);
            m.insert(15, 102_531_392.0);
            m.insert(16, 90_354_753.0);
            m.insert(17, 81_195_210.0);
            m.insert(18, 78_077_248.0);
            m.insert(19, 59_128_983.0);
            m.insert(20, 63_025_520.0);
            m.insert(21, 48_129_895.0);
            m.insert(22, 51_304_566.0);
            m.insert(23, 155_270_560.0);
            m.insert(24, 59_373_566.0);
            m.insert(25, 16_569.0);
            m
        }),
    }
}

/// Exposed for tests mirroring `tests/test_preprocess.py`.
pub fn canonical_length(build: &str, chrom: i32) -> Option<f64> {
    canonical_chr_lengths(build).get(&chrom).copied()
}

#[derive(Debug, Clone)]
pub struct GenomeLayout {
    pub build: String,
    pub offsets: BTreeMap<i32, f64>,
    pub chr_sizes: BTreeMap<i32, f64>,
    pub x_min: f64,
    pub x_max: f64,
}

#[derive(Debug, Clone)]
pub struct PlotDataset {
    pub layout: GenomeLayout,
    pub x: Vec<f64>,
    pub chrom: Vec<i32>,
    pub pos: Vec<f64>,
    pub p: Vec<f64>,
    pub mlog10p: Vec<f64>,
}

pub fn normalize_build(build: Option<&str>) -> String {
    let token = build.unwrap_or(DEFAULT_BUILD).trim().to_uppercase();
    match token.as_str() {
        "37" | "GRCH37" | "HG19" => "37".to_string(),
        "38" | "GRCH38" | "HG38" => "38".to_string(),
        _ => DEFAULT_BUILD.to_string(),
    }
}

pub fn build_genome_layout(
    clean: &CleanSumstats,
    build: Option<&str>,
    data_driven_lengths: bool,
) -> Result<GenomeLayout, String> {
    let resolved_build = normalize_build(build);
    let chr_max_pos = &clean.chr_max_pos;
    if chr_max_pos.is_empty() {
        return Err("No variants available to build genome layout.".to_string());
    }

    let canonical = canonical_chr_lengths(resolved_build.as_str());

    let mut offsets: BTreeMap<i32, f64> = BTreeMap::new();
    let mut chr_sizes: BTreeMap<i32, f64> = BTreeMap::new();
    let mut offset = 0.0_f64;
    for (&chrom, &max_pos) in chr_max_pos.iter() {
        offsets.insert(chrom, offset);
        let span = if data_driven_lengths {
            max_pos.max(1.0)
        } else {
            *canonical.get(&chrom).unwrap_or(&max_pos)
        };
        chr_sizes.insert(chrom, span);
        offset += span;
    }

    let sorted_chroms: Vec<i32> = chr_sizes.keys().copied().collect();
    let first_chrom = *sorted_chroms.first().unwrap();
    let last_chrom = *sorted_chroms.last().unwrap();
    let x_min = *offsets.get(&first_chrom).unwrap();
    let x_max = offsets.get(&last_chrom).unwrap() + chr_sizes.get(&last_chrom).unwrap();

    Ok(GenomeLayout {
        build: resolved_build,
        offsets,
        chr_sizes,
        x_min,
        x_max,
    })
}

pub fn prepare_plot_dataset(
    clean: &CleanSumstats,
    build: Option<&str>,
    data_driven_lengths: bool,
) -> Result<PlotDataset, String> {
    let layout = build_genome_layout(clean, build, data_driven_lengths)?;
    let n = clean.chr.len();
    if n == 0 {
        return Ok(PlotDataset {
            layout,
            x: vec![],
            chrom: vec![],
            pos: vec![],
            p: vec![],
            mlog10p: vec![],
        });
    }
    let mut x = Vec::with_capacity(n);
    for i in 0..n {
        let c = clean.chr[i];
        let off = *layout.offsets.get(&c).unwrap();
        x.push(clean.pos[i] + off);
    }
    Ok(PlotDataset {
        layout,
        x,
        chrom: clean.chr.clone(),
        pos: clean.pos.clone(),
        p: clean.p.clone(),
        mlog10p: clean.mlog10p.clone(),
    })
}

pub fn visible_mask(dataset: &PlotDataset, x_start: f64, x_end: f64) -> Vec<bool> {
    if dataset.x.is_empty() {
        return vec![];
    }
    dataset
        .x
        .iter()
        .map(|xv| *xv >= x_start && *xv <= x_end)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::load_sumstats;
    use crate::preprocess::preprocess_sumstats;
    use std::path::PathBuf;

    fn fixture(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/fixtures")
            .join(name)
    }

    #[test]
    fn data_driven_uses_observed_max() {
        let raw = crate::preprocess::RawSumstats {
            chrom_tokens: vec!["23".into(), "23".into()],
            pos_tokens: vec!["5000000".into(), "8000000".into()],
            p_tokens: vec!["1e-8".into(), "1e-8".into()],
            mlog10p_tokens: vec![],
            has_p_column: true,
        };
        let clean = preprocess_sumstats(raw, 5.0).unwrap();
        let ds = prepare_plot_dataset(&clean, Some("37"), true).unwrap();
        assert_eq!(ds.layout.chr_sizes.get(&23).copied().unwrap(), 8_000_000.0);
        assert!(ds.layout.chr_sizes[&23] < canonical_length("37", 23).unwrap());
    }

    #[test]
    fn canonical_lengths_match_python_table() {
        let raw = load_sumstats(
            &fixture("sumstats_small.tsv"),
            "\t",
            None,
            None,
            None,
            None,
        )
        .unwrap();
        let clean = preprocess_sumstats(raw, 3.0).unwrap();
        let ds = prepare_plot_dataset(&clean, Some("37"), false).unwrap();
        assert_eq!(
            ds.layout.chr_sizes.get(&1).copied().unwrap(),
            canonical_length("37", 1).unwrap()
        );
    }
}
