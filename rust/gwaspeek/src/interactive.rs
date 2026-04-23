use std::collections::{BTreeMap, HashMap};
use std::io::{self, BufRead, IsTerminal, Read, Write};
use std::path::Path;
use std::time::Duration;

use flate2::read::GzDecoder;
use std::fs::File;
use std::io::BufReader;

use crossterm::event::{self, Event, KeyCode, KeyEvent, KeyEventKind, MouseEventKind};
use crossterm::terminal::{disable_raw_mode, enable_raw_mode, size};

use crate::io::repo_data_dir;
use crate::manhattan::{
    cumulative_to_chr_pos, density_legend, render_manhattan, sig_threshold_legend,
    viewport_chr_label,
};
use crate::plot_state::{normalize_build, prepare_plot_dataset, visible_mask, PlotDataset};
use crate::preprocess::CleanSumstats;

const PAN_FRAC_FINE: f64 = 0.10;
const PAN_FRAC_COARSE: f64 = 0.20;
const ZOOM_IN_FINE: f64 = 1.25;
const ZOOM_OUT_FINE: f64 = 1.0 / ZOOM_IN_FINE;
const ZOOM_IN_COARSE: f64 = 2.0;
const ZOOM_OUT_COARSE: f64 = 0.5;

pub struct Viewport {
    pub global_min: f64,
    pub global_max: f64,
    pub start: f64,
    pub end: f64,
    pub min_window_bp: f64,
}

impl Viewport {
    pub fn new(global_min: f64, global_max: f64) -> Self {
        Self {
            global_min,
            global_max,
            start: global_min,
            end: global_max,
            min_window_bp: 10_000.0,
        }
    }

    pub fn width(&self) -> f64 {
        self.end - self.start
    }

    fn clamp(&mut self) {
        if self.start < self.global_min {
            let shift = self.global_min - self.start;
            self.start += shift;
            self.end += shift;
        }
        if self.end > self.global_max {
            let shift = self.end - self.global_max;
            self.start -= shift;
            self.end -= shift;
        }
        self.start = self.start.max(self.global_min);
        self.end = self.end.min(self.global_max);
        if self.end <= self.start {
            self.start = self.global_min;
            self.end = self.global_max;
        }
    }

    pub fn zoom(&mut self, factor: f64, anchor_ratio: f64) {
        let anchor_ratio = anchor_ratio.clamp(0.0, 1.0);
        let anchor = self.start + anchor_ratio * self.width();
        let mut new_width = self.width() / factor;
        let min_width = self.min_window_bp.max(1.0);
        let full_width = self.global_max - self.global_min;
        new_width = new_width.clamp(min_width, full_width);
        self.start = anchor - anchor_ratio * new_width;
        self.end = self.start + new_width;
        self.clamp();
    }

    pub fn pan(&mut self, delta: f64) {
        self.start += delta;
        self.end += delta;
        self.clamp();
    }

    pub fn set_window(&mut self, start: f64, end: f64) {
        if end <= start {
            return;
        }
        self.start = start;
        self.end = end;
        self.clamp();
    }

    pub fn reset(&mut self) {
        self.start = self.global_min;
        self.end = self.global_max;
    }
}

#[derive(Clone, Debug)]
pub struct FrameSummary {
    pub region: String,
    pub view_size: String,
    pub n_vars: usize,
    pub n_chrs: usize,
    pub lead_label: Option<String>,
    pub lead_gene: Option<String>,
    pub gene_panel_active: bool,
}

pub fn default_gtf_path() -> std::path::PathBuf {
    repo_data_dir().join("GRCh37_latest_genomic.gene_only.gtf.gz")
}

pub fn default_gtf38_path() -> std::path::PathBuf {
    repo_data_dir().join("GRCh38_latest_genomic.gene_only.gtf.gz")
}

fn parse_region(region: &str) -> Result<(String, i32, i32), String> {
    let region = region.trim();
    let Some((left, right)) = region.split_once(':') else {
        return Err("Region must match chr:start-end (example: 3:100000-200000)".to_string());
    };
    let Some((s, e)) = right.split_once('-') else {
        return Err("Region must match chr:start-end (example: 3:100000-200000)".to_string());
    };
    let mut chrom_token = left.trim().to_uppercase();
    if chrom_token.starts_with("CHR") {
        chrom_token = chrom_token[3..].to_string();
    }
    let start: i32 = s.trim().parse().map_err(|_| "invalid start".to_string())?;
    let end: i32 = e.trim().parse().map_err(|_| "invalid end".to_string())?;
    if end <= start {
        return Err("Region end must be greater than region start".to_string());
    }
    if chrom_token == "M" {
        chrom_token = "MT".to_string();
    }
    Ok((chrom_token, start, end))
}

fn region_to_window(
    region: &str,
    offsets: &BTreeMap<i32, f64>,
    chr_sizes: &BTreeMap<i32, f64>,
) -> Result<(f64, f64), String> {
    let (chrom_token, start, end) = parse_region(region)?;
    let chrom = match chrom_token.as_str() {
        "X" => 23,
        "Y" => 24,
        "MT" => 25,
        _ => chrom_token
            .parse::<i32>()
            .map_err(|_| format!("Invalid chromosome: {chrom_token}"))?,
    };
    if !offsets.contains_key(&chrom) || !chr_sizes.contains_key(&chrom) {
        return Err(format!("Chromosome {chrom_token} not present in data"));
    }
    let max_pos = *chr_sizes.get(&chrom).unwrap() as i32;
    if start > max_pos {
        return Err(format!(
            "Region start exceeds chromosome max position ({max_pos})"
        ));
    }
    let clamped_end = end.min(max_pos);
    if clamped_end <= start {
        return Err("Region has no visible span after clamping to chromosome bounds".to_string());
    }
    let base = *offsets.get(&chrom).unwrap();
    Ok((base + start as f64, base + clamped_end as f64))
}

fn chr_token(chrom: i32) -> String {
    match chrom {
        23 => "X".to_string(),
        24 => "Y".to_string(),
        25 => "MT".to_string(),
        _ => chrom.to_string(),
    }
}

fn format_bp_span(span: f64) -> String {
    let bp = (span.round() as i64).max(1);
    if bp >= 1_000_000 {
        format!("{:.2}Mb", bp as f64 / 1_000_000.0)
    } else if bp >= 1_000 {
        format!("{:.1}kb", bp as f64 / 1_000.0)
    } else {
        format!("{bp}bp")
    }
}

fn nc_prefix_to_chrom() -> &'static [(&'static str, i32)] {
    &[
        ("NC_000001", 1),
        ("NC_000002", 2),
        ("NC_000003", 3),
        ("NC_000004", 4),
        ("NC_000005", 5),
        ("NC_000006", 6),
        ("NC_000007", 7),
        ("NC_000008", 8),
        ("NC_000009", 9),
        ("NC_000010", 10),
        ("NC_000011", 11),
        ("NC_000012", 12),
        ("NC_000013", 13),
        ("NC_000014", 14),
        ("NC_000015", 15),
        ("NC_000016", 16),
        ("NC_000017", 17),
        ("NC_000018", 18),
        ("NC_000019", 19),
        ("NC_000020", 20),
        ("NC_000021", 21),
        ("NC_000022", 22),
        ("NC_000023", 23),
        ("NC_000024", 24),
        ("NC_012920", 25),
    ]
}

fn parse_gtf_seqname_to_chrom(seqname: &str) -> Option<i32> {
    let raw = seqname.trim();
    if raw.is_empty() {
        return None;
    }
    let mut token = raw.to_uppercase();
    if token.starts_with("CHR") {
        token = token[3..].to_string();
    }
    if token == "M" {
        token = "MT".to_string();
    }
    match token.as_str() {
        "X" => return Some(23),
        "Y" => return Some(24),
        "MT" => return Some(25),
        _ => {}
    }
    if let Ok(iv) = token.parse::<i32>() {
        return Some(iv);
    }
    if token.starts_with("NC_") || token.starts_with("NT_") {
        let base = token.split('.').next().unwrap_or("");
        for (p, c) in nc_prefix_to_chrom() {
            if base == *p {
                return Some(*c);
            }
        }
    }
    None
}

fn parse_gtf_attrs(attr_field: &str) -> HashMap<String, String> {
    let mut out = HashMap::new();
    for chunk in attr_field.split(';') {
        let token = chunk.trim();
        if token.is_empty() || !token.contains(' ') {
            continue;
        }
        let mut sp = token.splitn(2, ' ');
        let key = sp.next().unwrap();
        let raw = sp.next().unwrap_or("");
        let val = raw.trim().trim_matches('"').to_string();
        out.insert(key.to_string(), val);
    }
    out
}

fn load_protein_coding_genes(gtf_path: &Path) -> HashMap<i32, Vec<(i32, i32, String)>> {
    let mut genes: HashMap<i32, Vec<(i32, i32, String)>> = HashMap::new();
    if !gtf_path.exists() {
        return genes;
    }
    let file = match File::open(gtf_path) {
        Ok(f) => f,
        Err(_) => return genes,
    };
    let reader: Box<dyn std::io::Read> = if gtf_path.to_string_lossy().ends_with(".gz") {
        Box::new(GzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let br = BufReader::new(reader);
    for line in br.lines().map_while(Result::ok) {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 9 || parts[2] != "gene" {
            continue;
        }
        let Some(chrom) = parse_gtf_seqname_to_chrom(parts[0]) else {
            continue;
        };
        let attrs = parse_gtf_attrs(parts[8]);
        let gene_type = attrs
            .get("gene_type")
            .or_else(|| attrs.get("gene_biotype"))
            .map(|s| s.as_str());
        if gene_type != Some("protein_coding") {
            continue;
        }
        let start: i32 = parts[3].parse().unwrap_or(0);
        let end: i32 = parts[4].parse().unwrap_or(0);
        let name = attrs
            .get("gene_name")
            .or_else(|| attrs.get("gene_id"))
            .cloned()
            .unwrap_or_else(|| "gene".to_string());
        genes.entry(chrom).or_default().push((start, end, name));
    }
    for v in genes.values_mut() {
        v.sort_by_key(|g| g.0);
    }
    genes
}

pub struct GeneTrackStore {
    pub gtf37_path: std::path::PathBuf,
    pub gtf38_path: std::path::PathBuf,
    cache: HashMap<String, HashMap<i32, Vec<(i32, i32, String)>>>,
}

impl GeneTrackStore {
    pub fn new(gtf37: std::path::PathBuf, gtf38: std::path::PathBuf) -> Self {
        Self {
            gtf37_path: gtf37,
            gtf38_path: gtf38,
            cache: HashMap::new(),
        }
    }

    fn path_for(&self, mode: &str) -> &std::path::Path {
        if mode == "37" {
            &self.gtf37_path
        } else {
            &self.gtf38_path
        }
    }

    pub fn exists(&self, mode: &str) -> bool {
        mode != "off" && self.path_for(mode).exists()
    }

    pub fn loaded(&self, mode: &str) -> bool {
        self.cache.contains_key(mode)
    }

    pub fn get(&mut self, mode: &str) -> HashMap<i32, Vec<(i32, i32, String)>> {
        if mode == "off" {
            return HashMap::new();
        }
        if !self.cache.contains_key(mode) {
            let g = load_protein_coding_genes(self.path_for(mode));
            self.cache.insert(mode.to_string(), g);
        }
        self.cache.get(mode).cloned().unwrap_or_default()
    }
}

fn lead_variant_in_view(data: &PlotDataset, mask: &[bool]) -> Option<(f64, f64, String)> {
    let mut best_i: Option<usize> = None;
    let mut best_y = f64::NEG_INFINITY;
    for (i, m) in mask.iter().enumerate() {
        if !*m {
            continue;
        }
        let y = data.mlog10p[i];
        if y > best_y {
            best_y = y;
            best_i = Some(i);
        }
    }
    let i = best_i?;
    let chrom = data.chrom[i];
    let pos = data.pos[i] as i64;
    Some((
        data.x[i],
        data.mlog10p[i],
        format!("lead {}:{}", chr_token(chrom), pos),
    ))
}

fn nearest_gene_label_for_x(
    x_cum: f64,
    offsets: &BTreeMap<i32, f64>,
    chr_sizes: &BTreeMap<i32, f64>,
    genes_by_chr: &HashMap<i32, Vec<(i32, i32, String)>>,
) -> Option<String> {
    let sorted: Vec<i32> = chr_sizes.keys().copied().collect();
    let (chrom, pos) = cumulative_to_chr_pos(x_cum, &sorted, offsets, chr_sizes);
    let genes = genes_by_chr.get(&chrom)?;
    if genes.is_empty() {
        return None;
    }
    let mut best_name: Option<String> = None;
    let mut best_dist: Option<i32> = None;
    for (g_start, g_end, g_name) in genes {
        if *g_start <= pos && pos <= *g_end {
            return Some(g_name.clone());
        }
        let dist = if pos < *g_start {
            g_start - pos
        } else {
            pos - g_end
        };
        if best_dist.is_none() || dist < best_dist.unwrap() {
            best_dist = Some(dist);
            best_name = Some(g_name.clone());
        }
    }
    best_name
}

fn nearest_gene_for_chr_pos(
    chrom: i32,
    pos: i32,
    genes_by_chr: &HashMap<i32, Vec<(i32, i32, String)>>,
) -> Option<String> {
    let genes = genes_by_chr.get(&chrom)?;
    let mut best_name: Option<String> = None;
    let mut best_dist: Option<i32> = None;
    for (g_start, g_end, g_name) in genes {
        if *g_start <= pos && pos <= *g_end {
            return Some(g_name.clone());
        }
        let dist = if pos < *g_start {
            g_start - pos
        } else {
            pos - g_end
        };
        if best_dist.is_none() || dist < best_dist.unwrap() {
            best_dist = Some(dist);
            best_name = Some(g_name.clone());
        }
    }
    best_name
}

fn protein_coding_track_in_view(
    x_start: f64,
    x_end: f64,
    offsets: &BTreeMap<i32, f64>,
    chr_sizes: &BTreeMap<i32, f64>,
    genes_by_chr: &HashMap<i32, Vec<(i32, i32, String)>>,
) -> Vec<(f64, f64, String)> {
    let sorted: Vec<i32> = chr_sizes.keys().copied().collect();
    let (c0, p0) = cumulative_to_chr_pos(x_start, &sorted, offsets, chr_sizes);
    let (c1, p1) = cumulative_to_chr_pos(x_end, &sorted, offsets, chr_sizes);

    fn genes_on_chrom(
        chrom: i32,
        lo_bp: i32,
        hi_bp: i32,
        offsets: &BTreeMap<i32, f64>,
        genes_by_chr: &HashMap<i32, Vec<(i32, i32, String)>>,
    ) -> Vec<(f64, f64, String)> {
        let span_bp = (hi_bp - lo_bp).abs();
        if span_bp > 1_000_000 {
            return vec![];
        }
        let Some(genes) = genes_by_chr.get(&chrom) else {
            return vec![];
        };
        if genes.is_empty() {
            return vec![];
        }
        let lo = lo_bp.min(hi_bp);
        let hi = lo_bp.max(hi_bp);
        let mut out = Vec::new();
        let off = *offsets.get(&chrom).unwrap();
        for (g_start, g_end, g_name) in genes {
            if *g_end < lo || *g_start > hi {
                continue;
            }
            out.push((off + *g_start as f64, off + *g_end as f64, g_name.clone()));
            if out.len() >= 25 {
                break;
            }
        }
        out
    }

    if c0 == c1 {
        return genes_on_chrom(c0, p0, p1, offsets, genes_by_chr);
    }

    let x_lo = x_start.min(x_end);
    let x_hi = x_start.max(x_end);
    let x_mid = 0.5 * (x_lo + x_hi);
    let (c_mid, _) = cumulative_to_chr_pos(x_mid, &sorted, offsets, chr_sizes);
    let chr_lo = *offsets.get(&c_mid).unwrap();
    let chr_hi = chr_lo + *chr_sizes.get(&c_mid).unwrap();
    let vx0 = x_lo.max(chr_lo);
    let vx1 = x_hi.min(chr_hi);
    if vx1 <= vx0 {
        return vec![];
    }
    let mut p_lo = (vx0 - chr_lo).round() as i32;
    let mut p_hi = (vx1 - chr_lo).round() as i32;
    p_lo = p_lo.max(1);
    p_hi = p_hi.min(*chr_sizes.get(&c_mid).unwrap() as i32);
    genes_on_chrom(c_mid, p_lo, p_hi, offsets, genes_by_chr)
}

fn inspect_view(
    data: &PlotDataset,
    viewport: &Viewport,
    show_lead: bool,
    show_track: bool,
    genes_by_chr: &HashMap<i32, Vec<(i32, i32, String)>>,
    y_min: f64,
    non_human: bool,
) -> (
    Vec<bool>,
    Option<(f64, f64, String)>,
    Vec<(f64, f64, String)>,
    FrameSummary,
    String,
) {
    let show_track = if non_human { false } else { show_track };
    let mask = visible_mask(data, viewport.start, viewport.end);
    let mut n_vars = 0usize;
    let mut chrs = std::collections::HashSet::new();
    for (i, m) in mask.iter().enumerate() {
        if *m {
            n_vars += 1;
            chrs.insert(data.chrom[i]);
        }
    }
    let n_chrs = chrs.len();
    let lead_variant = if show_lead {
        lead_variant_in_view(data, &mask)
    } else {
        None
    };
    let mut lead_gene = None;
    if let Some(ref lv) = lead_variant {
        if !genes_by_chr.is_empty() {
            lead_gene = nearest_gene_label_for_x(
                lv.0,
                &data.layout.offsets,
                &data.layout.chr_sizes,
                genes_by_chr,
            );
        }
    }
    let track_allowed = show_track && viewport.width() <= 1_000_000.0;
    let gene_track = if track_allowed && !genes_by_chr.is_empty() {
        protein_coding_track_in_view(
            viewport.start,
            viewport.end,
            &data.layout.offsets,
            &data.layout.chr_sizes,
            genes_by_chr,
        )
    } else {
        vec![]
    };

    let region = viewport_chr_label(
        viewport.start,
        viewport.end,
        &data.layout.offsets,
        &data.layout.chr_sizes,
    );
    let view_size = format_bp_span(viewport.width());
    let mut title = format!(
        "Manhattan {region} | view {view_size} | vars {n_vars} | skip>={y_min:.1}"
    );
    if n_chrs > 1 {
        title = format!("{title} | chrs {n_chrs}");
    }
    if let Some(ref lv) = lead_variant {
        if let Some(ref g) = lead_gene {
            title = format!("{title} [{} | gene {g}]", lv.2);
        } else {
            title = format!("{title} [{}]", lv.2);
        }
    }

    let summary = FrameSummary {
        region: region.clone(),
        view_size,
        n_vars,
        n_chrs,
        lead_label: lead_variant.as_ref().map(|x| x.2.clone()),
        lead_gene,
        gene_panel_active: track_allowed,
    };
    (mask, lead_variant, gene_track, summary, title)
}

fn render_frame(
    data: &PlotDataset,
    viewport: &Viewport,
    show_lead: bool,
    show_track: bool,
    genes_by_chr: &HashMap<i32, Vec<(i32, i32, String)>>,
    width: usize,
    height: usize,
    sig_level: f64,
    ymax: Option<f64>,
    unicode: bool,
    y_min: f64,
    color: bool,
    non_human: bool,
    light_theme: bool,
) -> Result<(String, FrameSummary), String> {
    let (mask, lead_variant, gene_track, summary, title) = inspect_view(
        data,
        viewport,
        show_lead,
        show_track,
        genes_by_chr,
        y_min,
        non_human,
    );
    let frame = render_manhattan(
        data,
        width,
        height,
        sig_level,
        ymax,
        Some(viewport.start),
        Some(viewport.end),
        Some(title.as_str()),
        unicode,
        y_min,
        lead_variant,
        if gene_track.is_empty() {
            None
        } else {
            Some(gene_track.as_slice())
        },
        summary.gene_panel_active,
        Some(mask.as_slice()),
        color,
        light_theme,
    )?;
    Ok((frame, summary))
}

fn clear_screen(out: &mut impl Write) -> io::Result<()> {
    out.write_all(b"\x1b[H\x1b[2J")?;
    Ok(())
}

fn enter_alt_screen(out: &mut impl Write) -> io::Result<()> {
    out.write_all(b"\x1b[?1049h\x1b[?1000h\x1b[?1006h\x1b[?25l")?;
    Ok(())
}

fn leave_alt_screen(out: &mut impl Write) -> io::Result<()> {
    out.write_all(b"\x1b[?25h\x1b[?1006l\x1b[?1000l\x1b[?1049l")?;
    Ok(())
}

fn plot_anchor_ratio(col_1based: u16, row_1based: u16, frame_width: usize, frame_height: usize) -> Option<f64> {
    let row = row_1based as usize;
    if row < 1 || row > frame_height {
        return None;
    }
    let x_plot0 = 7usize;
    let w = frame_width.saturating_sub(x_plot0 + 1);
    if w == 0 {
        return None;
    }
    let x = col_1based.saturating_sub(1) as usize;
    if x < x_plot0 || x > x_plot0 + w {
        return None;
    }
    Some((x - x_plot0) as f64 / w as f64)
}

fn help_text() -> &'static str {
    "Keys: q quit  h help  v vars  t 37|38|off  g jump  l lead  m theme  r reset  a/d pan (A/D fast)  w/s zoom (W/S fast)  wheel@cursor  +/-"
}

fn render_help_screen() -> String {
    let lines = [
        "gwaspeek interactive help",
        "",
        "Navigation:",
        "  a/d          : pan (fine)",
        "  A/D          : pan (coarse)",
        "  w/s          : zoom out/in (fine)",
        "  W/S          : zoom out/in (coarse)",
        "  +/-          : zoom in/out (fine)",
        "  mouse wheel  : zoom at cursor position",
        "  g            : jump to region (chr:start-end)",
        "",
        "Toggles:",
        "  l : toggle lead variant annotation",
        "  t : cycle gene track mode (37 -> 38 -> off)",
        "  v : toggle variants-in-view list",
        "  h : toggle this help",
        "  m : toggle dark / light palette (plot + status colors)",
        "  r : reset to full genome view",
        "  q : quit",
        "",
        "Dense cells use heavier glyphs instead of dropping overlapping points.",
        "Press h (or v) again to return to plot.",
    ];
    lines.join("\n")
}

fn render_variants_view(
    data: &PlotDataset,
    _viewport: &Viewport,
    mask: &[bool],
    summary: &FrameSummary,
    genes_by_chr: &HashMap<i32, Vec<(i32, i32, String)>>,
    y_min: f64,
) -> String {
    let mut idx: Vec<usize> = mask
        .iter()
        .enumerate()
        .filter(|(_, m)| **m)
        .map(|(i, _)| i)
        .collect();
    if idx.is_empty() {
        return format!(
            "Variants in view: none\nRegion: {}\nSkip threshold: -log10P >= {y_min:.2}",
            summary.region
        );
    }
    idx.sort_by(|&a, &b| data.mlog10p[b].partial_cmp(&data.mlog10p[a]).unwrap());
    idx.truncate(25);
    let mut lines = vec![
        format!(
            "Variants in view ({} shown, sorted by -log10P):",
            idx.len()
        ),
        format!("Region: {}", summary.region),
        format!("Skip threshold: -log10P >= {y_min:.2}"),
        "CHR\tPOS\tP\t-log10P\tnearest_gene".to_string(),
    ];
    for i in idx {
        let chrom = data.chrom[i];
        let pos = data.pos[i] as i64;
        let nearest = nearest_gene_for_chr_pos(chrom, pos as i32, genes_by_chr)
            .unwrap_or_else(|| "-".to_string());
        lines.push(format!(
            "{chrom}\t{pos}\t{p:.3e}\t{m:.3}\t{nearest}",
            p = data.p[i],
            m = data.mlog10p[i]
        ));
    }
    lines.join("\n")
}

fn genes_for_view(
    store: &mut GeneTrackStore,
    track_mode: &str,
    viewport: &Viewport,
    view_mode: &str,
    non_human: bool,
) -> HashMap<i32, Vec<(i32, i32, String)>> {
    if non_human || track_mode == "off" {
        return HashMap::new();
    }
    if view_mode == "variants" || viewport.width() <= 1_000_000.0 {
        return store.get(track_mode);
    }
    HashMap::new()
}

fn ansi_enabled(enabled: bool) -> bool {
    let term = std::env::var("TERM").unwrap_or_default();
    io::stdout().is_terminal()
        && term != "" && term != "dumb"
        && std::env::var("NO_COLOR").is_err()
        && enabled
}

fn ansi(text: &str, code: &str, enabled: bool) -> String {
    if !enabled {
        return text.to_string();
    }
    format!("\x1b[{code}m{text}\x1b[0m")
}

fn truncate_line(text: &str, width: usize) -> String {
    if width == 0 {
        return String::new();
    }
    let chs: Vec<char> = text.chars().collect();
    if chs.len() <= width {
        return text.to_string();
    }
    if width <= 3 {
        return chs.into_iter().take(width).collect();
    }
    let head: String = chs.into_iter().take(width - 3).collect::<String>().trim_end().to_string();
    format!("{head}...")
}

fn format_track_state(
    track_mode: &str,
    store: &GeneTrackStore,
    summary: &FrameSummary,
    non_human: bool,
) -> String {
    if non_human {
        return "track off (non-human)".to_string();
    }
    if track_mode == "off" {
        return "track off".to_string();
    }
    if !store.exists(track_mode) {
        return format!("track {track_mode} missing");
    }
    if summary.gene_panel_active {
        let state = if store.loaded(track_mode) {
            "loaded"
        } else {
            "loading"
        };
        format!("track {track_mode} {state}")
    } else {
        format!("track {track_mode} ready (zoom <=1Mb)")
    }
}

fn status_line(
    summary: &FrameSummary,
    build: &str,
    track_mode: &str,
    store: &GeneTrackStore,
    notice: Option<&str>,
    width: usize,
    color: bool,
    light_theme: bool,
    non_human: bool,
) -> String {
    let mut parts = vec![
        format!("build {build}"),
        format!("region {}", summary.region),
        format!("view {}", summary.view_size),
        format!("vars {}", summary.n_vars),
        format!("chrs {}", summary.n_chrs),
        format_track_state(track_mode, store, summary, non_human),
    ];
    if let Some(ll) = &summary.lead_label {
        let mut lead = ll.clone();
        if let Some(g) = &summary.lead_gene {
            lead = format!("{lead} ({g})");
        }
        parts.push(lead);
    }
    if let Some(n) = notice {
        parts.push(n.to_string());
    }
    let line = truncate_line(&parts.join(" | "), width);
    let sgr = if notice.is_some() {
        if light_theme {
            "31"
        } else {
            "33"
        }
    } else if light_theme {
        "34"
    } else {
        "36"
    };
    ansi(&line, sgr, color)
}

fn footer_help_line(
    width: usize,
    unicode: bool,
    color: bool,
    sig_level: f64,
    light_theme: bool,
) -> String {
    let line1 = truncate_line(help_text(), width);
    let line2 = truncate_line(
        &format!(
            "{}  |  {}",
            density_legend(unicode),
            sig_threshold_legend(unicode, sig_level)
        ),
        width,
    );
    let text = format!("{line1}\n{line2}");
    let muted = if light_theme { "30" } else { "37" };
    ansi(&text, muted, color)
}

fn next_track_mode(track_mode: &str) -> &'static str {
    match track_mode {
        "37" => "38",
        "38" => "off",
        _ => "37",
    }
}

fn active_build(track_mode: &str, base_build: &str) -> String {
    if track_mode == "37" || track_mode == "38" {
        track_mode.to_string()
    } else {
        base_build.to_string()
    }
}

fn remap_cumulative_position(x_cum: f64, old_data: &PlotDataset, new_data: &PlotDataset) -> f64 {
    let sorted_old: Vec<i32> = old_data.layout.chr_sizes.keys().copied().collect();
    let (chrom, pos) = cumulative_to_chr_pos(
        x_cum,
        &sorted_old,
        &old_data.layout.offsets,
        &old_data.layout.chr_sizes,
    );
    if !new_data.layout.offsets.contains_key(&chrom) {
        return new_data.layout.x_min;
    }
    let max_pos = *new_data.layout.chr_sizes.get(&chrom).unwrap() as i32;
    let pos = pos.max(1).min(max_pos);
    new_data.layout.offsets.get(&chrom).unwrap() + pos as f64
}

fn remap_viewport_to_build(viewport: &mut Viewport, old_data: &PlotDataset, new_data: &PlotDataset) {
    let new_start = remap_cumulative_position(viewport.start, old_data, new_data);
    let new_end = remap_cumulative_position(viewport.end, old_data, new_data);
    viewport.global_min = new_data.layout.x_min;
    viewport.global_max = new_data.layout.x_max;
    viewport.start = new_start;
    viewport.end = new_end;
    viewport.clamp();
}

fn apply_key(
    key: char,
    viewport: &mut Viewport,
    show_lead: &mut bool,
    track_mode: &mut String,
    non_human: bool,
) -> bool {
    match key {
        'q' | 'Q' => return false,
        'r' | 'R' => viewport.reset(),
        'l' | 'L' => *show_lead = !*show_lead,
        't' | 'T' => {
            if non_human {
                *track_mode = "off".to_string();
            } else {
                *track_mode = next_track_mode(track_mode.as_str()).to_string();
            }
        }
        'a' => viewport.pan(-PAN_FRAC_FINE * viewport.width()),
        'A' => viewport.pan(-PAN_FRAC_COARSE * viewport.width()),
        'd' => viewport.pan(PAN_FRAC_FINE * viewport.width()),
        'D' => viewport.pan(PAN_FRAC_COARSE * viewport.width()),
        'w' | '-' | '_' => viewport.zoom(ZOOM_OUT_FINE, 0.5),
        'W' => viewport.zoom(ZOOM_OUT_COARSE, 0.5),
        's' | '+' | '=' => viewport.zoom(ZOOM_IN_FINE, 0.5),
        'S' => viewport.zoom(ZOOM_IN_COARSE, 0.5),
        _ => {}
    }
    true
}

pub fn run_interactive_manhattan(
    clean: CleanSumstats,
    width: usize,
    height: usize,
    sig_level: f64,
    ymax: Option<f64>,
    unicode: bool,
    y_min: f64,
    gtf_path: Option<&Path>,
    gtf38_path: Option<&Path>,
    build: &str,
    color: bool,
    non_human: bool,
) -> Result<(), String> {
    let gtf37 = gtf_path
        .map(|p| p.to_path_buf())
        .unwrap_or_else(default_gtf_path);
    let gtf38 = gtf38_path
        .map(|p| p.to_path_buf())
        .unwrap_or_else(default_gtf38_path);

    let base_build = normalize_build(Some(build));
    let ds37 = prepare_plot_dataset(&clean, Some("37"), non_human)?;
    let ds38 = prepare_plot_dataset(&clean, Some("38"), non_human)?;
    let mut datasets: HashMap<String, PlotDataset> = HashMap::new();
    datasets.insert("37".to_string(), ds37);
    datasets.insert("38".to_string(), ds38);

    let mut show_lead = true;
    let mut track_mode = if non_human {
        "off".to_string()
    } else {
        base_build.clone()
    };
    let mut view_mode = "plot".to_string();
    let mut light_theme = false;
    let mut notice: Option<String> = None;

    let data0 = datasets
        .get(&active_build(&track_mode, &base_build))
        .unwrap();
    let mut viewport = Viewport::new(data0.layout.x_min, data0.layout.x_max);
    let mut gene_store = GeneTrackStore::new(gtf37, gtf38);

    if !io::stdin().is_terminal() {
        loop {
            let active_build_s = active_build(&track_mode, &base_build);
            let data = datasets.get(&active_build_s).unwrap();
            let genes_by_chr = genes_for_view(
                &mut gene_store,
                &track_mode,
                &viewport,
                &view_mode,
                non_human,
            );
            let (mask, _, _, summary, _) = inspect_view(
                data,
                &viewport,
                show_lead,
                track_mode != "off",
                &genes_by_chr,
                y_min,
                non_human,
            );
            clear_screen(&mut io::stdout()).map_err(|e| e.to_string())?;
            if view_mode == "help" {
                println!("{}", render_help_screen());
            } else if view_mode == "variants" {
                println!("{}", render_variants_view(data, &viewport, &mask, &summary, &genes_by_chr, y_min));
            } else {
                let (frame, _) = render_frame(
                    data,
                    &viewport,
                    show_lead,
                    track_mode != "off",
                    &genes_by_chr,
                    width,
                    height,
                    sig_level,
                    ymax,
                    unicode,
                    y_min,
                    false,
                    non_human,
                    false,
                )?;
                println!("{frame}");
            }
            println!(
                "{}",
                status_line(
                    &summary,
                    &active_build_s,
                    &track_mode,
                    &gene_store,
                    notice.as_deref(),
                    width,
                    false,
                    false,
                    non_human,
                )
            );
            println!(
                "{}",
                footer_help_line(width, unicode, false, sig_level, false)
            );
            let mut buf = [0u8; 1];
            let n = io::stdin().read(&mut buf).map_err(|e| e.to_string())?;
            notice = None;
            if n == 0 {
                break;
            }
            let key = buf[0] as char;
            if key == 'h' || key == 'H' {
                view_mode = if view_mode == "help" {
                    "plot".to_string()
                } else {
                    "help".to_string()
                };
                continue;
            }
            if key == 'v' || key == 'V' {
                view_mode = if view_mode == "variants" {
                    "plot".to_string()
                } else {
                    "variants".to_string()
                };
                continue;
            }
            let old_build = active_build_s.clone();
            if !apply_key(key, &mut viewport, &mut show_lead, &mut track_mode, non_human) {
                break;
            }
            let new_build = active_build(&track_mode, &base_build);
            if new_build != old_build {
                remap_viewport_to_build(
                    &mut viewport,
                    datasets.get(&old_build).unwrap(),
                    datasets.get(&new_build).unwrap(),
                );
            }
        }
        return Ok(());
    }

    let mut stdout = io::stdout();
    let color_enabled = ansi_enabled(color);
    enable_raw_mode().map_err(|e| e.to_string())?;
    enter_alt_screen(&mut stdout).map_err(|e| e.to_string())?;
    let mut painted = false;
    let mut last_term: Option<(u16, u16)> = None;
    let mut last_ui_snap: Option<(String, bool, String, bool, f64, f64, f64)> = None;
    let mut quit = false;

    while !quit {
        let term_size = size().unwrap_or((width as u16, (height + 3) as u16));
        let frame_width = (term_size.0 as usize).max(20);
        let term_lines = term_size.1 as usize;
        let frame_height = (term_lines.saturating_sub(3)).max(8);
        let term_sig = (term_size.0, term_size.1);
        let resized = last_term.map(|t| t != term_sig).unwrap_or(false);
        last_term = Some(term_sig);

        let has_input = event::poll(Duration::from_millis(150)).map_err(|e| e.to_string())?;

        let mut key_char: Option<char> = None;
        if has_input {
            match event::read().map_err(|e| e.to_string())? {
                Event::Key(KeyEvent {
                    code,
                    kind: KeyEventKind::Press,
                    ..
                }) => {
                    if let KeyCode::Char(c) = code {
                        key_char = Some(c);
                    }
                }
                Event::Mouse(me) => {
                    let col = me.column.saturating_add(1);
                    let row = me.row.saturating_add(1);
                    match me.kind {
                        MouseEventKind::ScrollUp => {
                            if let Some(ar) =
                                plot_anchor_ratio(col, row, frame_width, frame_height)
                            {
                                viewport.zoom(ZOOM_IN_FINE, ar);
                            }
                        }
                        MouseEventKind::ScrollDown => {
                            if let Some(ar) =
                                plot_anchor_ratio(col, row, frame_width, frame_height)
                            {
                                viewport.zoom(ZOOM_OUT_FINE, ar);
                            }
                        }
                        _ => {}
                    }
                }
                _ => {}
            }
        }

        let active_build_s = active_build(&track_mode, &base_build);

        if let Some(k) = key_char {
            let data = datasets.get(&active_build_s).unwrap();
            let offsets_pre = &data.layout.offsets;
            let chr_sizes_pre = &data.layout.chr_sizes;

            if k == 'h' || k == 'H' {
                view_mode = if view_mode == "help" {
                    "plot".to_string()
                } else {
                    "help".to_string()
                };
            } else if k == 'v' || k == 'V' {
                view_mode = if view_mode == "variants" {
                    "plot".to_string()
                } else {
                    "variants".to_string()
                };
            } else if k == 'g' || k == 'G' {
                view_mode = "plot".to_string();
                disable_raw_mode().map_err(|e| e.to_string())?;
                leave_alt_screen(&mut stdout).map_err(|e| e.to_string())?;
                stdout
                    .write_all(b"\x1b[?25h\x1b[H\x1b[2J")
                    .map_err(|e| e.to_string())?;
                stdout
                    .write_all(b"Jump to region (chr:start-end, blank cancels): ")
                    .map_err(|e| e.to_string())?;
                stdout.flush().map_err(|e| e.to_string())?;
                let mut line = String::new();
                io::stdin()
                    .read_line(&mut line)
                    .map_err(|e| e.to_string())?;
                let region = line.trim();
                notice = if region.is_empty() {
                    Some("Jump cancelled".to_string())
                } else {
                    match region_to_window(region, offsets_pre, chr_sizes_pre) {
                        Ok((s, e)) => {
                            viewport.set_window(s, e);
                            Some(format!("Jumped to {region}"))
                        }
                        Err(e) => Some(format!("Jump failed: {e}")),
                    }
                };
                enable_raw_mode().map_err(|e| e.to_string())?;
                enter_alt_screen(&mut stdout).map_err(|e| e.to_string())?;
            } else if k == 'm' {
                light_theme = !light_theme;
            } else {
                let old_build = active_build_s.clone();
                if !apply_key(k, &mut viewport, &mut show_lead, &mut track_mode, non_human) {
                    quit = true;
                } else {
                    let new_build = active_build(&track_mode, &base_build);
                    if new_build != old_build {
                        remap_viewport_to_build(
                            &mut viewport,
                            datasets.get(&old_build).unwrap(),
                            datasets.get(&new_build).unwrap(),
                        );
                    }
                }
            }
        }

        if quit {
            break;
        }

        let active_build_s = active_build(&track_mode, &base_build);
        let data = datasets.get(&active_build_s).unwrap();

        let ui_snap = (
            view_mode.clone(),
            show_lead,
            track_mode.clone(),
            light_theme,
            viewport.start,
            viewport.end,
            viewport.width(),
        );
        if painted && !has_input && !resized && last_ui_snap.as_ref() == Some(&ui_snap) {
            continue;
        }

        let genes_by_chr = genes_for_view(
            &mut gene_store,
            &track_mode,
            &viewport,
            &view_mode,
            non_human,
        );
        let (mask, _, _, summary, _) = inspect_view(
            data,
            &viewport,
            show_lead,
            track_mode != "off",
            &genes_by_chr,
            y_min,
            non_human,
        );

        clear_screen(&mut stdout).map_err(|e| e.to_string())?;
        let body = if view_mode == "help" {
            render_help_screen()
        } else if view_mode == "variants" {
            render_variants_view(data, &viewport, &mask, &summary, &genes_by_chr, y_min)
        } else {
            let (frame, _) = render_frame(
                data,
                &viewport,
                show_lead,
                track_mode != "off",
                &genes_by_chr,
                frame_width,
                frame_height,
                sig_level,
                ymax,
                unicode,
                y_min,
                color_enabled,
                non_human,
                light_theme,
            )?;
            frame
        };
        write!(stdout, "{body}").map_err(|e| e.to_string())?;
        if !body.ends_with('\n') {
            writeln!(stdout).map_err(|e| e.to_string())?;
        }
        writeln!(
            stdout,
            "{}",
            status_line(
                &summary,
                &active_build_s,
                &track_mode,
                &gene_store,
                notice.as_deref(),
                frame_width,
                color_enabled,
                light_theme,
                non_human,
            )
        )
        .map_err(|e| e.to_string())?;
        writeln!(
            stdout,
            "{}",
            footer_help_line(
                frame_width,
                unicode,
                color_enabled,
                sig_level,
                light_theme
            )
        )
        .map_err(|e| e.to_string())?;
        stdout.flush().map_err(|e| e.to_string())?;

        painted = true;
        notice = None;
        last_ui_snap = Some(ui_snap);
    }

    let _ = disable_raw_mode();
    let _ = leave_alt_screen(&mut stdout);
    let _ = stdout.flush();

    Ok(())
}
