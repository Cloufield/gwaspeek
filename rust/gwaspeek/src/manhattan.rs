use std::collections::{BTreeMap, HashMap};

use crate::plot_state::PlotDataset;
use crate::terminal_canvas::{CanvasStyle, TerminalCanvas};
use crate::versioning::package_version;

pub fn stdout_color_supported() -> bool {
    use std::io::IsTerminal;
    let term = std::env::var("TERM").unwrap_or_default();
    std::io::stdout().is_terminal()
        && term != "" && term != "dumb"
        && std::env::var("NO_COLOR").is_err()
}

fn chr_token(chrom: i32) -> String {
    match chrom {
        23 => "X".to_string(),
        24 => "Y".to_string(),
        25 => "MT".to_string(),
        _ => chrom.to_string(),
    }
}

fn chr_alternating_sgr(chrom: i32, chr_order: &HashMap<i32, usize>, light_theme: bool) -> &'static str {
    let Some(idx) = chr_order.get(&chrom) else {
        return "39";
    };
    if light_theme {
        if idx % 2 == 0 {
            "34"
        } else {
            "35"
        }
    } else if idx % 2 == 0 {
        "36"
    } else {
        "37"
    }
}

fn ansi_paint_glyph(
    glyph: char,
    chrom: i32,
    chr_order: &HashMap<i32, usize>,
    color: bool,
    light_theme: bool,
) -> String {
    if !color {
        return glyph.to_string();
    }
    let code = chr_alternating_sgr(chrom, chr_order, light_theme);
    format!("\x1b[{code}m{glyph}\x1b[0m")
}

pub fn cumulative_to_chr_pos(
    x_cum: f64,
    sorted_chroms: &[i32],
    offsets: &BTreeMap<i32, f64>,
    chr_sizes: &BTreeMap<i32, f64>,
) -> (i32, i32) {
    let x_cum = x_cum;
    if sorted_chroms.is_empty() {
        return (1, 1);
    }
    let first = sorted_chroms[0];
    if x_cum <= *offsets.get(&first).unwrap() {
        return (first, 1);
    }
    let last = *sorted_chroms.last().unwrap();
    let last_lo = *offsets.get(&last).unwrap();
    let last_sz = *chr_sizes.get(&last).unwrap();
    if x_cum >= last_lo + last_sz {
        return (last, (last_sz as i32).max(1));
    }
    for (i, &c) in sorted_chroms.iter().enumerate() {
        let lo = *offsets.get(&c).unwrap();
        let next_lo = if i + 1 < sorted_chroms.len() {
            *offsets.get(&sorted_chroms[i + 1]).unwrap()
        } else {
            lo + chr_sizes.get(&c).unwrap() + 1.0
        };
        if lo <= x_cum && x_cum < next_lo {
            let mut pos = (x_cum - lo).round() as i32;
            let sz = *chr_sizes.get(&c).unwrap() as i32;
            pos = pos.max(1).min(sz);
            return (c, pos);
        }
    }
    let sz = *chr_sizes.get(&last).unwrap() as i32;
    (last, sz.max(1))
}

fn format_pos_compact(pos: i32) -> String {
    if pos >= 1_000_000 {
        format!("{:.1}M", pos as f64 / 1e6)
    } else if pos >= 10_000 {
        format!("{:.0}k", pos as f64 / 1e3)
    } else {
        pos.to_string()
    }
}

pub fn viewport_chr_label(
    x_lo: f64,
    x_hi: f64,
    offsets: &BTreeMap<i32, f64>,
    chr_sizes: &BTreeMap<i32, f64>,
) -> String {
    let sorted_chroms: Vec<i32> = chr_sizes.keys().copied().collect();
    let (ca, pa) = cumulative_to_chr_pos(x_lo, &sorted_chroms, offsets, chr_sizes);
    let (cb, pb) = cumulative_to_chr_pos(x_hi, &sorted_chroms, offsets, chr_sizes);
    format!("{}:{}-{}:{}", chr_token(ca), pa, chr_token(cb), pb)
}

fn format_y_tick(v: f64) -> String {
    format!("{v:.1}")
}

fn nice_step(rough_step: f64) -> f64 {
    if rough_step <= 0.0 || !rough_step.is_finite() {
        return 0.1;
    }
    let exp = rough_step.log10().floor();
    let base = rough_step / 10f64.powf(exp);
    let nice = if base <= 1.0 {
        1.0
    } else if base <= 2.0 {
        2.0
    } else if base <= 5.0 {
        5.0
    } else {
        10.0
    };
    nice * 10f64.powf(exp)
}

fn dynamic_y_ticks(y_min: f64, y_max: f64, plot_height: usize) -> Vec<f64> {
    if !y_min.is_finite() || !y_max.is_finite() {
        return if y_min.is_finite() { vec![y_min] } else { vec![0.0] };
    }
    if y_max <= y_min {
        let step = nice_step(0.1);
        return vec![y_min, y_min + step];
    }
    let span = y_max - y_min;
    let n_target = (plot_height / 4).max(3).min(9).max(3);
    let step = nice_step(span / (n_target.saturating_sub(1).max(2)) as f64);
    let top = y_min + ((span * 1.02) / step).ceil() * step;
    let mut ticks = Vec::new();
    let mut t = y_min;
    let mut guard = 0;
    while t <= top + 1e-9 && guard < 16 {
        ticks.push((t * 1e8).round() / 1e8);
        t += step;
        guard += 1;
    }
    if ticks.is_empty() {
        vec![y_min]
    } else {
        ticks
    }
}

fn chr_start_boundaries_in_view(
    x_min: f64,
    x_max: f64,
    sorted_chroms: &[i32],
    offsets: &BTreeMap<i32, f64>,
) -> Vec<(f64, i32)> {
    let mut out = Vec::new();
    for &c in sorted_chroms.iter().skip(1) {
        let bx = *offsets.get(&c).unwrap();
        if x_min < bx && bx < x_max {
            out.push((bx, c));
        }
    }
    out
}

fn boundary_cx(bx: f64, x_min: f64, x_max: f64, x_plot0: usize, w: usize) -> usize {
    let px = if x_max <= x_min {
        0.0
    } else {
        ((bx - x_min) / (x_max - x_min)).clamp(0.0, 1.0)
    };
    x_plot0 + (px * w as f64) as usize
}

fn draw_chr_transition_vlines(
    canvas: &mut TerminalCanvas,
    boundaries: &[(f64, i32)],
    x_min: f64,
    x_max: f64,
    x_plot0: usize,
    w: usize,
    y_ranges: &[(usize, usize)],
    unicode: bool,
) {
    let v = if unicode { '│' } else { '|' };
    for &(bx, _) in boundaries {
        let cx = boundary_cx(bx, x_min, x_max, x_plot0, w);
        for &(y_lo, y_hi) in y_ranges {
            let y0 = y_lo.max(0);
            let y1 = y_hi.min(canvas.height - 1);
            for y in y0..=y1 {
                canvas.set(cx, y, v);
            }
        }
    }
}

fn mark_chr_boundaries_on_axis(
    canvas: &mut TerminalCanvas,
    boundaries: &[(f64, i32)],
    x_min: f64,
    x_max: f64,
    x_plot0: usize,
    w: usize,
    axis_y: usize,
    unicode: bool,
) {
    let cross = if unicode { '┼' } else { '+' };
    for &(bx, _) in boundaries {
        let cx = boundary_cx(bx, x_min, x_max, x_plot0, w);
        if axis_y < canvas.height {
            canvas.set(cx, axis_y, cross);
        }
    }
}

fn draw_dynamic_x_ticks(
    canvas: &mut TerminalCanvas,
    x_min: f64,
    x_max: f64,
    global_x_min: f64,
    global_x_max: f64,
    sorted_chroms: &[i32],
    offsets: &BTreeMap<i32, f64>,
    chr_sizes: &BTreeMap<i32, f64>,
    x_plot0: usize,
    w: usize,
    axis_y: usize,
    label_y: usize,
) {
    let full_span = global_x_max - global_x_min;
    let vis_span = x_max - x_min;
    let frac = if full_span > 0.0 {
        vis_span / full_span
    } else {
        1.0
    };

    let mut visible_chroms: Vec<i32> = Vec::new();
    for &c in sorted_chroms {
        let lo = *offsets.get(&c).unwrap();
        let hi = lo + chr_sizes.get(&c).unwrap();
        if hi >= x_min && lo <= x_max {
            visible_chroms.push(c);
        }
    }

    let one_chrom = visible_chroms.len() == 1;
    let use_chr_centers = (!one_chrom) && frac >= 0.12;

    let mut ticks: Vec<(f64, String)> = Vec::new();
    if use_chr_centers {
        for &c in &visible_chroms {
            let lo = *offsets.get(&c).unwrap();
            let center = lo + chr_sizes.get(&c).unwrap() / 2.0;
            if x_min <= center && center <= x_max {
                ticks.push((center, chr_token(c)));
            }
        }
    }
    if ticks.is_empty() {
        let n = (w / 14).max(2).min(7).max(1);
        for i in 0..n {
            let t = if n > 1 {
                i as f64 / (n - 1) as f64
            } else {
                0.5
            };
            let x_c = x_min + t * (x_max - x_min);
            let (cc, pp) = cumulative_to_chr_pos(x_c, sorted_chroms, offsets, chr_sizes);
            ticks.push((x_c, format!("{}:{}", chr_token(cc), format_pos_compact(pp))));
        }
    }

    let mut used_cols = std::collections::HashSet::new();
    ticks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    for (x_c, label) in ticks {
        let px = if x_max <= x_min {
            0.0
        } else {
            ((x_c - x_min) / (x_max - x_min)).clamp(0.0, 1.0)
        };
        let tick_x = x_plot0 + (px * w as f64) as usize;
        if used_cols.contains(&tick_x) {
            continue;
        }
        used_cols.insert(tick_x);
        canvas.set(tick_x, axis_y, '+');
        let start_x = (tick_x.saturating_sub(label.chars().count() / 2)).max(x_plot0);
        for (j, ch) in label.chars().enumerate() {
            if start_x + j < canvas.width {
                canvas.set(start_x + j, label_y, ch);
            }
        }
    }
}

fn draw_y_tick_labels(
    canvas: &mut TerminalCanvas,
    ticks: &[f64],
    y_min: f64,
    y_scale: f64,
    y0: usize,
    h: usize,
) {
    // `h` is the vertical scale in *rows* from the top of the data area (row 1) to
    // the bottom (row y0), i.e. y0 - 1 when the title occupies row 0. Do not use
    // y0 as h: that maps py=1 to row 0 and misaligns ticks (clamped to row 1) vs points.
    let denom = y_scale - y_min;
    let mut last_cy: Option<usize> = None;
    for &yt in ticks {
        let py_ratio = if denom <= 0.0 {
            0.0
        } else {
            ((yt - y_min) / denom).clamp(0.0, 1.0)
        };
        let mut cy = if h == 0 {
            y0
        } else {
            y0.saturating_sub((py_ratio * h as f64) as usize)
        };
        cy = cy.max(1).min(y0);
        if let Some(lc) = last_cy {
            if lc == cy {
                continue;
            }
        }
        last_cy = Some(cy);
        let label = format_y_tick(yt);
        let label: String = format!("{:>6}", label.chars().take(6).collect::<String>());
        let label: String = label.chars().take(6).collect();
        for (i, ch) in label.chars().enumerate() {
            canvas.set(i, cy, ch);
        }
    }
}

fn density_glyph(count: usize, unicode: bool) -> char {
    if count <= 1 {
        if unicode {
            '●'
        } else {
            '*'
        }
    } else if count == 2 {
        if unicode {
            '◉'
        } else {
            'O'
        }
    } else if unicode {
        '█'
    } else {
        '#'
    }
}

pub fn density_legend(unicode: bool) -> &'static str {
    if unicode {
        "density 1x ●  2x ◉  3+x █"
    } else {
        "density 1x *  2x O  3+x #"
    }
}

pub fn sig_threshold_legend(unicode: bool, sig_level: f64) -> String {
    let dash = if unicode { '╌' } else { '=' };
    let p = format_sig_g(sig_level);
    format!("sig {} P={p} (genome-wide threshold)", dash.to_string().repeat(4))
}

fn format_sig_g(v: f64) -> String {
    // Match Python `f"{float(v):g}"` for the default GWAS threshold.
    if (v - 5e-8).abs() < 1e-30 {
        return "5e-08".to_string();
    }
    format!("{v:.0e}")
}

fn draw_gene_track(
    canvas: &mut TerminalCanvas,
    genes: &[(f64, f64, String)],
    x_min: f64,
    x_max: f64,
    x_plot0: usize,
    w: usize,
    lane_rows: &[usize],
    unicode: bool,
) {
    if genes.is_empty() || lane_rows.is_empty() {
        return;
    }
    let seg_ch = if unicode { '━' } else { '-' };
    let left_cap = if unicode { '├' } else { '|' };
    let right_cap = if unicode { '┤' } else { '|' };
    let lane_pad = 1usize;
    let short_label_min_inner = 4usize;

    let mut lane_last_end: Vec<i32> = vec![-1_000_000_000; lane_rows.len()];
    let mut label_end_by_row: HashMap<usize, i32> = HashMap::new();
    let mut pixel_genes: Vec<(usize, usize, String)> = Vec::new();

    for (g_start, g_end, g_name) in genes.iter().take(80) {
        if *g_end < x_min || *g_start > x_max {
            continue;
        }
        let px0 = if x_max <= x_min {
            0.0
        } else {
            (g_start.max(x_min) - x_min) / (x_max - x_min)
        };
        let px1 = if x_max <= x_min {
            0.0
        } else {
            (g_end.min(x_max) - x_min) / (x_max - x_min)
        };
        let x0 = x_plot0 + (px0.clamp(0.0, 1.0) * w as f64) as usize;
        let mut x1 = x_plot0 + (px1.clamp(0.0, 1.0) * w as f64) as usize;
        if x1 <= x0 {
            x1 = (x0 + 1).min(canvas.width - 1);
        }
        pixel_genes.push((x0, x1, g_name.clone()));
    }

    pixel_genes.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));

    let lane_min = *lane_rows.iter().min().unwrap_or(&0);
    let lane_max = *lane_rows.iter().max().unwrap_or(&0);

    for (x0, x1, g_name) in pixel_genes {
        let mut lane_idx = 0usize;
        while lane_idx < lane_rows.len() && (x0 as i32) <= lane_last_end[lane_idx] + lane_pad as i32 {
            lane_idx += 1;
        }
        if lane_idx >= lane_rows.len() {
            continue;
        }
        lane_last_end[lane_idx] = x1 as i32;
        let y = lane_rows[lane_idx];
        canvas.set(x0, y, left_cap);
        canvas.set(x1, y, right_cap);
        for x in (x0 + 1)..x1 {
            canvas.set(x, y, seg_ch);
        }
        let inner_start = x0 + 1;
        let inner_end = x1.saturating_sub(1);
        let inner_width = inner_end as isize - inner_start as isize + 1;
        if inner_width <= 0 {
            continue;
        }
        if inner_width >= short_label_min_inner as isize {
            let max_len = 12.min(inner_width as usize);
            let label: String = g_name.chars().take(max_len).collect();
            let label_x = inner_start + (inner_width as usize).saturating_sub(label.chars().count()) / 2;
            for (i, ch) in label.chars().enumerate() {
                let x = label_x + i;
                if x >= inner_start && x <= inner_end {
                    canvas.set(x, y, ch);
                }
            }
            continue;
        }

        let label: String = g_name.chars().take(12).collect();
        let center = (x0 + x1) / 2;
        let mut label_x = (x_plot0 + 1).max(center.saturating_sub(label.chars().count() / 2));
        label_x = label_x.min(canvas.width.saturating_sub(label.chars().count() + 2));
        let mut label_y: Option<usize> = None;
        for &cand in &[y.saturating_sub(1), y + 1] {
            if cand < lane_min || cand > lane_max {
                continue;
            }
            if label_x as i32 > *label_end_by_row.get(&cand).unwrap_or(&-1_000_000_000) {
                label_y = Some(cand);
                break;
            }
        }
        let Some(ly) = label_y else {
            continue;
        };
        for (i, ch) in label.chars().enumerate() {
            let x = label_x + i;
            if x >= x_plot0 + 1 && x <= canvas.width - 2 {
                canvas.set(x, ly, ch);
            }
        }
        label_end_by_row.insert(ly, (label_x + label.chars().count()).saturating_sub(1) as i32);
    }
}

fn draw_gene_panel_frame(
    canvas: &mut TerminalCanvas,
    x_plot0: usize,
    x_right: usize,
    panel_top: usize,
    panel_bottom: usize,
    unicode: bool,
) {
    if panel_top >= panel_bottom {
        return;
    }
    let h = if unicode { '─' } else { '-' };
    let v = if unicode { '│' } else { '|' };
    let tl = if unicode { '┌' } else { '+' };
    let tr = if unicode { '┐' } else { '+' };
    let bl = if unicode { '└' } else { '+' };
    let br = if unicode { '┘' } else { '+' };
    for x in x_plot0..=x_right {
        canvas.set(x, panel_top, h);
        canvas.set(x, panel_bottom, h);
    }
    for y in panel_top..=panel_bottom {
        canvas.set(x_plot0, y, v);
        canvas.set(x_right, y, v);
    }
    canvas.set(x_plot0, panel_top, tl);
    canvas.set(x_right, panel_top, tr);
    canvas.set(x_plot0, panel_bottom, bl);
    canvas.set(x_right, panel_bottom, br);
    for (i, ch) in " Genes ".chars().enumerate() {
        let pos = x_plot0 + 2 + i;
        if pos < x_right {
            canvas.set(pos, panel_top, ch);
        }
    }
}

fn draw_lead_annotation(
    canvas: &mut TerminalCanvas,
    lead_x: f64,
    lead_y_ratio: f64,
    lead_label: &str,
    x_min: f64,
    x_max: f64,
    x_plot0: usize,
    w: usize,
    y0: usize,
    h: usize,
) {
    let px = if x_max <= x_min {
        0.0
    } else {
        ((lead_x - x_min) / (x_max - x_min)).clamp(0.0, 1.0)
    };
    let cx = x_plot0 + (px * w as f64) as usize;
    let cy = if h == 0 {
        y0
    } else {
        y0.saturating_sub((lead_y_ratio.clamp(0.0, 1.0) * h as f64) as usize)
    };
    let label: String = lead_label.chars().take(20).collect();
    let label_y = cy.saturating_sub(1).max(1);
    let half = label.chars().count() / 2;
    let start_x = (cx.saturating_sub(half))
        .max(x_plot0)
        .min(x_plot0.max(canvas.width.saturating_sub(label.chars().count() + 1)));
    for (i, ch) in label.chars().enumerate() {
        if start_x + i < canvas.width {
            canvas.set(start_x + i, label_y, ch);
        }
    }
}

pub fn render_manhattan(
    prepared: &PlotDataset,
    width: usize,
    height: usize,
    sig_level: f64,
    ymax: Option<f64>,
    x_start: Option<f64>,
    x_end: Option<f64>,
    title: Option<&str>,
    unicode: bool,
    y_min: f64,
    lead_variant: Option<(f64, f64, String)>,
    gene_track: Option<&[(f64, f64, String)]>,
    force_gene_panel: bool,
    visible_rows: Option<&[bool]>,
    color: bool,
    light_theme: bool,
) -> Result<String, String> {
    let data = prepared;

    let offsets = &data.layout.offsets;
    let chr_sizes = &data.layout.chr_sizes;
    let global_x_min = data.layout.x_min;
    let global_x_max = data.layout.x_max;
    let global_y_hi = data
        .mlog10p
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .fold(f64::NAN, f64::max);
    let global_y_hi = if global_y_hi.is_nan() {
        y_min
    } else {
        global_y_hi
    };
    let sig_logp = -sig_level.log10();

    let mut x_min = x_start.unwrap_or(global_x_min).max(global_x_min);
    let mut x_max = x_end.unwrap_or(global_x_max).min(global_x_max);
    if x_max <= x_min {
        x_min = global_x_min;
        x_max = global_x_max;
    }

    let title = if let Some(t) = title {
        t.to_string()
    } else {
        format!(
            "Manhattan {} | skip>={y_min:.1}",
            viewport_chr_label(x_min, x_max, offsets, chr_sizes)
        )
    };

    let mask: Vec<bool> = if let Some(vr) = visible_rows {
        if vr.len() != data.x.len() {
            return Err("Visible row mask does not match prepared dataset length.".to_string());
        }
        vr.to_vec()
    } else {
        crate::plot_state::visible_mask(data, x_min, x_max)
    };

    let y_floor = if ymax.is_none() {
        let vis_hi = mask
            .iter()
            .enumerate()
            .filter(|(_, m)| **m)
            .map(|(i, _)| data.mlog10p[i])
            .filter(|v| v.is_finite())
            .fold(f64::NAN, f64::max);
        if vis_hi.is_nan() {
            global_y_hi
        } else {
            vis_hi
        }
    } else {
        ymax.unwrap()
    };
    let y_cap = y_floor.max(sig_logp).max(y_min + 1e-9);

    let style = CanvasStyle { unicode };
    let mut canvas = TerminalCanvas::new(width, height, style)?;
    canvas.draw_axes();
    canvas.label_top_pair(
        &title,
        &format!("gwaspeek v{}", package_version()),
        2,
    );

    let x_plot0 = 7usize;
    let axis_y = canvas.height - 2;
    let label_y = canvas.height - 1;
    let has_gene_track = force_gene_panel || gene_track.map(|g| !g.is_empty()).unwrap_or(false);
    let x_right = canvas.width - 1;
    let (panel_top, panel_bottom, gene_track_eff) = if has_gene_track && canvas.height >= 14 {
        let panel_bottom = axis_y - 1;
        let panel_top = panel_bottom - 3;
        (Some(panel_top), Some(panel_bottom), gene_track)
    } else {
        (None, None, None)
    };
    let _has_gene_panel_layout = panel_top.is_some();
    let y0 = if let Some(pt) = panel_top {
        pt - 2
    } else {
        axis_y - 1
    };
    // Vertical scale: data rows are 1..=y0 (row 0 is the title). Using h = y0 would
    // map the maximum y to row 0, misaligning y-tick labels (min row 1) and the plot.
    let h = y0.saturating_sub(1);
    let w = canvas.width - x_plot0 - 1;
    let y_ticks = dynamic_y_ticks(y_min, y_cap, y0);
    let y_scale = if !y_ticks.is_empty() {
        y_ticks
            .iter()
            .copied()
            .fold(y_cap, f64::max)
            .max(sig_logp)
    } else {
        y_cap.max(sig_logp)
    };
    let denom = y_scale - y_min;

    let sorted_chroms: Vec<i32> = chr_sizes.keys().copied().collect();
    let chr_order: HashMap<i32, usize> = sorted_chroms
        .iter()
        .enumerate()
        .map(|(i, c)| (*c, i))
        .collect();
    let chr_boundaries = chr_start_boundaries_in_view(x_min, x_max, &sorted_chroms, offsets);

    struct CellAgg {
        count: usize,
        chrom: i32,
    }
    let mut cells: HashMap<(usize, usize), CellAgg> = HashMap::new();

    #[derive(Clone)]
    struct CellPt {
        flat: i64,
        y_clip: f64,
        chrom: i32,
    }
    let mut pts: Vec<CellPt> = Vec::new();
    for (i, m) in mask.iter().enumerate() {
        if !*m {
            continue;
        }
        let vx = data.x[i];
        let vy = data.mlog10p[i];
        let vc = data.chrom[i];
        let vy_c = vy.min(y_scale);
        let px = if x_max == x_min {
            0.0
        } else {
            ((vx - x_min) / (x_max - x_min)).clamp(0.0, 1.0)
        };
        let py = if denom <= 0.0 {
            0.0
        } else {
            ((vy_c - y_min) / denom).clamp(0.0, 1.0)
        };
        let cx = x_plot0 + ((px * w as f64).floor() as usize);
        let cy = if h == 0 {
            y0
        } else {
            y0.saturating_sub((py * h as f64).floor() as usize)
        };
        let f = ((cy * canvas.width) + cx) as i64;
        pts.push(CellPt {
            flat: f,
            y_clip: vy_c,
            chrom: vc,
        });
    }
    if !pts.is_empty() {
        pts.sort_by(|a, b| {
            a.flat
                .cmp(&b.flat)
                .then_with(|| a.y_clip.partial_cmp(&b.y_clip).unwrap())
        });
        let mut i = 0usize;
        while i < pts.len() {
            let j0 = i;
            let f0 = pts[i].flat;
            while i < pts.len() && pts[i].flat == f0 {
                i += 1;
            }
            let last = &pts[i - 1];
            let count = i - j0;
            let cell_flat = f0 as usize;
            let cell_y = cell_flat / canvas.width;
            let cell_x = cell_flat % canvas.width;
            cells.insert(
                (cell_x, cell_y),
                CellAgg {
                    count,
                    chrom: last.chrom,
                },
            );
        }
    }

    for ((cx, cy), agg) in cells.iter() {
        let g = density_glyph(agg.count, unicode);
        let painted = ansi_paint_glyph(g, agg.chrom, &chr_order, color, light_theme);
        if color {
            canvas.set_cell(*cx, *cy, &painted);
        } else {
            canvas.set(*cx, *cy, g);
        }
    }

    if let Some((lead_x, lead_y, lead_label)) = lead_variant {
        let lead_py = if denom <= 0.0 {
            0.0
        } else {
            ((lead_y.min(y_scale) - y_min) / denom).clamp(0.0, 1.0)
        };
        draw_lead_annotation(
            &mut canvas,
            lead_x,
            lead_py,
            lead_label.as_str(),
            x_min,
            x_max,
            x_plot0,
            w,
            y0,
            h,
        );
    }

    let sig_ratio = if denom <= 0.0 {
        0.0
    } else {
        ((sig_logp - y_min) / denom).clamp(0.0, 1.0)
    };
    let line_y = if h == 0 {
        y0
    } else {
        y0.saturating_sub((sig_ratio * h as f64) as usize)
    };
    let line_ch = if unicode { '╌' } else { '=' };
    for x in x_plot0..canvas.width - 1 {
        canvas.set(x, line_y, line_ch);
    }

    if let (Some(pt), Some(pb)) = (panel_top, panel_bottom) {
        draw_gene_panel_frame(&mut canvas, x_plot0, x_right, pt, pb, unicode);
        let lane_rows: Vec<usize> = (pt + 1..pb).collect();
        if let Some(gt) = gene_track_eff {
            draw_gene_track(
                &mut canvas,
                gt,
                x_min,
                x_max,
                x_plot0,
                w,
                &lane_rows,
                unicode,
            );
        }
        draw_chr_transition_vlines(
            &mut canvas,
            &chr_boundaries,
            x_min,
            x_max,
            x_plot0,
            w,
            &[(pt + 1, pb - 1)],
            unicode,
        );
    }

    draw_dynamic_x_ticks(
        &mut canvas,
        x_min,
        x_max,
        global_x_min,
        global_x_max,
        &sorted_chroms,
        offsets,
        chr_sizes,
        x_plot0,
        w,
        axis_y,
        label_y,
    );
    mark_chr_boundaries_on_axis(
        &mut canvas,
        &chr_boundaries,
        x_min,
        x_max,
        x_plot0,
        w,
        axis_y,
        unicode,
    );

    draw_y_tick_labels(&mut canvas, &y_ticks, y_min, y_scale, y0, h);
    Ok(canvas.render())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::load_sumstats;
    use crate::plot_state::prepare_plot_dataset;
    use crate::preprocess::{preprocess_sumstats, RawSumstats};
    use std::path::PathBuf;

    fn fixture(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/fixtures")
            .join(name)
    }

    fn prepared_small() -> crate::plot_state::PlotDataset {
        let raw = load_sumstats(&fixture("sumstats_small.tsv"), "\t", None, None, None, None).unwrap();
        let clean = preprocess_sumstats(raw, 3.0).unwrap();
        prepare_plot_dataset(&clean, Some("37"), false).unwrap()
    }

    #[test]
    fn manhattan_render_smoke() {
        let p = prepared_small();
        let txt = render_manhattan(&p, 60, 18, 5e-8, None, None, None, None, false, 3.0, None, None, false, None, false, false).unwrap();
        assert!(txt.contains("Manhattan"));
        assert!(txt.contains("gwaspeek"));
        assert!(txt.contains("skip>="));
        assert!(txt.contains('*'));
    }

    #[test]
    fn two_hits_same_cell_ascii_o_glyph() {
        let raw = RawSumstats {
            chrom_tokens: vec!["1".into(), "1".into()],
            pos_tokens: vec!["100".into(), "100".into()],
            p_tokens: vec!["1e-8".into(), "1e-8".into()],
            mlog10p_tokens: vec![],
            has_p_column: true,
        };
        let clean = preprocess_sumstats(raw, 0.0).unwrap();
        let prep = prepare_plot_dataset(&clean, Some("37"), false).unwrap();
        let txt = render_manhattan(&prep, 40, 12, 5e-8, None, None, None, None, false, 0.0, None, None, false, None, false, false).unwrap();
        assert!(txt.contains('O'));
    }

    #[test]
    fn sig_legend_matches_plot_glyph_ascii() {
        let s = sig_threshold_legend(false, 5e-8);
        assert!(s.contains('='));
        assert!(s.contains("5e-08"));
    }
}
