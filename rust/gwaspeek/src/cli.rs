use std::io::Write;
use std::path::PathBuf;

use clap::Parser;

use crate::interactive::run_interactive_manhattan;
use crate::io::load_sumstats;
use crate::manhattan::{render_manhattan, stdout_color_supported};
use crate::plot_state::prepare_plot_dataset;
use crate::preprocess::preprocess_sumstats;
#[derive(Parser, Debug)]
#[command(
    name = "gwaspeek",
    version = env!("CARGO_PKG_VERSION"),
    about = "Draw GWAS Manhattan plots in the terminal: static snapshot (-s) or full-screen interactive pan/zoom.",
    after_help = "Examples:\n  gwaspeek tests/fixtures/sumstats_small.tsv\n  gwaspeek -i tests/fixtures/sumstats_small.tsv\n  gwaspeek -s tests/fixtures/sumstats_small.tsv\n"
)]
pub struct Args {
    /// GWAS summary statistics file (TSV/CSV). Same as -i FILE when this is the only input.
    #[arg(value_name = "FILE")]
    pub sumstats: Option<PathBuf>,

    /// Print one Manhattan frame to stdout and exit (non-interactive)
    #[arg(short = 's', long = "static", value_name = "FILE")]
    pub static_file: Option<PathBuf>,

    /// Open the interactive Manhattan viewer (explicit form of positional FILE)
    #[arg(short = 'i', long = "interactive", value_name = "FILE")]
    pub interactive_file: Option<PathBuf>,

    /// Field delimiter in the input file (default: tab)
    #[arg(long, default_value = "\t")]
    pub sep: String,

    /// Chromosome column name (auto-detect from bundled formatbook aliases if omitted)
    #[arg(long = "chr", value_name = "NAME")]
    pub chr_col: Option<String>,

    /// Genomic position column name (auto-detect if omitted)
    #[arg(long = "pos", value_name = "NAME")]
    pub pos_col: Option<String>,

    /// P-value column name (auto-detect; do not combine with --mlog10p)
    #[arg(long = "p", value_name = "NAME")]
    pub p_col: Option<String>,

    /// Per-variant -log10(P) column name (auto-detect if P is absent)
    #[arg(long = "mlog10p", value_name = "NAME")]
    pub mlog10p_col: Option<String>,

    /// Hide variants with -log10(P) below this threshold (also used as the y-axis floor)
    #[arg(long, default_value_t = 3.0)]
    pub skip: f64,

    /// Terminal width in characters for the plot frame
    #[arg(long, default_value_t = 100)]
    pub width: usize,

    /// Terminal height in lines for the plot frame
    #[arg(long, default_value_t = 28)]
    pub height: usize,

    /// Use ASCII drawing characters instead of Unicode box drawing and glyphs
    #[arg(long, action = clap::ArgAction::SetTrue)]
    pub ascii: bool,

    /// Reference assembly for cytoband lengths and cumulative genome layout
    #[arg(long, value_parser = ["37", "38"], default_value = "37")]
    pub build: String,

    /// Disable ANSI colors (interactive status/footer and static plot chromosome colors)
    #[arg(long = "no-color", action = clap::ArgAction::SetTrue)]
    pub no_color: bool,

    /// Non-human / non-GRCh layout
    #[arg(long = "nh", alias = "not-human", action = clap::ArgAction::SetTrue)]
    pub non_human: bool,

    /// Genome-wide significance P-value; drawn as a horizontal threshold in -log10(P) space
    #[arg(long = "sig-level", default_value_t = 5e-8)]
    pub sig_level: f64,

    /// Clamp the plot's maximum -log10(P) to this value (omit for data-driven ceiling)
    #[arg(long = "ymax", value_name = "L")]
    pub ymax: Option<f64>,

    /// GRCh37 protein-coding gene GTF(.gz) for the interactive gene track
    #[arg(long = "gtf")]
    pub gtf: Option<PathBuf>,

    /// GRCh38 protein-coding gene GTF(.gz) for the interactive gene track when build is 38
    #[arg(long = "gtf38")]
    pub gtf38: Option<PathBuf>,
}

fn resolve_input_and_mode(args: &Args) -> Result<(PathBuf, &'static str), String> {
    let pos = args.sumstats.as_ref();
    let static_f = args.static_file.as_ref();
    let inter = args.interactive_file.as_ref();
    let n = [pos, static_f, inter]
        .iter()
        .filter(|x| x.map(|p| !p.as_os_str().is_empty()).unwrap_or(false))
        .count();
    if n == 0 {
        return Err("Provide summary statistics: FILE, -s FILE, or -i FILE".to_string());
    }
    if n > 1 {
        return Err("Use only one of: positional FILE, -s FILE, or -i FILE".to_string());
    }
    if let Some(p) = inter {
        return Ok((p.clone(), "interactive"));
    }
    if let Some(p) = static_f {
        return Ok((p.clone(), "static"));
    }
    Ok((pos.unwrap().clone(), "interactive"))
}

pub fn run() -> Result<(), String> {
    let args = Args::parse();
    let (path, mode) = resolve_input_and_mode(&args)?;

    if mode == "interactive" {
        eprintln!(
            "[gwaspeek v{}] Loading {} ...",
            crate::versioning::package_version(),
            path.display()
        );
        let _ = std::io::stderr().flush();
    }

    let raw = load_sumstats(
        &path,
        &args.sep,
        args.chr_col.as_deref(),
        args.pos_col.as_deref(),
        args.p_col.as_deref(),
        args.mlog10p_col.as_deref(),
    )?;
    let clean = preprocess_sumstats(raw, args.skip)?;
    let use_unicode = !args.ascii;

    if mode == "static" {
        let prepared = prepare_plot_dataset(
            &clean,
            Some(args.build.as_str()),
            args.non_human,
        )?;
        let chr_color = stdout_color_supported() && !args.no_color;
        let out = render_manhattan(
            &prepared,
            args.width,
            args.height,
            args.sig_level,
            args.ymax,
            None,
            None,
            None,
            use_unicode,
            args.skip,
            None,
            None,
            false,
            None,
            chr_color,
            false,
        )?;
        println!("{out}");
        return Ok(());
    }

    run_interactive_manhattan(
        clean,
        args.width,
        args.height,
        args.sig_level,
        args.ymax,
        use_unicode,
        args.skip,
        args.gtf.as_deref(),
        args.gtf38.as_deref(),
        &args.build,
        !args.no_color,
        args.non_human,
    )?;

    Ok(())
}
