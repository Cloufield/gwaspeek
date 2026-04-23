use unicode_width::UnicodeWidthChar;

/// Display width in terminal columns, treating CSI sequences (e.g. `\x1b[36m`) as width 0.
fn display_width_ansi(s: &str) -> usize {
    let mut w = 0usize;
    let mut it = s.chars().peekable();
    while let Some(c) = it.next() {
        if c == '\x1b' && it.peek() == Some(&'[') {
            it.next(); // '['
            while let Some(c2) = it.next() {
                if ('\x40'..='\x7e').contains(&c2) {
                    break;
                }
            }
            continue;
        }
        w += UnicodeWidthChar::width(c).unwrap_or(0);
    }
    w
}

/// Trim trailing whitespace, then fit to `target` display columns (truncate if needed, pad with spaces).
fn fit_line_display_width(s: &str, target: usize) -> String {
    let trimmed = s.trim_end();
    let mut out = String::new();
    let mut w = 0usize;
    let mut it = trimmed.chars().peekable();
    while let Some(c) = it.next() {
        if c == '\x1b' && it.peek() == Some(&'[') {
            out.push(c);
            out.push(it.next().unwrap());
            while let Some(c2) = it.next() {
                out.push(c2);
                if ('\x40'..='\x7e').contains(&c2) {
                    break;
                }
            }
            continue;
        }
        let cw = UnicodeWidthChar::width(c).unwrap_or(0);
        if w + cw > target {
            break;
        }
        w += cw;
        out.push(c);
    }
    let dw = display_width_ansi(&out);
    for _ in 0..target.saturating_sub(dw) {
        out.push(' ');
    }
    out
}

#[derive(Clone, Debug)]
pub struct CanvasStyle {
    pub unicode: bool,
}

impl CanvasStyle {
    pub fn point(&self) -> char {
        if self.unicode {
            '●'
        } else {
            'o'
        }
    }

    pub fn axis_h(&self) -> char {
        if self.unicode {
            '─'
        } else {
            '-'
        }
    }

    pub fn axis_v(&self) -> char {
        if self.unicode {
            '│'
        } else {
            '|'
        }
    }
}

pub struct TerminalCanvas {
    pub width: usize,
    pub height: usize,
    pub style: CanvasStyle,
    grid: Vec<Vec<String>>,
}

impl TerminalCanvas {
    pub fn new(width: usize, height: usize, style: CanvasStyle) -> Result<Self, String> {
        if width < 20 || height < 8 {
            return Err("Canvas too small, use width>=20 and height>=8".to_string());
        }
        let grid = vec![vec![" ".to_string(); width]; height];
        Ok(Self {
            width,
            height,
            style,
            grid,
        })
    }

    pub fn set(&mut self, x: usize, y: usize, ch: char) {
        if x < self.width && y < self.height {
            self.grid[y][x] = ch.to_string();
        }
    }

    pub fn set_cell(&mut self, x: usize, y: usize, s: &str) {
        if x < self.width && y < self.height {
            self.grid[y][x] = s.to_string();
        }
    }

    pub fn draw_axes(&mut self) {
        let x0 = 6;
        let y0 = self.height - 2;
        for x in x0..self.width {
            self.set(x, y0, self.style.axis_h());
        }
        for y in 0..=y0 {
            self.set(x0, y, self.style.axis_v());
        }
        self.set(x0, y0, '+');
    }

    pub fn label_top_pair(&mut self, left: &str, right: &str, gap: usize) {
        let mut line: Vec<String> = vec![" ".to_string(); self.width];
        let right_label: String = right.chars().take(self.width).collect();
        let right_start = self.width.saturating_sub(right_label.chars().count());
        for (i, ch) in right_label.chars().enumerate() {
            if right_start + i < self.width {
                line[right_start + i] = ch.to_string();
            }
        }
        let max_left_width = right_start.saturating_sub(gap.max(0));
        let left_label: String = if max_left_width > 0 {
            left.chars().take(max_left_width).collect::<String>().trim_end().to_string()
        } else {
            String::new()
        };
        for (i, ch) in left_label.chars().enumerate() {
            line[i] = ch.to_string();
        }
        for i in 0..self.width {
            self.grid[0][i] = line[i].clone();
        }
    }

    pub fn render(&self) -> String {
        self.grid
            .iter()
            .map(|row| {
                let raw: String = row.iter().cloned().collect();
                // Per-cell layout assumes one column per cell, but ANSI is zero-width and some
                // glyphs are wide; pad/truncate to a fixed display width so rows don't wrap.
                fit_line_display_width(&raw, self.width)
            })
            .collect::<Vec<_>>()
            .join("\n")
    }
}

#[cfg(test)]
mod width_tests {
    use super::*;

    #[test]
    fn ansi_glyph_counts_as_one_column_for_fit() {
        let raw = format!("{}a\x1b[36m*\x1b[0mb", " ".repeat(7));
        let line = fit_line_display_width(&raw, 10);
        assert_eq!(display_width_ansi(&line), 10);
        assert!(line.contains('\x1b'));
    }

    #[test]
    fn render_row_pads_to_canvas_width() {
        let style = CanvasStyle { unicode: false };
        let mut c = TerminalCanvas::new(40, 10, style).unwrap();
        c.set(0, 5, 'x');
        let out = c.render();
        let row5 = out.lines().nth(5).unwrap();
        assert_eq!(display_width_ansi(row5), 40);
    }
}
