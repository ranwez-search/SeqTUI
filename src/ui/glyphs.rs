#[derive(Clone, Copy, Debug)]
pub struct Glyphs {
    pub dir_prefix: &'static str,
    pub file_prefix: &'static str,
    pub cursor: &'static str,
    pub h_separator: &'static str,
    pub arrow_left: &'static str,
    pub arrow_right: &'static str,
    pub arrow_up: &'static str,
    pub arrow_down: &'static str,
}

pub fn select(fancy_requested: bool) -> Glyphs {
    if fancy_requested {
        fancy()
    } else {
        ascii()
    }
}

fn ascii() -> Glyphs {
    Glyphs {
        dir_prefix: ">",
        file_prefix: " ",
        cursor: "|",
        h_separator: "-",
        arrow_left: "<",
        arrow_right: ">",
        arrow_up: "^",
        arrow_down: "v",
    }
}

fn fancy() -> Glyphs {
    Glyphs {
        dir_prefix: "📁 ",
        file_prefix: "📄 ",
        cursor: "█",
        h_separator: "─",
        arrow_left: "←",
        arrow_right: "→",
        arrow_up: "↑",
        arrow_down: "↓",
    }
}
