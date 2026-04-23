/// Match Python `gwaspeek.versioning.package_version`: installed package version, else `"dev"`.
pub fn package_version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}
