//! Smoke tests: run the real binary with verified argument patterns.

use std::process::Command;

fn bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_geoforge_v2"))
}

#[test]
fn format_text_writes_file() {
    let out = std::env::temp_dir().join("geoforge_v2_smoke_text");
    let _ = std::fs::remove_dir_all(&out);
    let status = bin()
        .args([
            "generate",
            "--seed",
            "42",
            "--format",
            "text",
            "-d",
            out.to_str().unwrap(),
        ])
        .status()
        .expect("spawn");
    assert!(status.success(), "generate failed: {status:?}");
    let txt = out.join("cosmology_seed42.txt");
    assert!(txt.is_file(), "missing {}", txt.display());
    let content = std::fs::read_to_string(&txt).unwrap();
    assert!(content.contains("Geoforge V2 cosmology report"));
}

#[test]
fn export_format_alias_writes_json() {
    let out = std::env::temp_dir().join("geoforge_v2_smoke_json");
    let _ = std::fs::remove_dir_all(&out);
    let status = bin()
        .args([
            "generate",
            "--seed",
            "99",
            "--export-format",
            "json",
            "-d",
            out.to_str().unwrap(),
        ])
        .status()
        .expect("spawn");
    assert!(status.success());
    assert!(out.join("cosmology_seed99.json").is_file());
}

#[test]
fn export_alias_works() {
    let out = std::env::temp_dir().join("geoforge_v2_smoke_export_alias");
    let _ = std::fs::remove_dir_all(&out);
    let status = bin()
        .args([
            "generate",
            "--seed",
            "7",
            "--export",
            "text",
            "-d",
            out.to_str().unwrap(),
        ])
        .status()
        .expect("spawn");
    assert!(status.success(), "export alias failed: {status:?}");
    assert!(out.join("cosmology_seed7.txt").is_file());
}
