use std::{
    collections::HashMap,
    ffi::{OsStr, OsString},
    io,
    path::Path,
    process::Command,
};

#[cfg(target_family = "unix")]
use std::os::unix::ffi::OsStrExt;

#[cfg(target_family = "windows")]
use std::os::windows::ffi::OsStringExt;

#[derive(Debug)]
struct ConfigVariables {
    map: HashMap<String, String>,
}

impl ConfigVariables {
    fn get_r_cmd_config(&self, key: &str) -> String {
        match self.map.get(key) {
            Some(value) => value.to_string(),
            None => String::from(""),
        }
    }
}

// frustratingly, something like the following does not exist in an
// OS-independent way in Rust
#[cfg(target_family = "unix")]
fn byte_array_to_os_string(bytes: &[u8]) -> OsString {
    let os_str = OsStr::from_bytes(bytes);
    os_str.to_os_string()
}

#[link(name = "kernel32")]
#[cfg(target_family = "windows")]
extern "system" {
    #[link_name = "GetConsoleCP"]
    fn get_console_code_page() -> u32;
    #[link_name = "MultiByteToWideChar"]
    fn multi_byte_to_wide_char(
        CodePage: u32,
        dwFlags: u32,
        lpMultiByteStr: *const u8,
        cbMultiByte: i32,
        lpWideCharStr: *mut u16,
        cchWideChar: i32,
    ) -> i32;
}

// convert bytes to wide-encoded characters on Windows
// from: https://stackoverflow.com/a/40456495/4975218
#[cfg(target_family = "windows")]
fn wide_from_console_string(bytes: &[u8]) -> Vec<u16> {
    assert!(bytes.len() < std::i32::MAX as usize);
    let mut wide;
    let mut len;
    unsafe {
        let cp = get_console_code_page();
        len = multi_byte_to_wide_char(
            cp,
            0,
            bytes.as_ptr() as *const u8,
            bytes.len() as i32,
            std::ptr::null_mut(),
            0,
        );
        wide = Vec::with_capacity(len as usize);
        len = multi_byte_to_wide_char(
            cp,
            0,
            bytes.as_ptr() as *const u8,
            bytes.len() as i32,
            wide.as_mut_ptr(),
            len,
        );
        wide.set_len(len as usize);
    }
    wide
}

#[cfg(target_family = "windows")]
fn byte_array_to_os_string(bytes: &[u8]) -> OsString {
    // first, use Windows API to convert to wide encoded
    let wide = wide_from_console_string(bytes);
    // then, use `std::os::windows::ffi::OsStringExt::from_wide()`
    OsString::from_wide(&wide)
}

// Execute an R CMD config and return the captured output
fn r_cmd_config<S: AsRef<OsStr>>(r_binary: S) -> io::Result<OsString> {
    let out = Command::new(r_binary)
        .args(&["CMD", "config", "--all"])
        .output()?;

    // if there are any errors we print them out, helps with debugging
    if !out.stderr.is_empty() {
        println!(
            "> {}",
            byte_array_to_os_string(&out.stderr)
                .as_os_str()
                .to_string_lossy()
        );
    }

    Ok(byte_array_to_os_string(&out.stdout))
}

fn build_r_cmd_configs() -> ConfigVariables {
    let r_configs = r_cmd_config("R");

    let mut rcmd_config_map = HashMap::new();
    match r_configs {
        Ok(configs) => {
            let input = configs.as_os_str().to_string_lossy();
            for line in input.lines() {
                // Ignore lines beyond comment marker
                if line.starts_with("##") {
                    break;
                }
                let parts: Vec<_> = line.split('=').map(str::trim).collect();
                if let [name, value] = parts.as_slice() {
                    rcmd_config_map.insert(name.to_string(), value.to_string());
                }
            }
        }
        _ => (),
    }
    // Return the struct
    ConfigVariables {
        map: rcmd_config_map,
    }
}

fn get_libs_and_paths(strings: Vec<String>) -> (Vec<String>, Vec<String>) {
    let mut paths: Vec<String> = Vec::new();
    let mut libs: Vec<String> = Vec::new();

    for s in &strings {
        let parts: Vec<&str> = s.split_whitespace().collect();
        for part in parts {
            if part.starts_with("-L") {
                paths.push(part[2..].to_string());
            } else if part.starts_with("-l") {
                libs.push(part[2..].to_string());
            }
        }
    }
    (paths, libs)
}

fn main() {
    let r_configs = build_r_cmd_configs();
    let (lib_paths, libs) = get_libs_and_paths(
        [
            r_configs.get_r_cmd_config("BLAS_LIBS"),
            r_configs.get_r_cmd_config("LAPACK_LIBS"),
            r_configs.get_r_cmd_config("FLIBS"),
        ]
        .to_vec(),
    );

    for path in lib_paths.iter() {
        // Some R builds (e.g. homebrew) contain hardwired gfortran12
        // paths, which may or may not exist if one has upgraded
        // gfortran. So filter out non-existent ones, so that cargo
        // doesn't complain.
        if Path::new(path).exists() {
            println!("cargo:rustc-link-search={}", path);
        }
    }

    for lib in libs.iter() {
        println!("cargo:rustc-link-lib=dylib={}", lib);
    }
    println!("cargo:rerun-if-changed=build.rs");
}
