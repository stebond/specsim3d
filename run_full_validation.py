#!/usr/bin/env python3
"""SPECSIM3D: complete validation + benchmark + figures in one go.

What it does (everything automatic, results in  <package>/validation/<timestamp>/ ):
  0. records the machine: CPU, RAM, cores, OS, Fortran compiler (ifx / ifort / gfortran) version,
     Python / numpy / scipy / matplotlib / GSTools versions, build date of specsim3d(.exe)  -> system_info.json
  1. tests that do not need GSTools: test_conditional, test_rotation, test_bias, test_vario_fortran
  2. tests that compare with GSTools (and gstlearn simuSpectral, exponential + Gaussian only), run TWICE in separate processes:
        a) "released"  : GSTools exactly as installed (e.g. 1.7.0)
        b) "modified"  : same GSTools + the author's spherical-3D patch (docs/gstools_patch)
     test_gauss, test_cases, compare_gstools  (their csv/png/json outputs are kept for each variant)
  3. make_figures.py (all manual/paper figures at 300 dpi + tables; --full benchmark unless --quick)
  4. summary.md / summary.json with PASS/FAIL, timings and the machine description (copy it in the paper)

Usage (from the package folder):
    py run_full_validation.py             # everything, full benchmark (several minutes, ~3 GB RAM)
    py run_full_validation.py --quick     # reduced benchmark (no 256^3 runs)
    py run_full_validation.py --skip-figures | --skip-gstools | --skip-tests
Close other heavy programs while it runs: the timings are part of the results.
Author: Stefano Bondua
"""
import argparse, datetime, json, os, platform, re, shutil, subprocess, sys, time

ROOT = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.join(ROOT, "tests")
DOCS = os.path.join(ROOT, "docs")
PATCH_DIR = os.path.join(DOCS, "gstools_patch")
EXE = os.path.join(ROOT, "specsim3d.exe" if os.name == "nt" else "specsim3d")

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--quick", action="store_true", help="do not run the full benchmark (make_figures without --full)")
ap.add_argument("--skip-tests", action="store_true", help="skip the tests that do not need GSTools")
ap.add_argument("--skip-gstools", action="store_true", help="skip the GSTools tests (released and modified)")
ap.add_argument("--skip-gstlearn", action="store_true", help="skip the gstlearn comparison (simuSpectral, exponential and Gaussian only)")
ap.add_argument("--skip-figures", action="store_true", help="skip make_figures.py")
ap.add_argument("--out", default=None, help="output folder (default: validation/<timestamp>)")
args = ap.parse_args()

STAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
OUT = args.out or os.path.join(ROOT, "validation", STAMP)
os.makedirs(OUT, exist_ok=True)
ENV = dict(os.environ, PYTHONIOENCODING="utf-8", MPLBACKEND="Agg")
RESULTS = []          # (stage, name, status, seconds)


def log(msg=""):
    print(msg, flush=True)


# ---------------------------------------------------------------------------- machine description
def run_quiet(cmd, timeout=30):
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, errors="replace")
        return (p.stdout + p.stderr).strip(), p.returncode
    except Exception:
        return "", -1


def cpu_name():
    if os.name == "nt":
        o, rc = run_quiet(["powershell", "-NoProfile", "-Command",
                           "(Get-CimInstance Win32_Processor | Select-Object -First 1).Name"])
        if o and rc == 0:
            return o.splitlines()[0].strip()
        o, rc = run_quiet(["wmic", "cpu", "get", "name"])
        lines = [l.strip() for l in o.splitlines() if l.strip() and l.strip().lower() != "name"]
        if lines:
            return lines[0]
    elif os.path.exists("/proc/cpuinfo"):
        m = re.search(r"model name\s*:\s*(.+)", open("/proc/cpuinfo").read())
        if m:
            return m.group(1).strip()
    elif sys.platform == "darwin":
        o, _ = run_quiet(["sysctl", "-n", "machdep.cpu.brand_string"])
        if o:
            return o
    return platform.processor() or "unknown"


def ram_gb():
    try:
        import psutil
        return round(psutil.virtual_memory().total / 2**30, 1)
    except Exception:
        pass
    if os.name == "nt":
        o, _ = run_quiet(["powershell", "-NoProfile", "-Command",
                          "(Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory"])
        if o.strip().isdigit():
            return round(int(o.strip()) / 2**30, 1)
    elif os.path.exists("/proc/meminfo"):
        return round(int(re.search(r"MemTotal:\s+(\d+)", open("/proc/meminfo").read()).group(1)) / 2**20, 1)
    return None


def physical_cores():
    try:
        import psutil
        return psutil.cpu_count(logical=False)
    except Exception:
        pass
    if os.name == "nt":
        o, rc = run_quiet(["powershell", "-NoProfile", "-Command",
                           "(Get-CimInstance Win32_Processor | Measure-Object -Property NumberOfCores -Sum).Sum"])
        if o.strip().isdigit():
            return int(o.strip())
    return None


def compiler_version():
    if os.environ.get("SPECSIM3D_COMPILER"):
        return os.environ["SPECSIM3D_COMPILER"]
    bi = os.path.join(ROOT, "build_info.txt")          # written by build_specsim3d.bat: compiler and flags used for THIS exe
    if os.path.exists(bi):
        txt = " ".join(l.strip() for l in open(bi, encoding="utf-8", errors="replace") if l.strip())
        if txt:
            return txt
    if os.name == "nt" and shutil.which("ifx") is None:        # not started from the oneAPI prompt: look in the install folders
        import glob
        pats = [os.path.join(os.environ.get(v, ""), "Intel", "oneAPI", "compiler", "*", "bin", "ifx.exe")
                for v in ("ProgramFiles(x86)", "ProgramFiles")]
        cands = sorted(sum((glob.glob(p) for p in pats), []))
        if cands:
            o, rc = run_quiet([cands[-1], "/QV"], 30)
            lines = [l.strip() for l in o.splitlines() if l.strip()]
            if lines:
                return lines[0] + "  [" + cands[-1] + "]"
    for cmd in (["ifx", "/QV"], ["ifx", "--version"], ["ifort", "/QV"], ["ifort", "--version"], ["gfortran", "--version"]):
        if shutil.which(cmd[0]) is None:
            continue
        o, rc = run_quiet(cmd, 20)
        if o:
            first = [l for l in o.splitlines() if l.strip()][0].strip()
            return first
    return "not found in PATH (run from the Intel oneAPI command prompt, or set SPECSIM3D_COMPILER)"


def pkg_version(mod):
    code = f"import {mod};print({mod}.__version__)"
    o, rc = run_quiet([sys.executable, "-c", code])
    return o.splitlines()[-1] if rc == 0 and o else None


def machine_info():
    exe_time = datetime.datetime.fromtimestamp(os.path.getmtime(EXE)).strftime("%Y-%m-%d %H:%M") if os.path.exists(EXE) else None
    return {
        "date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "cpu": cpu_name(), "cores_logical": os.cpu_count(), "cores_physical": physical_cores(),
        "ram_gb": ram_gb(), "os": platform.platform(),
        "omp_num_threads_env": os.environ.get("OMP_NUM_THREADS", "(not set: all logical cores)"),
        "fortran_compiler": compiler_version(), "executable": os.path.basename(EXE), "executable_built": exe_time,
        "python": sys.version.split()[0], "numpy": pkg_version("numpy"), "scipy": pkg_version("scipy"),
        "matplotlib": pkg_version("matplotlib"), "gstools": pkg_version("gstools"), "gstlearn": pkg_version("gstlearn"),
    }


# ---------------------------------------------------------------------------- stage runner
def run_stage(stage, name, cmd, logfile, cwd=ROOT):
    """Run cmd, stream its output to the console and to logfile; record PASS/FAIL."""
    log(f"\n=== [{stage}] {name} ===")
    t0 = time.time()
    with open(logfile, "w", encoding="utf-8", errors="replace") as fh:
        p = subprocess.Popen(cmd, cwd=cwd, env=ENV, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             text=True, encoding="utf-8", errors="replace")
        for line in p.stdout:
            fh.write(line)
            sys.stdout.write("    " + line)
        p.wait()
    dt = time.time() - t0
    status = "PASS" if p.returncode == 0 else "FAIL"
    log(f"--- {name}: {status} ({dt:.0f} s)")
    RESULTS.append((stage, name, status, round(dt, 1)))
    return p.returncode == 0


def clean_tests_work():
    """The GSTools tests cache their GSTools fields in tests/work: clear it between variants."""
    shutil.rmtree(os.path.join(TESTS, "work"), ignore_errors=True)


def collect_results(dest):
    src = os.path.join(TESTS, "results")
    if os.path.isdir(src):
        shutil.copytree(src, dest, dirs_exist_ok=True)


def modified_cmd(script):
    """Command that runs a test script with the author's GSTools patch applied (own process)."""
    code = ("import sys, runpy; sys.path.insert(0, r'%s'); import gstools_spherical3d_patch; "
            "s = r'%s'; sys.argv = [s]; runpy.run_path(s, run_name='__main__')" % (PATCH_DIR, script))
    return [sys.executable, "-c", code]


# ---------------------------------------------------------------------------- main
def main():
    log(f"SPECSIM3D full validation - output: {OUT}")
    if not os.path.exists(EXE):
        sys.exit(f"Executable not found: {EXE}\nBuild it first (build_specsim3d.bat).")
    info = machine_info()
    ENV["SPECSIM3D_COMPILER"] = info["fortran_compiler"]        # the figure scripts print it in the footers
    open(os.path.join(ROOT, "compiler.txt"), "w", encoding="utf-8").write(info["fortran_compiler"] + "\n")
    json.dump(info, open(os.path.join(OUT, "system_info.json"), "w"), indent=1)
    log(json.dumps(info, indent=1))
    for need in ("numpy", "scipy", "matplotlib"):
        if not info[need]:
            sys.exit(f"Python package '{need}' missing: py -m pip install numpy scipy matplotlib")
    have_gs = bool(info["gstools"])
    t00 = time.time()

    if not args.skip_tests:
        for t in ("test_conditional", "test_rotation", "test_bias", "test_vario_fortran"):
            run_stage("tests", t, [sys.executable, os.path.join(TESTS, t + ".py")], os.path.join(OUT, t + ".log"))
        collect_results(os.path.join(OUT, "results_nogstools"))

    if not args.skip_gstools:
        if not have_gs:
            log("\nGSTools not installed: skipping the GSTools tests (py -m pip install gstools)")
        else:
            for variant in ("released", "modified"):
                clean_tests_work()
                shutil.rmtree(os.path.join(TESTS, "results"), ignore_errors=True)
                for t in ("test_gauss", "test_cases", "compare_gstools"):
                    script = os.path.join(TESTS, t + ".py")
                    cmd = [sys.executable, script] if variant == "released" else modified_cmd(script)
                    run_stage(f"gstools-{variant}", t, cmd, os.path.join(OUT, f"{t}_{variant}.log"))
                collect_results(os.path.join(OUT, f"results_gstools_{variant}"))
            clean_tests_work()

    if not args.skip_gstlearn and not args.skip_gstools:
        if not info["gstlearn"]:
            log("\ngstlearn not installed: skipping the gstlearn comparison (py -m pip install gstlearn)")
        else:
            clean_tests_work()
            shutil.rmtree(os.path.join(TESTS, "results"), ignore_errors=True)
            run_stage("gstlearn", "compare_gstlearn (simuSpectral, exp+gauss)", [sys.executable, os.path.join(TESTS, "compare_gstlearn.py")],
                      os.path.join(OUT, "compare_gstlearn.log"))
            collect_results(os.path.join(OUT, "results_gstlearn"))
            clean_tests_work()

    if not args.skip_figures:
        cmd = [sys.executable, os.path.join(DOCS, "make_figures.py")] + ([] if args.quick else ["--full"])
        if not have_gs:
            cmd.append("--no-gstools")
        run_stage("figures", "make_figures" + ("" if args.quick else " --full"), cmd, os.path.join(OUT, "make_figures.log"))
        sysjson = os.path.join(DOCS, "figures", "system_info.json")
        if os.path.exists(sysjson):
            shutil.copy(sysjson, os.path.join(OUT, "system_info_make_figures.json"))
        figs = os.path.join(DOCS, "figures")
        if os.path.isdir(figs):
            shutil.copytree(figs, os.path.join(OUT, "figures"), dirs_exist_ok=True,
                            ignore=shutil.ignore_patterns("_work"))

    # ---- summary
    total = time.time() - t00
    nfail = sum(1 for r in RESULTS if r[2] == "FAIL")
    md = ["# SPECSIM3D validation summary", "", f"Run: {info['date']}   total time: {total / 60:.1f} min", "",
          "## Machine (for the paper, Section 4.6)", "",
          f"- CPU: {info['cpu']} ({info['cores_physical']} physical / {info['cores_logical']} logical cores)",
          f"- RAM: {info['ram_gb']} GB", f"- OS: {info['os']}", f"- Fortran compiler: {info['fortran_compiler']}",
          f"- Executable: {info['executable']} built {info['executable_built']}",
          f"- Python {info['python']}, numpy {info['numpy']}, scipy {info['scipy']}, matplotlib {info['matplotlib']}, GSTools {info['gstools']}, gstlearn {info['gstlearn']}",
          f"- OpenMP threads: {info['omp_num_threads_env']}", "", "## Results", "",
          "| stage | test | status | time (s) |", "|---|---|---|---|"]
    md += [f"| {s} | {n} | {st} | {t} |" for s, n, st, t in RESULTS]
    md += ["", f"**{len(RESULTS) - nfail} passed, {nfail} failed.**", "",
           "Logs: *.log; test outputs: results_*; figures and tables: figures/ (copy of docs/figures)."]
    open(os.path.join(OUT, "summary.md"), "w", encoding="utf-8").write("\n".join(md) + "\n")
    json.dump({"machine": info, "results": RESULTS, "total_s": round(total)}, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
    log("\n" + "\n".join(md))
    log(f"\nAll results in: {OUT}")
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main())
