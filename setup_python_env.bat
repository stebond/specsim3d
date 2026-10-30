@echo off
rem ============================================================================
rem  SPECSIM3D - prepares a Windows machine with NO Python installed:
rem    1. looks for a usable Python 3.9 or newer; if none, downloads and installs the
rem       newest tested Python (up to 3.14) from python.org (per-user install, no admin rights)
rem    2. installs/updates the packages needed by all tests, figures and benchmarks:
rem         numpy scipy matplotlib gstools psutil   (+ optional: gstlearn)
rem    3. checks the imports, saves the interpreter in specsim_python.txt (used by
rem       run_all_tests.bat and run_full_validation.bat) and writes python_env.txt
rem  If the newest Python has no binary wheels yet for a package (typical for
rem  brand-new Python versions), the previous minor versions are tried automatically.
rem
rem  Usage:   setup_python_env.bat            (reuse an existing Python if it works)
rem           setup_python_env.bat /force     (always install a fresh Python)
rem           add /nopause to avoid the final "press a key"
rem  Needs an internet connection (python.org and pypi.org).
rem  Not covered here: the Fortran compiler (Intel oneAPI ifx), see build_specsim3d.bat
rem ============================================================================
setlocal EnableExtensions EnableDelayedExpansion
cd /d "%~dp0"
set "PKGS=numpy scipy matplotlib gstools psutil"
set "IMPORTS=numpy,scipy,matplotlib,gstools,psutil"
set "FORCE=0"
set "NOPAUSE=0"
for %%A in (%*) do (
  if /i "%%A"=="/force" set "FORCE=1"
  if /i "%%A"=="--force" set "FORCE=1"
  if /i "%%A"=="/nopause" set "NOPAUSE=1"
)
set "PYEXE="
set "PKGOK=0"
set "VERFILE=%TEMP%\specsim_py_versions.txt"
set "INSTALLER=%TEMP%\specsim_python_installer.exe"

echo.
echo ==== SPECSIM3D: Python environment setup ====

rem ---- 1. existing Python? -----------------------------------------------------
if "%FORCE%"=="0" call :probe py -3
if "%FORCE%"=="0" if not defined PYEXE call :probe python
if defined PYEXE (
  echo Found a usable Python: !PYEXE!
  call :pkgs
  if "!PKGOK!"=="1" goto :success
  echo The existing Python could not get all the packages: installing a fresh Python.
)

rem ---- 2. download + install latest stable Python ----------------------------------
echo.
echo Looking for the latest Python on python.org ...
del "%VERFILE%" 2>nul
set "SS_MODE=list"
set "SS_OUT=%VERFILE%"
call :runps
if errorlevel 1 goto :fail_net
if not exist "%VERFILE%" goto :fail_net

for /f "usebackq delims=" %%V in ("%VERFILE%") do (
  call :installpy %%V
  if "!PKGOK!"=="1" goto :success
)
goto :fail

rem ---- success -------------------------------------------------------------------
:success
> "%~dp0specsim_python.txt" echo "!PYEXE!"
"!PYEXE!" -c "import sys,numpy,scipy,matplotlib,gstools,psutil;print('Python',sys.version.split()[0]);print('numpy',numpy.__version__);print('scipy',scipy.__version__);print('matplotlib',matplotlib.__version__);print('gstools',gstools.__version__);print('psutil',psutil.__version__)" > "%~dp0python_env.txt"
"!PYEXE!" -c "import gstlearn;print('gstlearn',gstlearn.__version__)" >> "%~dp0python_env.txt" 2>nul
echo.
echo ==== READY ====
echo Interpreter: !PYEXE!
type "%~dp0python_env.txt"
where ifx >nul 2>&1
if errorlevel 1 (
  echo.
  echo NOTE: ifx not found in this prompt. Open the "Intel oneAPI command prompt" before
  echo       build_specsim3d.bat / run_full_validation.bat so that the compiler is recorded.
)
echo.
echo Next steps:   build_specsim3d.bat     then     run_full_validation.bat
if "%NOPAUSE%"=="0" pause
exit /b 0

:fail_net
echo.
echo ERROR: could not read the list of Python versions from python.org.
echo        Check the internet connection / proxy and run this script again.
if "%NOPAUSE%"=="0" pause
exit /b 1

:fail
echo.
echo ERROR: no Python version could install all the packages (%PKGS%).
echo        Read the pip messages above; install Python 3.12 or 3.13 manually from
echo        python.org (tick "Add python.exe to PATH") and run this script again.
if "%NOPAUSE%"=="0" pause
exit /b 1

rem ---- subroutines ---------------------------------------------------------------

:probe
rem  probe [command]: sets PYEXE if the command is a real Python 3.9 or newer
set "P1="
set "P2="
set "TMPV=%TEMP%\specsim_py_probe.txt"
del "%TMPV%" 2>nul
%* -c "import sys;print(sys.executable);print(sys.version_info[0]*100+sys.version_info[1])" >"%TMPV%" 2>nul
if errorlevel 1 exit /b 1
for /f "usebackq delims=" %%L in ("%TMPV%") do (
  if not defined P1 (set "P1=%%L") else set "P2=%%L"
)
if not defined P2 exit /b 1
if %P2% LSS 309 exit /b 1
set "PYEXE=%P1%"
exit /b 0

:pkgs
rem  installs or updates the packages (binary wheels only: fails fast if none)
set "PKGOK=0"
"!PYEXE!" -m ensurepip --upgrade >nul 2>&1
"!PYEXE!" -m pip install --upgrade pip
echo Installing packages: %PKGS%
"!PYEXE!" -m pip install --upgrade --only-binary=:all: %PKGS%
"!PYEXE!" -c "import %IMPORTS%" >nul 2>&1 && set "PKGOK=1"
if "!PKGOK!"=="1" (
  echo Installing the optional package gstlearn ^(comparison with gstlearn; failure is not an error^)
  "!PYEXE!" -m pip install --upgrade --only-binary=:all: gstlearn
)
exit /b 0

:installpy
rem  installpy [x.y.z]: downloads, installs (per user, silent) and tests that version
set "PKGOK=0"
set "PYEXE="
set "VER=%1"
echo.
echo ---- Python %VER% ----
del "%INSTALLER%" 2>nul
set "SS_MODE=get"
set "SS_VER=%VER%"
set "SS_OUT=%INSTALLER%"
call :runps
if errorlevel 1 (echo Download of Python %VER% failed. & exit /b 1)
if not exist "%INSTALLER%" (echo Download of Python %VER% failed. & exit /b 1)
echo Installing Python %VER% (silent, current user only) ...
start /wait "" "%INSTALLER%" /quiet InstallAllUsers=0 PrependPath=1 Include_launcher=1 Include_pip=1 Include_test=0 Include_doc=0 Include_tcltk=1 Shortcuts=0
for /f "tokens=1,2 delims=." %%a in ("%VER%") do set "MM=%%a%%b"
for /d %%D in ("%LOCALAPPDATA%\Programs\Python\Python!MM!*") do (
  if exist "%%D\python.exe" set "PYEXE=%%D\python.exe"
)
if not defined PYEXE (echo Installed Python %VER% not found under %LOCALAPPDATA%\Programs\Python & exit /b 1)
call :probe "!PYEXE!"
if errorlevel 1 (echo Python %VER% does not start. & exit /b 1)
call :pkgs
if "!PKGOK!"=="0" echo Packages could not be installed on Python %VER%: trying the previous version.
exit /b 0

:runps
rem  runs the PowerShell code stored at the end of this file (after the marker line)
powershell -NoProfile -ExecutionPolicy Bypass -Command "$s=[IO.File]::ReadAllText('%~f0'); iex $s.Substring($s.IndexOf('#PS-' + 'BEGIN'))"
exit /b %errorlevel%

rem ---- end of the batch part (nothing below is executed by cmd) -----------------------
exit /b 0
#PS-BEGIN
$ErrorActionPreference = 'Stop'
try {
  [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
  try { [Net.WebRequest]::DefaultWebProxy.Credentials = [Net.CredentialCache]::DefaultCredentials } catch {}
  $arch = if ($env:PROCESSOR_ARCHITECTURE -eq 'ARM64') { 'arm64' } else { 'amd64' }
  $base = 'https://www.python.org/ftp/python'
  if ($env:SS_MODE -eq 'list') {
    # newest stable patch release (3.9+) that has a Windows installer, for the 3 newest minor versions
    $html = (Invoke-WebRequest -UseBasicParsing "$base/").Content
    $all = [regex]::Matches($html, 'href="3\.(\d+)\.(\d+)/"') | ForEach-Object {
      [pscustomobject]@{ Minor = [int]$_.Groups[1].Value; Patch = [int]$_.Groups[2].Value } } | Where-Object { $_.Minor -ge 9 }
    $out = @()
    # prefer the newest minor version up to 3.14 (binary wheels of numpy, scipy, gstools exist and are tested),
    # then the newer ones, which are tried last because their wheels often do not exist yet
    $mins = $all | Select-Object -ExpandProperty Minor -Unique
    $pref = @($mins | Where-Object { $_ -le 14 } | Sort-Object -Descending) + @($mins | Where-Object { $_ -gt 14 } | Sort-Object -Descending)
    foreach ($m in $pref) {
      if ($out.Count -ge 3) { break }
      foreach ($p in ($all | Where-Object { $_.Minor -eq $m } | Sort-Object Patch -Descending)) {
        $v = "3.$m.$($p.Patch)"
        try {
          $r = Invoke-WebRequest -UseBasicParsing -Method Head "$base/$v/python-$v-$arch.exe"
          if ($r.StatusCode -eq 200) { $out += $v; break }
        } catch {}
      }
    }
    if ($out.Count -eq 0) { throw 'no Python installer found' }
    $out | Set-Content -Encoding ASCII $env:SS_OUT
  }
  elseif ($env:SS_MODE -eq 'get') {
    $v = $env:SS_VER
    Invoke-WebRequest -UseBasicParsing -Uri "$base/$v/python-$v-$arch.exe" -OutFile $env:SS_OUT
  }
  exit 0
} catch {
  Write-Host ("PowerShell error: " + $_.Exception.Message)
  exit 1
}
