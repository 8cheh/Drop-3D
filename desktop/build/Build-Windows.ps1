# Build the Drop-3D desktop application for Windows.
#
#   powershell -ExecutionPolicy Bypass -File desktop\build\Build-Windows.ps1
#
# Produces desktop\dist\Drop-3D\ (a folder you can double-click into) and
# desktop\dist\Drop-3D-<version>-win64.zip.
#
# Two deliberate choices, both learned the hard way:
#
#   * This file is pure ASCII and takes no policy from $ErrorActionPreference.
#     Windows PowerShell 5.1 reads .ps1 as ANSI unless it carries a BOM, so any
#     non-ASCII path written here would be mangled; and treating native stderr as
#     terminating turns PyInstaller's own deprecation notice into a build failure.
#     All paths are derived from the script location instead of being written out.
#
#   * The build virtual environment is created WITHOUT --system-site-packages.
#     Sharing the machine's site-packages makes PyInstaller bundle everything it
#     can see: the first attempt here produced 932 MB including torch and OpenCV.
#
# Extra arguments are passed through to the application (for example --acrylic).

param(
    [switch]$SkipTests,
    [switch]$NoZip
)

$ErrorActionPreference = 'Continue'

$BuildDir = $PSScriptRoot                      # desktop/build
$Desktop  = Split-Path $BuildDir -Parent       # desktop
$Repo     = Split-Path $Desktop -Parent        # repository root
$Venv     = Join-Path $Desktop '.build-venv'
$Python   = Join-Path $Venv 'Scripts\python.exe'
$Dist     = Join-Path $Desktop 'dist'
$Work     = Join-Path $Desktop 'build\work'
$Req      = Join-Path $Desktop 'requirements-desktop.txt'
$Spec     = Join-Path $BuildDir 'Drop3D.spec'

Write-Host '=== Drop-3D desktop build ===' -ForegroundColor Cyan
Write-Host "repository : $Repo"
Write-Host "output     : $Dist"
Write-Host ''

# ---------------------------------------------------------------- environment
if (-not (Test-Path $Python)) {
    Write-Host '[1/5] creating a clean build virtual environment'
    python -m venv $Venv
    if ($LASTEXITCODE -ne 0) { Write-Error 'failed to create the virtual environment'; exit 1 }
} else {
    Write-Host '[1/5] reusing the existing build virtual environment'
}

Write-Host '[2/5] installing pinned dependencies'
& $Python -m pip install --quiet --upgrade pip
& $Python -m pip install --quiet -r $Req
if ($LASTEXITCODE -ne 0) { Write-Error 'dependency install failed'; exit 1 }

# Make sure the frozen application is not built from a stale copy of the library.
Write-Host '       installing the repository itself'
& $Python -m pip install --quiet --no-deps --force-reinstall $Repo
if ($LASTEXITCODE -ne 0) {
    # A plain install fails when the project metadata lives in a checkout with no
    # build backend available; falling back to a source path keeps the build going.
    Write-Host '       falling back to installing from src/'
    & $Python -m pip install --quiet --no-deps --force-reinstall (Join-Path $Repo 'src')
}

Write-Host '       frozen versions:'
& $Python -m pip list --format=freeze | Select-String -Pattern '^(numpy|scipy|pillow|pywebview|pythonnet|clr-loader|pyinstaller)='

# ------------------------------------------------------------------- tests
if ($SkipTests) {
    Write-Host '[3/5] tests skipped (-SkipTests)'
} else {
    Write-Host '[3/5] running the library test suite before freezing'
    Push-Location $Repo
    & $Python -m pip install --quiet pytest
    & $Python -m pytest -q
    $testExit = $LASTEXITCODE
    Pop-Location
    if ($testExit -ne 0) { Write-Error 'tests failed; refusing to build'; exit 1 }
}

# -------------------------------------------------------------------- freeze
Write-Host '[4/5] freezing'
foreach ($p in @($Dist, $Work)) {
    if (Test-Path $p) { Remove-Item -Recurse -Force $p }
}
# No --specpath here: it is a makespec option and PyInstaller rejects it when a
# .spec file is given (the spec already knows where it lives).
& $Python -m PyInstaller --noconfirm --clean `
    --distpath $Dist --workpath $Work $Spec
if ($LASTEXITCODE -ne 0) { Write-Error 'PyInstaller failed'; exit 1 }

$AppDir = Join-Path $Dist 'Drop-3D'
$Exe = Join-Path $AppDir 'Drop-3D.exe'
if (-not (Test-Path $Exe)) { Write-Error "expected $Exe to exist"; exit 1 }

$files = Get-ChildItem $AppDir -Recurse -File
$sizeMb = [math]::Round((($files | Measure-Object -Property Length -Sum).Sum / 1MB), 1)
Write-Host ("       bundle: {0} MB across {1} files" -f $sizeMb, $files.Count)

# ----------------------------------------------------------------- zip
if ($NoZip) {
    Write-Host '[5/5] zip skipped (-NoZip)'
} else {
    Write-Host '[5/5] packaging the zip'
    $version = (& $Python -c "import drop3d; print(drop3d.__version__)").Trim()
    if (-not $version) { $version = '0.0.0' }
    $zip = Join-Path $Dist ("Drop-3D-{0}-win64.zip" -f $version)
    if (Test-Path $zip) { Remove-Item -Force $zip }
    Compress-Archive -Path $AppDir -DestinationPath $zip -CompressionLevel Optimal
    $zipMb = [math]::Round((Get-Item $zip).Length / 1MB, 1)
    Write-Host ("       {0}  ({1} MB)" -f $zip, $zipMb)
}

Write-Host ''
Write-Host 'Done. Run the application with:' -ForegroundColor Green
Write-Host "  $Exe"
Write-Host ''
Write-Host 'It needs the Microsoft Edge WebView2 runtime, which ships with Windows 11'
Write-Host 'and most Windows 10 installations.'
