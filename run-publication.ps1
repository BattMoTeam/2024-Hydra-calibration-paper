param(
    [switch]$IncludeBpx,
    [string]$MatlabExe = "",
    [string]$PythonExe = ""
)

$ErrorActionPreference = "Stop"

$repoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $repoRoot

function Resolve-MatlabExecutable {
    param([string]$Override)

    if ($Override) {
        return $Override
    }

    if ($env:MATLAB_EXE -and (Test-Path $env:MATLAB_EXE)) {
        return $env:MATLAB_EXE
    }

    $cmd = Get-Command matlab.exe -ErrorAction SilentlyContinue
    if ($cmd) {
        return $cmd.Source
    }

    $candidates = @(
        "C:\Program Files\MATLAB\R2025b\bin\matlab.exe",
        "C:\Program Files\MATLAB\R2024b\bin\matlab.exe",
        "C:\Program Files\MATLAB\R2024a\bin\matlab.exe",
        "C:\Program Files\MATLAB\R2023b\bin\matlab.exe"
    )

    foreach ($candidate in $candidates) {
        if (Test-Path $candidate) {
            return $candidate
        }
    }

    throw "Could not locate MATLAB. Pass -MatlabExe or set MATLAB_EXE."
}

function Resolve-PythonExecutable {
    param([string]$Override)

    if ($Override) {
        return $Override
    }

    $venvPython = Join-Path $repoRoot ".venv\Scripts\python.exe"
    if (Test-Path $venvPython) {
        return $venvPython
    }

    if ($env:PYTHON_EXE) {
        return $env:PYTHON_EXE
    }

    $cmd = Get-Command python.exe -ErrorAction SilentlyContinue
    if ($cmd) {
        return $cmd.Source
    }

    throw "Could not locate Python. Pass -PythonExe or create .venv."
}

function Invoke-Checked {
    param(
        [string]$FilePath,
        [string[]]$Arguments
    )

    & $FilePath @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "Command failed: $FilePath $($Arguments -join ' ')"
    }
}

$resolvedMatlab = Resolve-MatlabExecutable -Override $MatlabExe
$resolvedPython = Resolve-PythonExecutable -Override $PythonExe

Write-Host "Repository root: $repoRoot"
Write-Host "MATLAB: $resolvedMatlab"
Write-Host "Python: $resolvedPython"
Write-Host "Running primary BattMo publication workflow..."

$matlabBatch = "startup; run(fullfile('scripts','exportValidationReference.m')); run(fullfile('scripts','exportRateStudyReference.m')); run(fullfile('scripts','exportPublicationFigures.m'));"
Invoke-Checked -FilePath $resolvedMatlab -Arguments @("-batch", $matlabBatch)

Invoke-Checked -FilePath $resolvedPython -Arguments @("scripts\plot_battmo_validation.py")

if ($IncludeBpx) {
    Write-Host "Running optional BPX / PyBaMM FAIR interoperability workflow..."
    Invoke-Checked -FilePath $resolvedPython -Arguments @("scripts\export_bpx.py")
    Invoke-Checked -FilePath $resolvedPython -Arguments @("scripts\verify_bpx.py", "--output", "codex\figures\bpx_verification_summary.json")
    Invoke-Checked -FilePath $resolvedPython -Arguments @("scripts\compare_battmo_pybamm.py")
}

Write-Host ""
Write-Host "Primary BattMo outputs:"
Write-Host "  figures\battmo-validation-reference.json"
Write-Host "  figures\figure-12-cell-balancing-under-equilibrium-assumption.fig"
Write-Host "  figures\figure-12-cell-balancing-under-equilibrium-assumption.png"
Write-Host "  figures\figure-13-high-rate-calibration-at-2C.fig"
Write-Host "  figures\figure-13-high-rate-calibration-at-2C.png"
Write-Host "  figures\figure-14-experimental-voltages-and-p2d-results.fig"
Write-Host "  figures\figure-14-experimental-voltages-and-p2d-results.png"
Write-Host "  figures\supporting\..."
Write-Host "  figures\publication\INP5-70-120-H0B_graphite-lnmo_schmitt-2026_battmo-vs-experiment-summary.json"
Write-Host "  figures\publication\INP5-70-120-H0B_graphite-lnmo_schmitt-2026_battmo-vs-experiment.png"
Write-Host "  figures\rate-study\battmo-rate-study-reference.json"

if ($IncludeBpx) {
    Write-Host ""
    Write-Host "Optional BPX / PyBaMM outputs:"
    Write-Host "  parameters\INP5-70-120-H0B_graphite-lnmo_schmitt-2026_validation.bpx.json"
    Write-Host "  codex\figures\bpx_verification_summary.json"
    Write-Host "  codex\figures\battmo-vs-pybamm-bpx-summary.json"
    Write-Host "  codex\figures\battmo-vs-pybamm-bpx.png"
}
