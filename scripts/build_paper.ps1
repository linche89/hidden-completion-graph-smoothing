param(
    [switch]$SkipTests,
    [switch]$SkipData
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$repoRoot = Split-Path -Parent $PSScriptRoot
$figureDirectory = Join-Path $repoRoot 'paper/figures/tikz'
$paperDirectory = Join-Path $repoRoot 'paper'
$paperBuildDirectory = Join-Path $paperDirectory 'build'
# Freeze generated-PDF metadata at the paper baseline date so rebuilding an
# unchanged TikZ source does not dirty Git with timestamp-only differences.
$env:SOURCE_DATE_EPOCH = '1789344000'
$env:FORCE_SOURCE_DATE = '1'
$figures = @(
    'theorem_structure',
    'budget_pareto',
    'certificate_sampling_failure',
    'local_method_tradeoff',
    'scenario_scaling',
    'cpu_gpu_crossover',
    'external_stress'
)

function Invoke-Checked {
    param(
        [Parameter(Mandatory = $true)][string]$Program,
        [Parameter(Mandatory = $true)][string[]]$ArgumentList
    )
    & $Program @ArgumentList
    if ($LASTEXITCODE -ne 0) {
        throw "$Program failed with exit code $LASTEXITCODE"
    }
}

Push-Location $repoRoot
try {
    if (-not $SkipTests) {
        Invoke-Checked -Program 'python' -ArgumentList @(
            '-m', 'unittest', 'discover', '-s', 'experiments',
            '-p', 'test_*.py', '-v'
        )
    }

    if (-not $SkipData) {
        Invoke-Checked -Program 'python' -ArgumentList @(
            '-m', 'experiments.smoother.paper_data',
            '--p3-p4', 'experiments/results/p3_p4_smoke.json',
            '--external', 'experiments/results/p4_external_scale.json',
            '--gpu', 'experiments/results/smoother_gpu_scale.json',
            'experiments/results/smoother_gpu_scale_recheck1.json',
            'experiments/results/smoother_gpu_scale_recheck2.json'
        )
    }

    Push-Location $figureDirectory
    try {
        foreach ($figure in $figures) {
            Invoke-Checked -Program 'xelatex' -ArgumentList @(
                '-interaction=nonstopmode', '-halt-on-error', "$figure.tex"
            )
            Invoke-Checked -Program 'pdftoppm' -ArgumentList @(
                '-png', '-r', '200', '-singlefile', "$figure.pdf", $figure
            )
        }
    }
    finally {
        Pop-Location
    }

    New-Item -ItemType Directory -Force -Path $paperBuildDirectory | Out-Null
    Push-Location $paperDirectory
    try {
        Invoke-Checked -Program 'latexmk' -ArgumentList @(
            '-pdf', '-interaction=nonstopmode', '-halt-on-error',
            '-outdir=build', 'main.tex'
        )
    }
    finally {
        Pop-Location
    }

    Write-Host "Built $paperBuildDirectory/main.pdf"
}
finally {
    Pop-Location
}
