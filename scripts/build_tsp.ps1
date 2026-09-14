param(
    [switch]$SkipTests,
    [switch]$SkipData
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$repoRoot = Split-Path -Parent $PSScriptRoot
$figureDirectory = Join-Path $repoRoot 'paper/figures/tikz'
$tspDirectory = Join-Path $repoRoot 'paper/tsp'
$tspFigureDirectory = Join-Path $tspDirectory 'figures'
$tspBuildDirectory = Join-Path $tspDirectory 'build'
$outputDirectory = Join-Path $repoRoot 'output/pdf'
$env:SOURCE_DATE_EPOCH = '1789344000'
$env:FORCE_SOURCE_DATE = '1'

$figures = @(
    'theorem_structure',
    'budget_pareto',
    'scenario_scaling',
    'core_design_benchmark',
    'sampling_budget_curve',
    'path_tightness',
    'local_method_tradeoff',
    'external_stress',
    'cpu_gpu_crossover'
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
            '--core-evidence', 'experiments/results/p4_core_evidence.json',
            '--solver-crosscheck', 'experiments/results/p4_solver_crosscheck.json',
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

    New-Item -ItemType Directory -Force -Path $tspFigureDirectory | Out-Null
    New-Item -ItemType Directory -Force -Path $tspBuildDirectory | Out-Null
    New-Item -ItemType Directory -Force -Path $outputDirectory | Out-Null
    foreach ($figure in $figures) {
        Copy-Item -LiteralPath (Join-Path $figureDirectory "$figure.pdf") `
            -Destination (Join-Path $tspFigureDirectory "$figure.pdf") -Force
    }

    Push-Location $tspDirectory
    try {
        Invoke-Checked -Program 'latexmk' -ArgumentList @(
            '-pdf', '-interaction=nonstopmode', '-halt-on-error',
            '-outdir=build', 'main.tex'
        )
        Invoke-Checked -Program 'latexmk' -ArgumentList @(
            '-pdf', '-interaction=nonstopmode', '-halt-on-error',
            '-outdir=build', 'supplement.tex'
        )
    }
    finally {
        Pop-Location
    }

    Copy-Item -LiteralPath (Join-Path $tspBuildDirectory 'main.pdf') `
        -Destination (Join-Path $outputDirectory 'tsp_manuscript.pdf') -Force
    Copy-Item -LiteralPath (Join-Path $tspBuildDirectory 'supplement.pdf') `
        -Destination (Join-Path $outputDirectory 'tsp_supplement.pdf') -Force

    Write-Host "Built $outputDirectory/tsp_manuscript.pdf"
    Write-Host "Built $outputDirectory/tsp_supplement.pdf"
}
finally {
    Pop-Location
}
