Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$repoRoot = [System.IO.Path]::GetFullPath((Split-Path -Parent $PSScriptRoot))
$releaseRoot = [System.IO.Path]::GetFullPath((Join-Path $repoRoot 'output/release'))
$stageRoot = [System.IO.Path]::GetFullPath((Join-Path $releaseRoot 'tsp_reproducibility_artifact'))
$zipPath = [System.IO.Path]::GetFullPath((Join-Path $releaseRoot 'tsp_reproducibility_artifact.zip'))
$requiredPrefix = $releaseRoot.TrimEnd([System.IO.Path]::DirectorySeparatorChar) + `
    [System.IO.Path]::DirectorySeparatorChar

if (-not $stageRoot.StartsWith($requiredPrefix, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to stage outside the release directory: $stageRoot"
}
if (-not $zipPath.StartsWith($requiredPrefix, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to write an archive outside the release directory: $zipPath"
}

New-Item -ItemType Directory -Force -Path $releaseRoot | Out-Null
if (Test-Path -LiteralPath $stageRoot) {
    Remove-Item -LiteralPath $stageRoot -Recurse -Force
}
if (Test-Path -LiteralPath $zipPath) {
    Remove-Item -LiteralPath $zipPath -Force
}
New-Item -ItemType Directory -Force -Path $stageRoot | Out-Null

function Copy-ArtifactPath {
    param([Parameter(Mandatory = $true)][string]$RelativePath)

    $source = Join-Path $repoRoot $RelativePath
    if (-not (Test-Path -LiteralPath $source)) {
        throw "Missing artifact input: $RelativePath"
    }
    $destination = Join-Path $stageRoot $RelativePath
    $destinationParent = Split-Path -Parent $destination
    New-Item -ItemType Directory -Force -Path $destinationParent | Out-Null
    Copy-Item -LiteralPath $source -Destination $destination -Recurse -Force
}

$artifactPaths = @(
    '.gitattributes',
    'LICENSE',
    'README.md',
    'pyproject.toml',
    'REPRODUCIBILITY.md',
    'hidden_completion',
    'datasets/raw',
    'experiments/__init__.py',
    'experiments/test_exact_scenario_solver.py',
    'experiments/test_smoother_experiments.py',
    'experiments/requirements-smoother.txt',
    'experiments/requirements-repro-cpu.txt',
    'experiments/smoother',
    'experiments/configs',
    'experiments/results/p3_p4_smoke.json',
    'experiments/results/p4_core_evidence.json',
    'experiments/results/p4_solver_crosscheck.json',
    'experiments/results/p4_parameter_sweep.json',
    'experiments/results/p4_external_scale.json',
    'experiments/results/smoother_gpu_scale.json',
    'experiments/results/smoother_gpu_scale_recheck1.json',
    'experiments/results/smoother_gpu_scale_recheck2.json',
    'experiments/scalable',
    'paper/figures/tikz',
    'paper/tsp/main.tex',
    'paper/tsp/supplement.tex',
    'paper/tsp/references.bib',
    'paper/tsp/README.md',
    'paper/tsp/tables',
    'scripts/build_tsp.ps1',
    'scripts/package_tsp_artifact.ps1'
)

foreach ($relativePath in $artifactPaths) {
    Copy-ArtifactPath -RelativePath $relativePath
}

$generatedPatterns = @('__pycache__', 'build')
foreach ($directoryName in $generatedPatterns) {
    Get-ChildItem -LiteralPath $stageRoot -Directory -Recurse -Force |
        Where-Object { $_.Name -eq $directoryName } |
        Sort-Object FullName -Descending |
        ForEach-Object { Remove-Item -LiteralPath $_.FullName -Recurse -Force }
}
Get-ChildItem -LiteralPath $stageRoot -File -Recurse -Force |
    Where-Object { $_.Extension -in @('.aux', '.log', '.out', '.fls', '.fdb_latexmk', '.pyc') } |
    ForEach-Object { Remove-Item -LiteralPath $_.FullName -Force }

$manifestPath = Join-Path $stageRoot 'MANIFEST.sha256'
$manifestLines = Get-ChildItem -LiteralPath $stageRoot -File -Recurse |
    Where-Object { $_.FullName -ne $manifestPath } |
    Sort-Object FullName |
    ForEach-Object {
        $relative = $_.FullName.Substring($stageRoot.Length + 1).Replace('\', '/')
        $hash = (Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
        "$hash  $relative"
    }
$manifestLines | Set-Content -LiteralPath $manifestPath -Encoding utf8

Compress-Archive -LiteralPath $stageRoot -DestinationPath $zipPath -CompressionLevel Optimal
$zipHash = (Get-FileHash -LiteralPath $zipPath -Algorithm SHA256).Hash
$fileCount = (Get-ChildItem -LiteralPath $stageRoot -File -Recurse).Count
$zipSize = (Get-Item -LiteralPath $zipPath).Length

Write-Host "Packaged $fileCount files"
Write-Host "Archive: $zipPath"
Write-Host "Bytes: $zipSize"
Write-Host "SHA-256: $zipHash"
