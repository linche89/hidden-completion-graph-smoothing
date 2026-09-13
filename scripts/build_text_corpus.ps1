[CmdletBinding()]
param(
    [switch]$ForceExtract
)

$ErrorActionPreference = 'Stop'
$projectRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$literatureRoot = Join-Path $projectRoot 'literature'
$corpusRoot = Join-Path $projectRoot 'research\corpus'
$manifestPath = Join-Path $corpusRoot 'manifest.tsv'

New-Item -ItemType Directory -Path $corpusRoot -Force | Out-Null

$pdfs = Get-ChildItem -LiteralPath $literatureRoot -Recurse -File -Filter '*.pdf' | Sort-Object FullName
$manifestRows = [System.Collections.Generic.List[string]]::new()
$manifestRows.Add("category`tfile`tbytes`tpages`tsha256`ttext_file")

foreach ($pdf in $pdfs) {
    $relativePath = [System.IO.Path]::GetRelativePath($literatureRoot, $pdf.FullName)
    $category = Split-Path -Path $relativePath -Parent
    $textRelativePath = [System.IO.Path]::ChangeExtension($relativePath, '.txt')
    $textPath = Join-Path $corpusRoot $textRelativePath
    $textDirectory = Split-Path -Path $textPath -Parent
    New-Item -ItemType Directory -Path $textDirectory -Force | Out-Null

    if ($ForceExtract -or -not (Test-Path -LiteralPath $textPath) -or $pdf.LastWriteTimeUtc -gt (Get-Item -LiteralPath $textPath).LastWriteTimeUtc) {
        Write-Host "TEXT  $relativePath"
        & pdftotext.exe -layout -nopgbrk -enc UTF-8 $pdf.FullName $textPath
        if ($LASTEXITCODE -ne 0) {
            throw "pdftotext failed for $($pdf.FullName)"
        }
    } else {
        Write-Host "SKIP  $relativePath"
    }

    $pageLine = (& pdfinfo.exe $pdf.FullName | Select-String -Pattern '^Pages:\s+(\d+)$').Matches
    $pages = if ($pageLine.Count -gt 0) { $pageLine[0].Groups[1].Value } else { '' }
    $hash = (Get-FileHash -Algorithm SHA256 -LiteralPath $pdf.FullName).Hash.ToLowerInvariant()
    $manifestRows.Add("$category`t$($pdf.Name)`t$($pdf.Length)`t$pages`t$hash`t$textRelativePath")
}

[System.IO.File]::WriteAllLines($manifestPath, $manifestRows, [System.Text.UTF8Encoding]::new($false))
Write-Host ""
Write-Host "PDFs: $($pdfs.Count); manifest: $manifestPath"

