<#
.SYNOPSIS
    Fetch the public input data for G4_Promoter_Pipeline from NCBI GEO.

.DESCRIPTION
    Downloads, into the repository-local layout expected by config/config.yaml:

      data/bigwig/   28 CUT&Tag coverage tracks (GSE269084; ~3.4 GB total)
      data/rnaseq/   GSE269081_count_table.tsv  (RNA-seq counts, ~1.9 MB unpacked)

    With -Refs it also fetches the two reference annotations the workflow would
    otherwise download itself on first run:

      data/gencode.vM25.annotation.gtf.gz   (~28 MB, GENCODE)
      data/rmsk_mm10.txt.gz                 (~142 MB, UCSC RepeatMasker)

    The bigWig list, expected byte sizes and sample metadata come from
    data/geo_manifest.tsv. Existing files whose size already matches the manifest
    are skipped, and partial downloads are resumed, so the script is safe to
    re-run after an interruption.

    All 28 bigWigs are needed for a complete run: the ERCC tracks feed the
    ERCC peak sets built by scripts/07_call_peaks.R, which `rule all` requires.
    Use -SkipErcc only if you also drop ERCCWT/ERCCKO from
    `peak_calling.genotypes` in config/config.yaml.

.PARAMETER Refs
    Also download the GENCODE vM25 GTF and the mm10 RepeatMasker table.

.PARAMETER SkipErcc
    Skip the 8 ERCC control tracks (~1.3 GB). See the note above.

.PARAMETER WhatIf
    List every URL and what would happen, without downloading anything.

.EXAMPLE
    .\scripts\download_data.ps1
    Download the 28 bigWigs and the RNA-seq count table.

.EXAMPLE
    .\scripts\download_data.ps1 -Refs -WhatIf
    Show the full download plan, including reference annotations.

.NOTES
    Already have the bigWigs elsewhere? Don't copy them — point the workflow at
    them with a local overlay instead:
      cp config/config.local.example.yaml config/config.local.yaml   # edit paths
      snakemake --cores 4 --configfile config/config.local.yaml
#>
[CmdletBinding(SupportsShouldProcess = $true)]
param(
    [switch]$Refs,
    [switch]$SkipErcc
)

$ErrorActionPreference = 'Stop'

# Repository root = parent of the directory holding this script.
$RepoRoot   = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$Manifest   = Join-Path $RepoRoot 'data\geo_manifest.tsv'
$BigwigDir  = Join-Path $RepoRoot 'data\bigwig'
$RnaseqDir  = Join-Path $RepoRoot 'data\rnaseq'
$DataDir    = Join-Path $RepoRoot 'data'

$GeoSamples = 'https://ftp.ncbi.nlm.nih.gov/geo/samples'
$CountsUrl  = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE269nnn/GSE269081/suppl/GSE269081_count_table.tsv.gz'
$GtfUrl     = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz'
$RmskUrl    = 'https://hgdownload.gi.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz'

$Curl = (Get-Command curl.exe -ErrorAction SilentlyContinue).Source
if (-not $Curl) {
    throw "curl.exe not found. It ships with Windows 10 1803+ and Windows 11; on older systems install curl or use scripts/download_data.sh under WSL."
}
if (-not (Test-Path $Manifest)) { throw "Manifest not found: $Manifest" }

function Get-GsmUrl {
    # GEO buckets samples by GSM accession with the last three digits masked:
    #   GSM8305989 -> .../samples/GSM8305nnn/GSM8305989/suppl/<filename>
    param([string]$Gsm, [string]$FileName)
    $bucket = $Gsm.Substring(0, $Gsm.Length - 3) + 'nnn'
    "$GeoSamples/$bucket/$Gsm/suppl/$FileName"
}

function Get-RemoteFile {
    <# Download $Url to $Dest, resuming partial transfers. Skips the download when
       $ExpectedBytes is given and the local file already has exactly that size. #>
    param([string]$Url, [string]$Dest, [long]$ExpectedBytes = 0)

    $name = Split-Path -Leaf $Dest
    if ($ExpectedBytes -gt 0 -and (Test-Path $Dest) -and (Get-Item $Dest).Length -eq $ExpectedBytes) {
        Write-Host ("  [skip] {0} (already complete)" -f $name)
        return $true
    }
    if (-not $PSCmdlet.ShouldProcess($Dest, "download $Url")) {
        Write-Host ("  [plan] {0}`n         <- {1}" -f $name, $Url)
        return $true
    }

    Write-Host ("  [get ] {0}" -f $name)
    & $Curl --fail --location --retry 3 --retry-delay 5 --continue-at - `
            --progress-bar --output $Dest $Url
    if ($LASTEXITCODE -ne 0) {
        Write-Warning ("  download failed (curl exit {0}): {1}" -f $LASTEXITCODE, $Url)
        return $false
    }
    if ($ExpectedBytes -gt 0) {
        $actual = (Get-Item $Dest).Length
        if ($actual -ne $ExpectedBytes) {
            Write-Warning ("  size mismatch for {0}: got {1} bytes, expected {2}" -f $name, $actual, $ExpectedBytes)
            return $false
        }
    }
    return $true
}

New-Item -ItemType Directory -Force -Path $BigwigDir, $RnaseqDir | Out-Null

# --- CUT&Tag bigWigs --------------------------------------------------------
$rows = Import-Csv -Path $Manifest -Delimiter "`t"
if ($SkipErcc) { $rows = $rows | Where-Object { $_.group -ne 'ercc' } }
$totalGb = [math]::Round((($rows | Measure-Object -Property bytes -Sum).Sum) / 1GB, 2)

Write-Host ""
Write-Host ("CUT&Tag bigWigs: {0} files, {1} GB  ->  data\bigwig\" -f $rows.Count, $totalGb)
$failed = @()
foreach ($r in $rows) {
    $url = Get-GsmUrl -Gsm $r.gsm -FileName $r.filename
    $ok  = Get-RemoteFile -Url $url -Dest (Join-Path $BigwigDir $r.filename) -ExpectedBytes ([long]$r.bytes)
    if (-not $ok) { $failed += $r.filename }
}

# --- RNA-seq count table ----------------------------------------------------
Write-Host ""
Write-Host "RNA-seq count table (GSE269081)  ->  data\rnaseq\"
$gz  = Join-Path $RnaseqDir 'GSE269081_count_table.tsv.gz'
$tsv = Join-Path $RnaseqDir 'GSE269081_count_table.tsv'
if (Test-Path $tsv) {
    Write-Host "  [skip] GSE269081_count_table.tsv (already present)"
} elseif (Get-RemoteFile -Url $CountsUrl -Dest $gz) {
    if ($PSCmdlet.ShouldProcess($tsv, 'gunzip')) {
        # .NET GZipStream: no external gzip needed on Windows.
        $in  = [System.IO.File]::OpenRead($gz)
        $out = [System.IO.File]::Create($tsv)
        try {
            $gzs = New-Object System.IO.Compression.GZipStream($in, [System.IO.Compression.CompressionMode]::Decompress)
            $gzs.CopyTo($out)
            $gzs.Dispose()
        } finally { $out.Dispose(); $in.Dispose() }
        Remove-Item $gz
        Write-Host "  [ok  ] unpacked GSE269081_count_table.tsv"
    }
} else { $failed += 'GSE269081_count_table.tsv.gz' }

# --- Optional reference annotations ----------------------------------------
if ($Refs) {
    Write-Host ""
    Write-Host "Reference annotations  ->  data\"
    if (-not (Get-RemoteFile -Url $GtfUrl  -Dest (Join-Path $DataDir 'gencode.vM25.annotation.gtf.gz'))) { $failed += 'gencode.vM25.annotation.gtf.gz' }
    if (-not (Get-RemoteFile -Url $RmskUrl -Dest (Join-Path $DataDir 'rmsk_mm10.txt.gz')))               { $failed += 'rmsk_mm10.txt.gz' }
}

# --- Upstream model repository ---------------------------------------------
Write-Host ""
Write-Host "G4ShapePredictor (topology model)  ->  external\"
$g4sp = Join-Path $RepoRoot 'external\G4ShapePredictor'
if (Test-Path (Join-Path $g4sp '.git')) {
    Write-Host "  [skip] external\G4ShapePredictor (already cloned)"
} elseif (Get-Command git -ErrorAction SilentlyContinue) {
    if ($PSCmdlet.ShouldProcess($g4sp, 'git clone https://github.com/donn-liew/G4ShapePredictor')) {
        New-Item -ItemType Directory -Force -Path (Join-Path $RepoRoot 'external') | Out-Null
        git clone --depth 1 https://github.com/donn-liew/G4ShapePredictor $g4sp
    }
} else {
    Write-Warning "  git not found - clone https://github.com/donn-liew/G4ShapePredictor into external\ manually"
}

Write-Host ""
if ($failed.Count -gt 0) {
    Write-Warning ("{0} item(s) did not complete: {1}" -f $failed.Count, ($failed -join ', '))
    Write-Warning "Re-run this script to resume; downloads continue where they stopped."
    exit 1
}
Write-Host "Done. Next: conda activate G4; snakemake -n" -ForegroundColor Green
