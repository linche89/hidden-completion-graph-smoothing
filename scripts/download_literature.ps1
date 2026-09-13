[CmdletBinding()]
param(
    [switch]$ForceDownload,
    [string[]]$Category,
    [ValidateRange(10, 600)]
    [int]$MaxSeconds = 75
)

$ErrorActionPreference = 'Stop'
$projectRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$libraryRoot = Join-Path $projectRoot 'literature'

$papers = @(
    # Direct lineage of the seed WLS/state-estimation paper.
    @{ Category='01_direct_lineage'; File='2013_tai_lin_fu_sun_distributed_state_estimation_acc.pdf'; Url='https://www.eng.newcastle.edu.au/~mf140/home/Papers/ACC_2013_2.pdf' },
    @{ Category='01_direct_lineage'; File='2015_marelli_fu_distributed_wls_automatica.pdf'; Url='https://www.eng.newcastle.edu.au/~mf140/home/Papers/Automatica2014_2.pdf' },
    @{ Category='01_direct_lineage'; File='2015_marelli_ninness_fu_power_networks_sysid.pdf'; Url='https://www.eng.newcastle.edu.au/~mf140/home/Papers/SYSID_2015.pdf' },
    @{ Category='01_direct_lineage'; File='2018_sui_marelli_fu_lu_loopy_wls_accuracy.pdf'; Url='https://arxiv.org/pdf/1806.09104' },
    @{ Category='01_direct_lineage'; File='2020_yang_zhang_fu_wls_edge_measurements.pdf'; Url='https://arxiv.org/pdf/2002.11221' },
    @{ Category='01_direct_lineage'; File='2020_marelli_et_al_gaussian_bp_static_estimator.pdf'; Url='https://arxiv.org/pdf/2004.01969' },
    @{ Category='01_direct_lineage'; File='2012_pasqualetti_carli_bullo_iterative_projections.pdf'; Url='https://arxiv.org/pdf/1103.0579' },
    @{ Category='01_direct_lineage'; File='2018_wu_fu_xu_lu_distributed_kalman_loop_removal.pdf'; Url='https://www.eng.newcastle.edu.au/~mf140/home/Papers/Automatica_2018_3.pdf' },
    @{ Category='01_direct_lineage'; File='2023_cai_zhang_fu_belief_propagation_computations.pdf'; Url='https://www.sciopen.com/local/article_pdf/10.1016/j.jai.2023.06.003.pdf' },

    # Gaussian graphical models, correlation decay, and localized inference.
    @{ Category='02_graphical_models_local_inference'; File='2001_weiss_freeman_gaussian_bp_arbitrary_topology.pdf'; Url='https://proceedings.neurips.cc/paper_files/paper/1999/file/10c272d06794d3e5785d5e7c5356e9ff-Paper.pdf' },
    @{ Category='02_graphical_models_local_inference'; File='2006_malioutov_johnson_willsky_walk_sums_gaussian_bp.pdf'; Url='https://www.jmlr.org/papers/volume7/malioutov06a/malioutov06a.pdf' },
    @{ Category='02_graphical_models_local_inference'; File='2016_giscard_et_al_path_sums_exact_gaussian_inference.pdf'; Url='https://www.jmlr.org/papers/volume17/14-445/14-445.pdf' },
    @{ Category='02_graphical_models_local_inference'; File='2018_du_et_al_vector_gaussian_bp_convergence.pdf'; Url='https://www.jmlr.org/papers/volume18/16-556/16-556.pdf' },
    @{ Category='02_graphical_models_local_inference'; File='2018_feng_yin_local_distributed_sampling_counting.pdf'; Url='https://arxiv.org/pdf/1802.06686' },
    @{ Category='02_graphical_models_local_inference'; File='2018_chen_peng_liu_efficient_localized_inference.pdf'; Url='https://www.ijcai.org/proceedings/2018/0692.pdf' },
    @{ Category='02_graphical_models_local_inference'; File='2012_hatami_lovasz_szegedy_local_global_graph_limits.pdf'; Url='https://arxiv.org/pdf/1205.4356' },

    # Sparse inverse decay, graph filters, and operator locality.
    @{ Category='03_sparse_inverse_graph_filters'; File='1984_demko_moss_smith_decay_inverse_band_matrices.pdf'; Url='https://www.ams.org/mcom/1984-43-168/S0025-5718-1984-0758197-9/S0025-5718-1984-0758197-9.pdf' },
    @{ Category='03_sparse_inverse_graph_filters'; File='1990_jaffard_matrices_bien_localisees.pdf'; Url='https://www.numdam.org/article/AIHPC_1990__7_5_461_0.pdf' },
    @{ Category='03_sparse_inverse_graph_filters'; File='2007_benzi_razouk_decay_sparse_matrix_functions.pdf'; Url='https://etna.ricam.oeaw.ac.at/vol.28.2007-2008/pp16-39.dir/pp16-39.pdf' },
    @{ Category='03_sparse_inverse_graph_filters'; File='2011_shuman_et_al_distributed_chebyshev_approximation.pdf'; Url='https://arxiv.org/pdf/1111.5239' },
    @{ Category='03_sparse_inverse_graph_filters'; File='2015_segarra_marques_ribeiro_distributed_linear_network_operators.pdf'; Url='https://arxiv.org/pdf/1510.03947' },
    @{ Category='03_sparse_inverse_graph_filters'; File='2020_emirov_cheng_jiang_sun_polynomial_inverse_graph_filter.pdf'; Url='https://arxiv.org/pdf/2003.11152' },
    @{ Category='03_sparse_inverse_graph_filters'; File='2019_spakula_zhang_quasi_locality_property_a.pdf'; Url='https://arxiv.org/pdf/1809.00532' },
    @{ Category='03_sparse_inverse_graph_filters'; File='2020_li_spakula_zhang_measured_asymptotic_expanders_roe.pdf'; Url='https://arxiv.org/pdf/2010.10749' },
    @{ Category='03_sparse_inverse_graph_filters'; File='2023_ozawa_matrix_embeddings_quasi_local_not_uniform_roe.pdf'; Url='https://arxiv.org/pdf/2310.03677' },
    @{ Category='03_sparse_inverse_graph_filters'; File='2025_lu_wang_zhang_asymptotic_expansion_groupoids_roe.pdf'; Url='https://arxiv.org/pdf/2506.16772' },
    @{ Category='03_sparse_inverse_graph_filters'; File='2015_cheng_jiang_sun_spatially_distributed_sampling_reconstruction.pdf'; Url='https://arxiv.org/pdf/1511.08541' },
    @{ Category='03_sparse_inverse_graph_filters'; File='2014_motee_sun_sparsity_spatial_localization.pdf'; Url='https://arxiv.org/pdf/1402.4148' },
    @{ Category='03_sparse_inverse_graph_filters'; File='2026_li_zhang_zhu_k_theoretic_roe_quasi_local_projections.pdf'; Url='https://arxiv.org/pdf/2608.22439' },

    # Exact finite-time computation of linear maps and linear systems.
    @{ Category='04_distributed_linear_computation'; File='2015_costello_egerstedt_global_linear_to_local_edge_rules.pdf'; Url='https://repository.gatech.edu/bitstreams/7a449bad-91ec-4f0d-a1bf-bad9e0062941/download' },
    @{ Category='04_distributed_linear_computation'; File='2021_kar_pueschel_moura_finite_time_linear_transforms.pdf'; Url='https://arxiv.org/pdf/2104.01502' },
    @{ Category='04_distributed_linear_computation'; File='2018_yang_et_al_distributed_least_squares_solver.pdf'; Url='https://arxiv.org/pdf/1810.00156' },
    @{ Category='04_distributed_linear_computation'; File='2020_gade_liu_vaidya_finite_time_linear_equations.pdf'; Url='https://arxiv.org/pdf/2004.04680' },
    @{ Category='04_distributed_linear_computation'; File='2020_lei_peng_shi_anderson_finite_data_rates.pdf'; Url='https://openresearch-repository.anu.edu.au/server/api/core/bitstreams/23281126-0f32-4578-9eb4-6024ed079988/content' },

    # Localized control, Kalman filtering, and functional observability.
    @{ Category='05_localized_control_observability'; File='2015_wang_you_matni_localized_distributed_kalman_filters.pdf'; Url='https://nikolaimatni.github.io/papers/Necsys15_LDKF.pdf' },
    @{ Category='05_localized_control_observability'; File='2014_park_martins_lti_distributed_observers_iff.pdf'; Url='https://arxiv.org/pdf/1401.0926' },
    @{ Category='05_localized_control_observability'; File='2017_mitra_sundaram_distributed_observers_lti.pdf'; Url='https://arxiv.org/pdf/1608.01429' },
    @{ Category='05_localized_control_observability'; File='2017_mitra_sundaram_distributed_functional_observers.pdf'; Url='https://arxiv.org/pdf/1705.10891' },
    @{ Category='05_localized_control_observability'; File='2019_anderson_et_al_system_level_synthesis.pdf'; Url='https://arxiv.org/pdf/1904.01634' },
    @{ Category='05_localized_control_observability'; File='2022_kjellqvist_yu_infinite_horizon_sls.pdf'; Url='https://arxiv.org/pdf/2210.15815' },
    @{ Category='05_localized_control_observability'; File='2022_functional_observability_target_state_estimation.pdf'; Url='https://arxiv.org/pdf/2201.07256' },
    @{ Category='05_localized_control_observability'; File='2025_arbelaiz_et_al_how_far_to_share_measurements.pdf'; Url='https://arxiv.org/pdf/2406.14781' },
    @{ Category='05_localized_control_observability'; File='2026_zhao_et_al_discrete_lti_design_trilemma.pdf'; Url='https://arxiv.org/pdf/2603.20144' },
    @{ Category='05_localized_control_observability'; File='2026_fattore_et_al_discrete_lti_jordan_observer.pdf'; Url='https://arxiv.org/pdf/2603.10656' },
    @{ Category='05_localized_control_observability'; File='2026_jaiswal_berger_tomar_partial_state_estimation.pdf'; Url='https://arxiv.org/pdf/2604.00680' },

    # Communication, rate-distortion, and energy/accuracy tradeoffs.
    @{ Category='06_resource_tradeoffs'; File='2017_scaman_et_al_optimal_distributed_optimization.pdf'; Url='https://proceedings.mlr.press/v70/scaman17a/scaman17a.pdf' },
    @{ Category='06_resource_tradeoffs'; File='2010_tavildar_viswanath_wagner_gaussian_many_help_one.pdf'; Url='https://arxiv.org/pdf/0805.1857' },
    @{ Category='06_resource_tradeoffs'; File='2016_yang_grover_kar_lossy_in_network_computation.pdf'; Url='https://arxiv.org/pdf/1601.06224' },
    @{ Category='06_resource_tradeoffs'; File='2016_braverman_et_al_communication_lower_bounds_estimation.pdf'; Url='https://arxiv.org/pdf/1506.07216' },
    @{ Category='06_resource_tradeoffs'; File='2000_heinzelman_chandrakasan_balakrishnan_leach.pdf'; Url='https://pdos.csail.mit.edu/archive/decouto/papers/heinzelman00.pdf' },
    @{ Category='06_resource_tradeoffs'; File='2008_estimation_wireless_sensor_networks_tradeoffs.pdf'; Url='https://people.kth.se/~kallej/papers/wsn_estimation_ifac08.pdf' },
    @{ Category='06_resource_tradeoffs'; File='2026_perez_salesa_et_al_event_triggered_state_estimation.pdf'; Url='https://zaguan.unizar.es/record/172015/files/texto_completo.pdf' },

    # Sheaf/topological formulations of consistency and distributed data fusion.
    @{ Category='07_sheaves_topology'; File='2016_robinson_sheaves_sensor_integration.pdf'; Url='https://arxiv.org/pdf/1603.01446' },
    @{ Category='07_sheaves_topology'; File='2018_hansen_ghrist_spectral_theory_cellular_sheaves.pdf'; Url='https://arxiv.org/pdf/1808.01513' }
)

if ($Category.Count -gt 0) {
    $papers = @($papers | Where-Object { $_.Category -in $Category })
}

New-Item -ItemType Directory -Path $libraryRoot -Force | Out-Null

$downloaded = 0
$skipped = 0
$failed = [System.Collections.Generic.List[object]]::new()

foreach ($paper in $papers) {
    $categoryPath = Join-Path $libraryRoot $paper.Category
    New-Item -ItemType Directory -Path $categoryPath -Force | Out-Null
    $targetPath = Join-Path $categoryPath $paper.File
    $partPath = "$targetPath.part"

    if ((Test-Path -LiteralPath $targetPath) -and -not $ForceDownload) {
        $header = [System.IO.File]::ReadAllBytes($targetPath)[0..4]
        if ([System.Text.Encoding]::ASCII.GetString($header) -eq '%PDF-') {
            Write-Host "SKIP  $($paper.File)"
            $skipped++
            continue
        }
    }

    if (Test-Path -LiteralPath $partPath) {
        Remove-Item -LiteralPath $partPath -Force
    }

    Write-Host "GET   $($paper.File)"
    & curl.exe --location --fail --silent --show-error --retry 2 --retry-delay 2 --connect-timeout 15 --max-time $MaxSeconds --user-agent 'Mozilla/5.0 (research literature archive)' --output $partPath $paper.Url
    $curlExit = $LASTEXITCODE

    if ($curlExit -ne 0 -or -not (Test-Path -LiteralPath $partPath)) {
        $failed.Add([pscustomobject]@{ File=$paper.File; Url=$paper.Url; Reason="curl exit $curlExit" })
        continue
    }

    $bytes = [System.IO.File]::ReadAllBytes($partPath)
    $isPdf = $bytes.Length -ge 5 -and [System.Text.Encoding]::ASCII.GetString($bytes[0..4]) -eq '%PDF-'
    if (-not $isPdf) {
        $invalidPath = "$targetPath.invalid-download"
        Move-Item -LiteralPath $partPath -Destination $invalidPath -Force
        $failed.Add([pscustomobject]@{ File=$paper.File; Url=$paper.Url; Reason='response was not a PDF' })
        continue
    }

    Move-Item -LiteralPath $partPath -Destination $targetPath -Force
    $downloaded++
}

Write-Host ""
Write-Host "Downloaded: $downloaded; skipped: $skipped; failed: $($failed.Count)"
if ($failed.Count -gt 0) {
    $failed | Format-Table -AutoSize
    exit 1
}
