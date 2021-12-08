devtools::load_all()
t1<-read_computed_matrix("data/cgi_sorted_mainchr_kate_chip.mtx.gz")
t1$metadata
t1$calculate_mean_coverage_each_sample()
t1$sample_mean_coverage
t1$calculate_total_coverage_each_region()
t1$region_total_coverage
