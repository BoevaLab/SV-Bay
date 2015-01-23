SV-Bay contains 3 main steps.

1)Main clustering. On this step SV-Bay calculates all the statistics about fragment length distribution and clusters abnormal fragments.

python -B main_clustering.py -c config.yaml

2)Main probability. On this step SV-Bay calculates probability for each cluster to associate it with noise or real SV.

python -B main_probability.py -c config.yaml

3)Assembling. On this step SV-Bay assembles clusters to complex and simple SVs. And outputs final results.

python assemly_links1.py -c config/config.yaml > results

Assemling step can also exclude all germline mutations, if the data is provided.

To exclude normal reads:	run SV-Bay with normal reads, run assembly step with flag -n and name of the folder with germ-line clusters.

python -B main_clustering.py -c config_germ.yaml

python assemly_links1.py -c config/config.yaml -n 'cluster_files_germ/' > results


For each step should be the same config file. Config example chech in config/config.yaml.