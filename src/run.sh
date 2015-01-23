#run instructions
#python vesion 2.7
#!/bin/sh
python -B main_clustering.py -c config.yaml
python -B main_probabilities.py
python assemly_links.py 
