/data/bnf/sw/nextflow/nextflow run main.nf -resume \
			       -params-file params_small.json \
				   --csv shards.csv \
			       -with-singularity /data/bnf/dev/bjorn/nextflow_test/container_2019-04-18.sif \
			       -with-dag test.dag.png
 
