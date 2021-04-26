## Alignment

Aligment is performed by sentieon, either sharded or non-sharded. In the sharded setup the alignment is split up into a predefined amount of shards, see nextflow.config -> params -> bwa_shards = 8. Sharding alignment speeds up alignment for a specific sample and cen be enabled by adding the flag "--sharded" to the start-command. Amount of shards are preferable set to the maximum number of nodes for a specific cluster. For a 30x human WGS sample on 8 nodes (56 cores each) this would take around 20 minutes. Unsharded sentieon bwa takes around 2-3 hours.

The final bwa is then sent to distributed locuscollector, dedup and base quality recalibration scoring. 