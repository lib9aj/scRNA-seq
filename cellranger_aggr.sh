#bin/bash
module load cellranger/7.0.1
cellranger aggr --id=aggregated_output --csv=/scratch/lib9aj/aggregation.csv --normalize=mapped
