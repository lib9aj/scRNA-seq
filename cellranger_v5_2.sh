# bsub -L /bin/bash -W 40:00 -n 4 -R "span[ptile=4]" -M 128000 -e %J.err -o %J.out ./cellranger_v5.sh
module load cellranger/5.0.1
cellranger count --id=CRR794362 --transcriptome=/database/cellranger/mm10-v2 --fastqs=./fastqs
