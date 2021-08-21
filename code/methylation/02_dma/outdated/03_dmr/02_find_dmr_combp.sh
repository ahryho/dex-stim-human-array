#!/bin/bash

src_data_pre="/Users/anastasiia_hry/bio/datasets/methylation/20_DMA/02_dmp"
dmr_pre="/Users/anastasiia_hry/bio/datasets/methylation/20_DMA/03_dmr"

mkdir -p ${dmr_pre}/comb-p


comb-p pipeline \
    -c 4 \
    --seed 0.00005 \
    --dist 5000 \
    --acf-dist 1000 \
    --step 200 \
    -p ${dmr_pre}/comb-p/dex \
    ${src_data_pre}/dmps_svs_for_combp.bed
    
    # --annotate hg19 \
    # --table refGene \
    # --region-filter-p 0.05 \
    

# awk 'NR == 1 || ($5 > 3 && $7 < 0.05)' ${PRE}.regions-p.bed > ${PRE}.regions.filtered.bed

# python -m cruzdb hg19 refGene ${PRE}.regions.filtered.bed cpgIslandExt > ${PRE}.anno.hg19.regions.filtered.bed
