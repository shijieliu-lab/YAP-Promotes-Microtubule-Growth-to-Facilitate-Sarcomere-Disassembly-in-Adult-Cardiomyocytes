/bin/bash
#BSUB -W 5:00
#BSUB -n 2
#BSUB -M 20000
#BSUB -e %J.err
#BSUB -o %J.out

#module load deeptools
source ~/.bash_profile
module load deeptools/3.5.5

#bigwig samples--------------------------------
YAP5SA="/data/Molcardiolab/Liu_lab_Fengz/project/011625NRCM_YAP5SA_CUTRUN/bigwigs/NRCM_YAP5SA_rabbit_YAP_rn7.sorted.multi.Q30.dedup.bw"
YAPS94A="/data/Molcardiolab/Liu_lab_Fengz/project/011625NRCM_YAP5SA_CUTRUN/bigwigs/NRCM_YAPS94A_rabbit_YAP_rn7.sorted.multi.Q30.dedup.bw"

BED1="/data/Molcardiolab/Liu_lab_Fengz/project/011625NRCM_YAP5SA_CUTRUN/Heatmap/Bed_cluster/250310.YAP.rabbit.mat.cluster1.bed"
BED2="/data/Molcardiolab/Liu_lab_Fengz/project/011625NRCM_YAP5SA_CUTRUN/Heatmap/Bed_cluster/250310.YAP.rabbit.mat.cluster2.bed"
BED3="/data/Molcardiolab/Liu_lab_Fengz/project/011625NRCM_YAP5SA_CUTRUN/Heatmap/Bed_cluster/250310.YAP.rabbit.mat.cluster3.bed"

computeMatrix reference-point \
	-S $YAPS94A $YAP5SA \
	-R $BED1 $BED2 $BED3 \
	-a 2000 -b 2000 \
	--outFileName 250902.YAP.rabbit.3clusters.mat.gz \
	--sortRegions 'descend' \
	--samplesLabel   'YAP5SA_S94A' 'YAP5SA'  \
	--outFileNameMatrix  250902.YAP.rabbit.3clusters.matrix \
	--missingDataAsZero \
	--sortUsingSamples 1 

        plotHeatmap --matrixFile 250902.YAP.rabbit.3clusters.mat.gz \
            --outFileName 250902.YAP.rabbit.mat.3clusters.pdf \
            --outFileSortedRegions 250902.YAP.rabbit.mat.3clusters.bed \
            --dpi 300 \
            --sortRegions "descend" \
            --colorMap 'RdYlBu_r' \
            --boxAroundHeatmaps no \
            --legendLocation "lower-center" \
            --sortUsingSamples 1 \
            --zMax 35 35 35 \
