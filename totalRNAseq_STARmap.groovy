config {
	executor="lsf"
	commands {    
        	Fastqc {
           		lsf_request_options="-W 01:00 -M 10000 -n 2"
        	}
	
		Trimgalore {
	   		lsf_request_options="-W 04:00 -M 10000 -n 2"  	

		}	
		
		hisat2_map {
			lsf_request_options="-W 08:00 -M 20000 -n 5"
		}
	
		STAR_map {
			lsf_request_options="-W 04:00 -M 30000 -n 5"
		}
		
		sortByName {
			lsf_request_options="-W 03:00 -M 10000 -n 3"
		}
				
		sortAndIndex {
			lsf_request_options="-W 03:00 -M 10000 -n 3"
		}

		uniqueAlign {
			lsf_request_options="-W 01:00 -M 10000 -n 3"
		}

		stringTie {
			lsf_request_options="-W 02:00 -M 10000 -n 3"
		}
					
		bamCoverage {
           		lsf_request_options="-W 01:00 -M 20000 -n 5"
        	}

		TEcount {
			lsf_request_options="-W 06:00 -M 40000 -n 3"
		}

		bowtie2_map {
			lsf_request_options="-W 01:00 -M 10000 -n 6"
		}
		
		ERCC_count {
			lsf_request_options="-W 01:00 -M 10000 -n 3"
		}
     }	
}

title: "RNA-Seq analysis pipeline (total RNA stranded)"
DataDir = "/path/to/01.RawData"

def branches = [

        YAP5SA_3        : ["${DataDir}/YAP5SA_3/YAP5SA_3_1.fq.gz",
                           "${DataDir}/YAP5SA_3/YAP5SA_3_2.fq.gz"],

        YAP5SA_4        : ["${DataDir}/YAP5SA_4/YAP5SA_4_1.fq.gz",
                            "${DataDir}/YAP5SA_4/YAP5SA_4_2.fq.gz"],

        YAPS94A_4       : ["${DataDir}/YAPS94A_4/YAPS94A_4_1.fq.gz",
                            "${DataDir}/YAPS94A_4/YAPS94A_4_2.fq.gz"],

        YAPS94A_5       : ["${DataDir}/YAPS94A_5/YAPS94A_5_1.fq.gz",
                            "${DataDir}/YAPS94A_5/YAPS94A_5_2.fq.gz"],
]


Fastqc ={
      doc title: "QC using fastQC",
      desc: "Using fastQC 0.11.4 to generate quality control reports",      
      author: "zchen0423@gmail.com"

    var SOUT : "./FastQC"
    var threds: 2
    def basename = new File(input).getName().prefix;

	exec """
      module load fastqc/0.11.7 &&
      mkdir -p ${SOUT}/${branch.name} &&
	  fastqc --threads ${threds} --extract  --outdir ${SOUT}/${branch.name} $inputs
      ""","Fastqc"
	forward inputs
}

Trimgalore = {
      def basename1 = new File(input1).getName().prefix.replaceAll(".fq","");
      def basename2 = new File(input2).getName().prefix.replaceAll(".fq","");

    produce("./Trimgalore/${branch.name}_trimgalore/${basename1}_val_1.fq.gz","./Trimgalore/${branch.name}_trimgalore/${basename2}_val_2.fq.gz"){
     exec """
       module load trimgalore/0.6.6 &&
       mkdir -p Trimgalore &&
        mkdir -p Trimgalore/${branch.name}_trimgalore &&
       trim_galore  --paired --clip_R2 5  --output_dir "Trimgalore/${branch.name}_trimgalore" $inputs
     ""","Trimgalore"
    }
}

hisat2_map = {
       def hist2_index="/users/zha5dp/RNA_Seq/Rat_annotation/hist2_index_rn7/hist2_index"
       def nbcpu = 5
       def outdir = "${branch.name}_hisat2"
       def alignment_splicesites = "/users/zha5dp/RNA_Seq/Rat_annotation/GCF_015227675.2_mRatBN7.2.splice.site.txt"
       branch.hisat2_dir = outdir

       produce("${outdir}/${branch.name}.bam"){
              exec """
              module load hisat2/2.2.1 samtools/1.9.0 &&              
              mkdir -p ${outdir} &&
              hisat2 -x ${hist2_index} 
                     -1 ${input1} 
                     -2 ${input2} 
                     --known-splicesite-infile $alignment_splicesites 
                     --sp 1000,1000
                     --no-mixed                     
                     --no-discordant
                     -p ${nbcpu} 
                     -k 100
                     --un-conc-gz ${outdir}/${branch.name}.unmapped.hisat2.gz
                     --met-stderr 
                     --new-summary 
                     --summary-file ${outdir}/${branch.name}.hisat2_summary.txt 
                     | samtools view -bS -F 4 -F 8 - > ${outdir}/${branch.name}.bam 
              ""","hisat2_map"
       }
}


STAR_map = {	
	def STAR_index="/users/zha5dp/RNA_Seq/Rat_annotation/rn7_STAR_index"
	def nbcpu = 5
	def outdir = "${branch.name}_1st_STAR"
	
	produce ("${outdir}/${branch.name}Aligned.out.bam"){
		exec """
		module load STAR/2.7.9 samtools/1.9.0 &&
		mkdir -p ${outdir} && 
		STAR --runThreadN ${nbcpu}
			--genomeDir ${STAR_index}
			--readFilesIn $inputs
			--outFileNamePrefix ${outdir}/${branch.name}
			--outSAMtype BAM Unsorted 
			--outFilterMultimapNmax 100
			--winAnchorMultimapNmax 100 
			--outSAMattributes NH HI NM MD XS AS 
			--genomeLoad NoSharedMemory
			--alignIntronMax 1000000
			--alignMatesGapMax 1000000
			--limitBAMsortRAM 30000000000  
			--outFilterType BySJout
			--readFilesCommand zcat	
		""","STAR_map"
	}
}

TEcount = {
	//def TEcount="/data/ZYChenlab/Zhiyuan/seqprg/TEtranscripts-master/bin/TEcount"
	def GTF="/data/ZYChenlab/Zhiyuan/genomes_annotations/mm10/annotations/gencode.vM25.annotation.gtf"
	def REPEAT="/data/ZYChenlab/Zhiyuan/genomes_annotations/mm10/annotations/mm10_rmsk_TE.gtf"
	def outdir = "${branch.name}_TEcount"
	produce ("${outdir}/${branch.name}.cntTable"){
		exec"""
		source ~/.bash_profile &&
		conda activate tetranscripts &&
		mkdir -p ${outdir} &&
		TEcount -b $input 
			 --GTF $GTF
			 --TE $REPEAT
			 --project ${branch.name}
			 --outdir $outdir
			 --stranded reverse
		""","TEcount"
	}	
	forward input
}

bw2_ERCC_genome = "/data/ZYChenlab/Zhiyuan/genomes_annotations/ERCC/ERCC92"
bw2_dCas9_KRAB_genome = "/data/ZYChenlab/Zhiyuan/genomes_annotations/dCas9_KRAB/dCas9_KRAB"

bowtie2_map = {
	var genome : bw2_ERCC_genome
	var organism : "ERCC"
	var additional_flags : ""

   	def minisize=0
   	def maxisize=1000
   	def fqDir = "./Trimgalore/${branch.name}_trimgalore"

	def r1 = new File(fqDir).list().find{it=~/_val_1.fq.gz/}
   	def r2 = new File(fqDir).list().find{it=~/_val_2.fq.gz/}

	produce("./bowtie2_mapping_${organism}/${branch.name}_${organism}.bam"){
       		exec """
         		module load gcc/9.3.0 boost/1.79.0 samtools/1.9.0 bowtie2/2.4.2 &&
         		mkdir -p bowtie2_mapping_${organism} &&
         		(bowtie2 -x ${genome}
         		-1 ${fqDir}/$r1
         		-2 ${fqDir}/$r2
         		-p 6 --no-unal --no-mixed --no-discordant
         		-I $minisize
         		-X $maxisize
         		$additional_flags | samtools view -bS - > ./bowtie2_mapping_${organism}/${branch.name}_${organism}.bam) >& "./bowtie2_mapping_${organism}/${branch.name}_${organism}_bowtie2_log.txt" &
       		""","bowtie2_map"
     }
}

sortAndIndex = {
       def nbcpu = 2
       produce("${input.prefix}.sorted.bam"){
              exec """
              module load  samtools/1.9.0 &&
              samtools sort $input -@ ${nbcpu} -o ${input.prefix}.sorted.bam &&
              samtools index ${input.prefix}.sorted.bam
              ""","sortAndIndex"
       }
}

sortByName  = {
	def nbcpu = 2
	produce("${input.prefix}.namesorted.bam"){
		exec """
              		module load  samtools/1.9.0 &&
              		samtools sort $input -@ ${nbcpu} -n -o ${input.prefix}.namesorted.bam   		
		""","sortByName"
	}
	forward input
}

uniqueAlign = {
	def basename = new File(input).getName().prefix
	produce("./uniqueAlign/${basename}_uniq.bam"){
		exec"""
			module load samtools/1.9.0 sambamba/0.6.8 &&
			mkdir -p ./uniqueAlign &&
			sambamba view -p -t 3 -f bam -F "mapping_quality >=60" $input > ./uniqueAlign/${basename}_uniq.bam
		""","uniqueAlign"
	}
}

stringTie = {
       def outdir = "stringTie/${branch.name}"
       def gtf = "/users/zha5dp/RNA_Seq/Rat_annotation/GCF_015227675.2_mRatBN7.2_genomic.gff"
       def nbcpu = 5
       produce("$outdir/${branch.name}.gene_abund.txt"){
       exec """
                module load stringtie/2.2.1 &&
		mkdir  -p $outdir &&
                stringtie $input 
                -o $outdir/${branch.name}_transcripts.gtf 
                -v 
                -G $gtf 
                -A $outdir/${branch.name}.gene_abund.txt
                -C $outdir/${branch.name}.cov_refs.gtf 
                -b $outdir/${branch.name}_ballgown 
                -e
                -p $nbcpu
		--rf
       ""","stringTie"   
    }
    forward input
}


bamCoverage = {
	def basename = new File(input).getName().prefix  
   	produce("./bigwigs/${basename}.bw"){
     		exec """
      			mkdir -p ./bigwigs &&
       			module load deeptools/3.5.5 samtools/1.9.0 && 
       			samtools index $input.bam &&
			bamCoverage --bam $input.bam  --binSize 25 
                   	--normalizeUsing RPKM --outFileFormat bigwig
                   	--scaleFactor 1 --numberOfProcessors 3
                   	-o ./bigwigs/${basename}.bw
     		""","bamCoverage"
 	}
 	forward input
}

ERCC_count = {
	var organism : "ERCC"
	produce("./bowtie2_mapping_${organism}/${branch.name}.cntTable"){
		exec """
			module load samtools/1.9.0 &&
			mkdir -p bowtie2_mapping_${organism} &&
			samtools idxstats $input | cut -f 1,3 > ./bowtie2_mapping_${organism}/${branch.name}.cntTable 
		""","ERCC_count"
	}
} 
run {    
   branches * [
   	       Fastqc +
               Trimgalore +
	       hisat2_map +
//	       STAR_map +
//	       TEcount +  
	       sortAndIndex +
	       uniqueAlign + 
               stringTie + 
               bamCoverage 
//	       bowtie2_map.using(genome: bw2_dCas9_KRAB_genome, organism: "dCas9_KRAB") + 
//	       bowtie2_map.using(genome: bw2_ERCC_genome, organism: "ERCC") +  
//	       ERCC_count
	
             ]
}
