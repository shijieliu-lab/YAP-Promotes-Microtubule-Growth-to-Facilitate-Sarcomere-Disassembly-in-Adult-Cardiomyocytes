config {
    executor="lsf"
    // queue="normal"
    // lsf_request_options="-R 'rusage[mem=4000]' -n 2 -P CR_Mapping -J chipcustomMap.sjcab_mapping"

     commands {    
        fastqc {
           lsf_request_options="-W 01:00 -R 'rusage[mem=5000]' -n 6"
        }
	
	Trimgalore {
	   //lsf_request_options="-W 02:00 -R 'rusage[mem=5000]' -n 2 -P CR_Mapping -J chipcustomMap.sjcab_Trimgalore"
	   //lsf_request_options="-W 01:00 -M 1000 -n 2 -P CR_Mapping -J chipcustomMap.sjcab_Trimgalore"
	   lsf_request_options="-W 04:00 -M 1000 -n 2"  	
	}	

	Bowtie2 {
	   lsf_request_options="-W 06:00 -M 15000 -n 6"
	}

	Sortbam {
           lsf_request_options="-W 01:00 -M 10000 -n 2"
        }			

	FilterbyMappingQual {
           lsf_request_options="-W 01:00 -M 10000 -n 3"
        }

	Dedup {
           lsf_request_options="-W 01:00 -M 50000 -n 1"
        }

	Macs2 {
           lsf_request_options="-W 01:00 -M 5000 -n 3"
        }

	BamCoverage {
           lsf_request_options="-W 01:00 -M 10000 -n 5"
        }

	TEenrich {
           lsf_request_options="-W 02:00 -M 10000 -n 2"
        }

	//SplitAlleles {
        //   lsf_request_options="-W 04:00 -R 'rusage[mem=10000]' -n 4 -P CR_Mapping -J chipcustomMap.sjcab_SplitAlleles"
        //}

	//getUnmapped {
        //   lsf_request_options="-W 03:00 -R 'rusage[mem=5000]' -n 2 -P CR_Mapping -J chipcustomMap.sjcab_getUnmapped"
        //}
     }
}

DataDir = "/path/to/01.RawData"

def branches = [

    NRCM_YAP5SA_mouse_YAP     :  ["${DataDir}/NRCM_YAP5SA_mouse_YAP/NRCM_YAP5SA_mouse_YAP_CKDL250000333-1A_22M5WNLT4_L2_1.fq.gz",
                                    "${DataDir}/NRCM_YAP5SA_mouse_YAP/NRCM_YAP5SA_mouse_YAP_CKDL250000333-1A_22M5WNLT4_L2_2.fq.gz"],

	NRCM_YAP5SA_rabbit_YAP     :  ["${DataDir}/NRCM_YAP5SA_rabbit_YAP/NRCM_YAP5SA_rabbit_YAP_CKDL250000333-1A_22M5WNLT4_L2_1.fq.gz",
                                    "${DataDir}/NRCM_YAP5SA_rabbit_YAP/NRCM_YAP5SA_rabbit_YAP_CKDL250000333-1A_22M5WNLT4_L2_2.fq.gz"],

	NRCM_YAPS94A_mouse_YAP     :  ["${DataDir}/NRCM_YAPS94A_mouse_YAP/NRCM_YAPS94A_mouse_YAP_CKDL250000333-1A_22M5WNLT4_L2_1.fq.gz",
                                    "${DataDir}/NRCM_YAPS94A_mouse_YAP/NRCM_YAPS94A_mouse_YAP_CKDL250000333-1A_22M5WNLT4_L2_2.fq.gz"],

	NRCM_YAPS94A_rabbit_YAP     :  ["${DataDir}/NRCM_YAPS94A_rabbit_YAP/NRCM_YAPS94A_rabbit_YAP_CKDL250000333-1A_22M5WNLT4_L2_1.fq.gz",
                                    "${DataDir}/NRCM_YAPS94A_rabbit_YAP/NRCM_YAPS94A_rabbit_YAP_CKDL250000333-1A_22M5WNLT4_L2_2.fq.gz"],

]

fastqc ={
	doc title: "QC using fastQC",
	desc: "Using fastQC to generate quality control reports",
	author: "zchen0423@gmail.com"

	var SOUT : "./FastQC"
	var threads: 6

	def basename1 = new File(input1).getName().prefix.replaceAll(".fq","");
	def basename2 = new File(input2).getName().prefix.replaceAll(".fq","");
 	produce("${SOUT}/${branch.name}/${basename1}_fastqc.html","${SOUT}/${branch.name}/${basename2}_fastqc.html"){
       	exec """
        	module load fastqc/0.11.7 &&
         	mkdir -p ${SOUT} && 
         	mkdir -p ${SOUT}/${branch.name} &&
         	fastqc --threads ${threads} --extract  --outdir "${SOUT}/${branch.name}" $inputs
	""","fastqc"
    }
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
      		trim_galore --paired  --fastqc --output_dir Trimgalore/${branch.name}_trimgalore $inputs
	""","Trimgalore"    
   }
}

bw2_genome = "/users/zha5dp/project/01162025NRCM_YAP5SA_CUTRUN/rat/rn7"

Bowtie2 = {
	var genome : bw2_genome
   	var organism : "rn7"
   	var additional_flags : ""

   	def minisize=0
   	def maxisize=1000
   	def fqDir = "./Trimgalore/${branch.name}_trimgalore"

   	def r1 = new File(fqDir).list().find{it=~/_val_1.fq/}
   	def r2 = new File(fqDir).list().find{it=~/_val_2.fq/}

   	if(organism != "rn7"){
       		branch.spikein = "./bowtie2_mapping_${organism}/${branch.name}_${organism}.bam"
       		//pattern = "./unmapped_rn7/${branch.name}*.fq"
       		fqDir = "./unmapped_rn7"

       		r1 = "${branch.name}_1.fq" 
       		r2 = "${branch.name}_2.fq" 

       		maxisize = 700
       		minisize = 0
   	}

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
         		$additional_flags | samtools view -bS - > ./bowtie2_mapping_${organism}/${branch.name}_${organism}.bam) >& "./bowtie2_mapping_${organism}/${branch.name}_${organism}_bowtie2_log.txt"
       		""","Bowtie2"
     }
}

// RemoveUnmapped = {
// 	def basename = new File(input).getName().prefix
// 	produce("./mapped_only/${basename}.mapped.bam"){
//		exec """
//         		module load samtools/1.9.0 &&
//       		mkdir -p mapped_only &&
//	     		samtools view -F 4 -b $input > "./mapped_only/${basename}.mapped.bam"
//		""","RemoveUnmapped"
//	}
//}

FilterbyMappingQual = {
	def basename = new File(input).getName().prefix
    	def qual = 30
	produce("./filter_qual_Q${qual}/${basename}.Q${qual}.bam"){
		exec """
			module load samtools/1.9.0 sambamba/0.6.8 &&
           		mkdir -p ./filter_qual_Q${qual} &&
			sambamba view -p -t 3 -f bam -F "mapping_quality >= ${qual} and not (unmapped or mate_is_unmapped)" $input > ./filter_qual_Q${qual}/${basename}.Q${qual}.bam
		""","FilterbyMappingQual"
	}
}

Sortbam = {    
	produce("${input.prefix}.sorted.multi.bam"){
	exec """
		module load samtools/1.9.0 &&
		samtools sort $input -o "${input.prefix}.sorted.multi.bam"
	""","Sortbam"
	}    
}

PICARD="/usr/local/picard/2.18.22/picard.jar"
Dedup = {
	def basename = new File(input).getName().prefix
 	branch.basename = basename
	produce("./deduplicated/${basename}.dedup.bam"){
		exec """
			module load java/1.7.0u40 picard/2.18.22 &&
			mkdir -p deduplicated &&
			java -jar $PICARD MarkDuplicates REMOVE_DUPLICATES=true 
			input=$input 
			output=./deduplicated/${basename}.dedup.bam
			METRICS_FILE=./deduplicated/${branch.name}_metrics.txt
			VALIDATION_STRINGENCY=SILENT
			CREATE_INDEX=true
		""","Dedup"
	//	exec """
	//		module load java/1.7.0u40 picard/2.18.22 &&
	//		mkdir -p deduplicated &&
	//		java -jar $PICARD MarkDuplicates -REMOVE_DUPLICATES true 
	//		-input $input -output ./deduplicated/${basename}.dedup.bam
	//		-METRICS_FILE ./deduplicated/${branch.name}_metrics.txt
	//		-VALIDATION_STRINGENCY SILENT -CREATE_INDEX true
	//	""","Dedup"
	}
}

Macs2 = {
	produce("./macs2/${branch.name}_peaks.narrowPeak"){
	exec """
		module load MACS/2.1.4 &&
		mkdir -p macs2 &&
      		macs2 callpeak -t $input -f BAMPE -g mm --outdir ./macs2 -n ${branch.name} -B -q 0.00001 --nolambda --nomodel
	""","Macs2"
	forward input 
  }
}

// SamToBam = {
//	doc title : "Convert a sam to Bam"
//	transform("sam") to ("bam") {
//	exec """
//		module load gcc/6.2.0 boost/1.62.0 samtools/1.9 &&
//		samtools view -Sb $input > $output
//   	""","SamToBam"
// }
//}

// IndexBam = {
//	produce("${input}.bai"){
//    	exec """
//		module load samtools/1.9.0 &&
//     		samtools index $input "${input}.bai"
//    	"""
//   	}	
//   	forward input
//}

BamCoverage = {
	doc title : "Generate bigwig file"
        author: "zchen0423@gmail.com"
   	def basename = new File(input).getName().prefix
   	def scale_unit=10000
   	def effGenomeSize=2678524473 //https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html (2nd table)
   	produce("./bigwigs/${basename}.bw"){
     		exec """
			module load deeptools/2.0.0 && 
			mkdir -p ./bigwigs &&              
       			bamCoverage --bam $input.bam  
                   	--binSize 25 
                   	--extendReads 150
                   	--normalizeUsing RPKM 
                   	--outFileFormat bigwig                   
                   	--scaleFactor 1
                   	--numberOfProcessors 5
                   	-o ./bigwigs/${basename}.bw
    		 ""","BamCoverage"
	}
  	forward input
}

chrom_size = "/data/ZYChenlab/Zhiyuan/genomes_annotations/mm10/genomes/chrNameLength.txt"
te_bed = "/data/ZYChenlab/Zhiyuan/genomes_annotations/mm10/annotations/mm10_rmsk_TE.bed"

TEenrich = {

        def peakDir = "./macs2"
        def peak = new File(peakDir).list().find{it=~/${branch.name}_multi_peaks.narrowPeak/}
        produce("./TEenrich/${branch.name}_TEenrich.txt"){
                exec """
                        module load bedtools/2.27.0 &&
                        mkdir -p ./TEenrich &&
                        /data/ZYChenlab/Zhiyuan/seqprg/TEenrich-FET.pl 
                        /data/ZYChenlab/Zhiyuan/genomes_annotations/mm10/genomes/chrNameLength.txt
                        /data/ZYChenlab/Zhiyuan/genomes_annotations/mm10/annotations/mm10_rmsk_TE2.bed 
                        ${peakDir}/$peak > ./TEenrich/${branch.name}_TEenrich.txt
                ""","TEenrich"
        }
}

// SplitAleles = {
//	def SNPFile = "/n/groups/zhanglab/nadhir/genomes/mm10_BDF1_PWK/snp/all_SNPs_PWKPhJ_GRCm38.txt.gz"
//   	def dedup =  "./deduplicated/${branch.basename}.dedup.bam"
//   	def basename = new File(input).getName().prefix
//   	def parentDir = new File(dedup).getParent()
//
//	if(! new File("./allelic_data/${basename}.genome1.bam").exists() && ! new File("./allelic_data/${basename}.genome2.bam").exists() ){
//   	produce("./allelic_data/${basename}.genome1.bam","./allelic_data/${basename}.genome2.bam"){
//     	exec """           
//        	module load gcc/6.2.0 boost/1.62.0 samtools/1.9 &&
//       	mkdir -p ./allelic_data &&
//       	(SNPsplit --snp_file $SNPFile  --paired  $input.bam) >& ${branch.name}.snpsplit.log &&
//       	mv ${parentDir}/${basename}*.genome*.bam ./allelic_data/ &&
//        	mv ${parentDir}/${basename}*sortedByName.bam ./allelic_data/ &&
//       	mv ${parentDir}/${basename}*allele_flagged.bam ./allelic_data/ &&
//        	mv ${parentDir}/${basename}*.unassigned.bam ./allelic_data/ &&
//       	mv ${parentDir}/${basename}*SNPsplit*.txt ./allelic_data/ 
//     	""","SplitAleles"
//    }
//   }
//  	forward "./allelic_data/${basename}.genome1.bam","./allelic_data/${basename}.genome2.bam"
//}

//fixEcoliBam = {
//	produce("${input.prefix}.fixed.bam"){
//        exec """
//             samtools view -h ${input} | sed "s/U00096.3_Escherichia_coli_str._K-12_substr._MG1655,_complete/chr1/g" | samtools view -b -o ${input.prefix}.fixed.bam  -
//            """
//    }
//}

//getUnmapped = {   
//	var organism : "mm10"
//   	def bam = "./bowtie2_mapping_${organism}/${branch.name}_${organism}.bam"   
//   	produce("./unmapped_${organism}/${branch.name}_1.fq", "./unmapped_${organism}/${branch.name}_2.fq"){
//       		from(glob("./Trimgalor/${branch.name}_trimgalor/*_val*.fq.gz")){
//       			exec """
//        			module load gcc boost python/3.6.0 &&
//          			source ~/python-3.6.0/bin/activate &&
//          			mkdir -p unmapped_${organism} &&
//         			python getUmappedReads.py -b $bam -f $input1 -s $input2 -p ./unmapped_${organism}/${branch.name}
//       			""","getUnmapped"
//     }
//   }
//}

run {
  branches * [    fastqc + 
	          Trimgalore +
                  Bowtie2.using(genome: bw2_genome, organism: "rn7") +                
                  Sortbam + 
                  FilterbyMappingQual +            
                  Dedup + 
		  Macs2 + 
		  BamCoverage 
		  //+ SplitAleles + "%.bam" * [ Sortbam + IndexBam + BamCoverage ]// +            
                 // getUnmapped.using(organism: "mm10") +
                 // Bowtie2.using(genome: bw2_yeast_genome, organism: "yeast") + 
		 // Bowtie2.using(genome: bw2_ecoli_genome, organism: "ecoli")
             ]
}
