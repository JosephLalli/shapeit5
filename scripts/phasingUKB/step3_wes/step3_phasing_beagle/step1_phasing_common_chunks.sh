#!/bin/bash




threads=32


for CHR in 2; do

	CHUNKS=../step1_chunking/robin_chunks_exome/chunks_chr${CHR}.txt

	while read LINE; do

		SCAFFOLD_REG=$(echo $LINE | awk '{ print $3; }')
		SCAFFOLD_REG_START=$(echo ${SCAFFOLD_REG} | cut -d":" -f 2 | cut -d"-" -f1)
		SCAFFOLD_REG_END=$(echo ${SCAFFOLD_REG} | cut -d":" -f 2 | cut -d"-" -f2)
		SCAFFOLD_REG_NAME=${CHR}_${SCAFFOLD_REG_START}_${SCAFFOLD_REG_END}

		INPUT_REG=$(echo $LINE | awk '{ print $4; }')
		INPUT_REG_START=$(echo ${INPUT_REG} | cut -d":" -f 2 | cut -d"-" -f1)
		INPUT_REG_END=$(echo ${INPUT_REG} | cut -d":" -f 2 | cut -d"-" -f2)
		INPUT_REG_NAME=${CHR}_${INPUT_REG_START}_${INPUT_REG_END}


		ODIR=UKB_PHASING_EXOME_ARRAY/step0_merge/chr${CHR}/
	
	
		VCF=/mnt/project/UKB_PHASING_EXOME_ARRAY/step0_merge/chr${CHR}/UKB.chr${CHR}.exome_array.WO_parents.sorted.vcf.gz
		BGL=/mnt/project/docker/beagle.19Apr22.7c0.jar
		MAP=/mnt/project/data/plink_maps/plink.prefix.chr${CHR}.GRCh38.map	
		CHUNKS=../step1_chunking/chunks.chr${CHR}.txt
		CUT=0.01
		THREADS=64
		#ODIR=UKB_PHASING_EXOME_ARRAY/step5_phasing_beagle_common/chr${CHR}

		dx mkdir -p ${ODIR}

	
	
		OUT=UKB_chr${CHR}.exome_array.WO_parents.beagle.common.${SCAFFOLD_REG}.vcf.gz
		LOG=UKB_chr${CHR}.exome_array.WO_parents.beagle.common.${SCAFFOLD_REG}.log
		TIM=UKB_chr${CHR}.exome_array.WO_parents.beagle.common.${SCAFFOLD_REG}.time
	
	
		dx run app-swiss-army-knife --folder="./$ODIR/" -icmd="/usr/bin/time -vo $TIM java -Xmx256G -jar $BGL gt=$VCF map=$MAP out=$OUT nthreads=${threads} chrom=${SCAFFOLD_REG} && bcftools index -f $OUT\.vcf.gz --threads 32" --tag benchWGS --tag beagle5.4 --tag chr20 --instance-type mem2_ssd1_v2_x32 --priority normal --name beagle_chr${CHR}_${SCAFFOLD_REG} -y
	
	
	done < $CHUNKS

done
