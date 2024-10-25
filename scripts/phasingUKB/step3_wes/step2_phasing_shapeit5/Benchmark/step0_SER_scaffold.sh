#!/bin/bash

threads=8

for CHR in {1..22}; do


	CUT=0.01
	ODIR=UKB_PHASING_EXOME_ARRAY/step4_validation/Scaffold/chr${CHR}

	dx mkdir -p ${ODIR}

	
	
	INPUT=/mnt/project/UKB_PHASING_EXOME_ARRAY/step2_phasing_common/chr${CHR}/UKB_chr${CHR}.exome_array.FULL.shapeit5.common.bcf
	families=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/samples.families.txt
	BCF_estimation=UKB_chr${CHR}.trios_kids.bcf
	

	#dx run app-swiss-army-knife -icmd="cat ${families} | grep -v 'NA' > trios.families.txt && cat trios.families.txt | cut -f 1 > trios.kids.txt && bcftools view --threads ${threads} --force-samples -S trios.kids.txt -Ob -o ${BCF_estimation} ${INPUT} && bcftools index ${BCF_estimation} --threads ${threads}" --tag chr${CHR} --tag bench_${chunk} --instance-type mem1_ssd1_v2_x8 --folder="./${ODIR}" --name bench1_${CHR} --priority normal -y


	
	SUB=/mnt/project/Phasing/PhasingWGS/step1_preparedata/validation.samples.txt
	OUT=validation_chr${CHR}.bcf
	BCF=/mnt/project/UKB_PHASING_EXOME_ARRAY/step0_merge/chr${CHR}/UKB.chr${CHR}.exome_array.full.sorted.bcf
	
	dx run app-swiss-army-knife --folder="./${ODIR}" -icmd="bcftools view --threads ${threads} -Ob -o $OUT -S $SUB --force-samples $BCF && bcftools index $OUT --threads ${threads}" --tag subsetV --tag benchWGS --instance-type mem1_ssd1_v2_x8 --name benchWGS_subsetV --priority normal -y
		
		



	BCF=/mnt/project/UKB_PHASING_EXOME_ARRAY/step4_validation/Scaffold/chr${CHR}/${BCF_estimation}
	VAL=/mnt/project/UKB_PHASING_EXOME_ARRAY/step4_validation/Scaffold/chr${CHR}/validation_chr${CHR}.bcf
	FRQ=/mnt/project/UKB_PHASING_EXOME_ARRAY/step1_chunks/chr${CHR}/UKB.chr${CHR}.exome_array.sorted.sites.bcf
	PED=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/samples.families.txt
	OUT=SER.chr${CHR}.fq
	LOG=SER.chr${CHR}.fq.log

	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="./${ODIR}/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $BCF --frequency $FRQ --pedigree $PED --region chr$CHR --output $OUT --log $LOG --thread $threads" --tag benchWGS --tag switchSHP4 --instance-type mem2_ssd2_v2_x8 --priority low --name benchWGS_chr${CHR} -y


		

done







