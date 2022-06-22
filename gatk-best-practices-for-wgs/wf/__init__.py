"""
Call variants using GATK best-practices recommendations
"""

from cgitb import small
import subprocess
from pathlib import Path
from latch import small_task, workflow
from latch.types import LatchFile, LatchDir

@small_task
def samtools_faidx_task(reference: LatchFile) -> True:

    #output
    reference_index = Path("reference.fa.fai").resolve()

    _samtools_faidx_cmd = [
        "samtools",

        "faidx",
        reference.local_path,

    ]

    subprocess.run(_samtools_faidx_cmd, check=True)

    _mv_cmd = [         # testing this mv command in 0.1.23
        "mv",
        "*.fai",
        "latch:///reference.fa.fai",
    ]

    subprocess.run(_mv_cmd)
    return True


@small_task
def gatk_indexGenome_task(reference: LatchFile) -> LatchFile:

    #output
    gatk_genomeIndex = Path("reference.dict").resolve()

    _gatk_index_cmd = [
        "gatk",
        "CreateSequenceDictionary",

        "--REFERENCE",
        reference.local_path,

        "--OUTPUT",
        str(gatk_genomeIndex)

    ]

    subprocess.run(_gatk_index_cmd)

    return LatchFile(str(gatk_genomeIndex), "latch:///reference.dict")

@small_task
def gatk_indexVcf_task(vcf: LatchFile) -> LatchFile:

    #output
    gatk_vcfIndex = Path("known_sites.vcf.idx").resolve()

    _gatk_indexFeatureFile_cmd = [

        "gatk",
        "IndexFeatureFile",

        "-I",
        vcf.local_path,

        "-O",
        str(gatk_vcfIndex)
    ]

    subprocess.run(_gatk_indexFeatureFile_cmd)

    return LatchFile(str(gatk_vcfIndex), "latch:///known_sites.vcf.idx")


@small_task
def samtools_sort_task(bam_file: LatchFile) -> LatchFile:

    #output
    sorted_bam = Path("sorted.bam").resolve()

    _samtools_sort_cmd = [
        "samtools", 
        "sort",
        "-o",
        str(sorted_bam),
        "-O",
        "bam",
        bam_file.local_path,
    ]

    subprocess.run(_samtools_sort_cmd)

    return LatchFile(str(sorted_bam), "latch:///sorted.bam")

@small_task 
def mark_duplicates_task(sorted_bam: LatchFile) -> LatchFile:

    #using coordinate-sorted data to increase processing speed

    # A reference to our output.
    sorted_dedup_bam = Path("sorted_dedup.bam").resolve()

    duplication_metrics_file = Path("duplication_metrics.txt")

    _gatk_markDuplicates_cmd = [
        "gatk",
        "MarkDuplicates",

        "--INPUT",
        sorted_bam.local_path,

        "--OUTPUT",
        str(sorted_dedup_bam),

        "--METRICS_FILE",
        str(duplication_metrics_file),

        "--REMOVE_DUPLICATES",
        "true",

        "--REMOVE_SEQUENCING_DUPLICATES",
        "true",

        "--CREATE_INDEX",
        "true",
    ]

    subprocess.run(_gatk_markDuplicates_cmd)

    return LatchFile(str(sorted_dedup_bam), "latch:///sorted_dedup.bam")

@small_task
def samtools_index_task(bam_input: LatchFile) -> LatchFile:

    #output
    bam_index_output = Path("indexed_output.bam.bai").resolve()

    _samtools_index_cmd = [
        "samtools",
        "index",

        bam_input.local_path,
        str(bam_index_output),
    ]

    subprocess.run(_samtools_index_cmd)

    return LatchFile(str(bam_index_output), "latch:///indexed_output.bam.bai")

@small_task
def baseRecalibrator_task(sorted_dedup_bam: LatchFile, 
                            reference: LatchFile, 
                            known_sites: LatchFile) -> LatchFile:

    # output
    recal_table = Path("recal_data.table").resolve()

    _gatk_Recalibrator_cmd = [
        "gatk",
        "BaseRecalibrator",

        "--input",
        sorted_dedup_bam.local_path,

        "--reference",
        reference.local_path,

        "--output",
        str(recal_table),

        "--known-sites",
        known_sites.local_path,
    ]

    subprocess.run(_gatk_Recalibrator_cmd)
    return LatchFile(str(recal_table), "latch:///recal_data.table")

@small_task
def applybqsr_task(sorted_dedup_bam: LatchFile, 
                    reference: LatchFile, 
                    recal_file: LatchFile) -> LatchFile:

    #output
    sorted_dedup_bqsr_bam = Path("sorted_dedup_bqsr_recal.bam").resolve()

    _gatk_applyBQSR_cmd = [
        "gatk",
        "ApplyBQSR",

        "--bqsr-recal-file",
        recal_file.local_path,

        "--input",
        sorted_dedup_bam.local_path,

        "--ouput",
        str(sorted_dedup_bqsr_bam),

        "--reference",
        reference.local_path,
    ]

    subprocess.run(_gatk_applyBQSR_cmd)
    return LatchFile(sorted_dedup_bqsr_bam, "latch:///sorted_dedup_bqsr_recal.bam")

# include these steps in a later version

# #Post BQSR recall table
# gatk BaseRecalibrator \
#         --input "${sample}"_sorted_dedup_BQSR_recal.bam \
#         --output "${sample}"_post_recal_data.table \
#         --reference "${reference}"  \
#         --known-sites "${vcf}" \
#         --tmp-dir tmp 


# # Analyse covariates (compare before and After BQSR)
# # Evaluate and compare base quality score recalibration tables

# gatk AnalyzeCovariates \
#         -before "${sample}"_recal_data.table \
#         -after "${sample}"_post_recal_data.table \
#         -plots "${sample}"_AnalyzeCovariates.pdf


#  Variant calling HaplotypeCaller
# Call germline SNPs and indels via local re-assembly of haplotypes
# call potential variant sites per sample and save results in GVCF format

@small_task
def haplotypecaller_task(sorted_dedup_bqsr_bam: LatchFile, 
                            reference: LatchFile) -> LatchFile:

    #output
    gvcf = Path("sample.g.vcf.gz").resolve()

    _gatk_haplotype_caller_cmd = [
        "gatk",
        "HaplotypeCaller",

        "--reference",
        reference.local_path,

        "--input",
        sorted_dedup_bqsr_bam.local_path,

        "--output",
        str(gvcf),

        "--ERC",
        "GVCF",
    ]

    subprocess.run(_gatk_haplotype_caller_cmd)
    return LatchFile(str(gvcf), "latch:///sample.g.vcf.gz")


@workflow
def gatk_best_practices(
    bam_file: LatchFile,
    reference: LatchFile,
    known_sites: LatchFile
    ) -> LatchFile:

    """Description...

    GATK best-practices workflow
    
    ----

    This workflow uses GATK-best practices to call GVCFs from genomes 
    GATK-best practises recommendations are used for quality control of variants called by Haplotyper.

    __metadata__:
        display_name: GATK best-practices workflow
        author:
            name: Ibra Lujumba
            email: ibra.lujumba@gmail.com
            github:
        repository:
        license:
            id: MIT

    Args:

        bam_file:
          BAM file.

          __metadata__:
            display_name: BAM file

        reference:
          Reference FASTA genome of the organism.

          __metadata__:
            display_name: Reference Genome      
            
        known_sites:
         VCF with known variant sites

          __metadata__:
            display_name: Known variant sites

    """
    # index reference genome
    samtools_index = samtools_faidx_task(reference = reference)

    gatk_genomeIndex = gatk_indexGenome_task(reference = reference)

    gatk_vcfIndex = gatk_indexVcf_task(vcf = known_sites)

    sorted_bam = samtools_sort_task(bam_file = bam_file)

    sorted_bam_index = samtools_index_task(bam_input = sorted_bam)

    sorted_dedup_bam = mark_duplicates_task(sorted_bam = sorted_bam)

    sorted_dedup_bam_index = samtools_index_task(bam_input = sorted_dedup_bam)


    #base recalibration steps
    data_table = baseRecalibrator_task(sorted_dedup_bam = sorted_dedup_bam,
                                        reference = reference,
                                        known_sites = known_sites)

    sorted_dedup_bqsr_bam = applybqsr_task(sorted_dedup_bam = sorted_dedup_bam,
                                            reference = reference, 
                                            recal_file = data_table)

    sorted_dedup_bqsr_bam_index = samtools_index_task(bam_input = sorted_dedup_bqsr_bam)

    #include analyse_covariates steps
    # include clean up steps

    # get gvcf
    return haplotypecaller_task(sorted_dedup_bqsr_bam = sorted_dedup_bqsr_bam, reference = reference)
    return baseRecalibrator_task(sorted_dedup_input = sorted_dedup_bam, 
                                    reference = reference, 
                                    known_sites = known_sites)


# include steps for genotyping human genomes using machine learning algorithms (CNN)
# include steps for genotyping non-human genomes