# Align reads from a pair of fastq files using bwa mem
# and produce an aligned, duplicate-marked, unsorted BAM file.
task fastq_bwa_mem_paired {

	File ref_fasta
	# NB  64-bit bwa index files should be used
	File ref_fasta_amb
	File ref_fasta_ann
	File ref_fasta_bwt
	File ref_fasta_pac
	File ref_fasta_sa
	# optional list of contigs that are "alternative"
	File? ref_fasta_alt
	Array[File] fastqs_r1
	Array[File] fastqs_r2
	# rg header file must use `\t` and `\n` literally
	File rg_header
	String sample_id

	Int cpu
	Int memory_gb
	Int preemptible_tries

	Int diskspace_gb = ceil(4 * (size(fastqs_r1[0], "GB") + size(fastqs_r2[0], "GB")) * length(fastqs_r1))

	command <<<
		set -euo pipefail
		
		if [[ ! -f ${ref_fasta_sa} ]]; then
			#	index the reference fasta if it is missing
			bwa index -a bwtsw ${ref_fasta}
		fi
		
		# align reads, mark duplicates, and create unsorted bam file with fast compression
		# difference in duplicate marking between samblaster and Picard:
		#  the first encountered read-pair of a duplicate set will considered as the 
		#  prototype instead of the 'best' read-pair
		bwa mem -Y -t ${cpu} \
			-R $(cat ${rg_header}) \
			${ref_fasta} <(cat ${sep=' ' fastqs_r1}) <(cat ${sep=' ' fastqs_r2}) \
			| samblaster | samtools view -b1 - \
			> ${sample_id}_aligned.bam
	>>>

	output {
		File bam = "${sample_id}_aligned.bam"
	}

	runtime {
		# docker: "djhshih/seqkit:0.1"
		memory: "${memory_gb} GB"
		cpu: "${cpu}"
		disks: "local-disk ${diskspace_gb} HDD"
		preemptible: preemptible_tries
	}

}

# Sort BAM file by coordinate and index.
task bam_sort_coord {

	File input_bam
	String sample_id

	# Sorting is done is two phases
	# 1. Sort chunks in parallel
	# 2. Merge chunks
	# Phase 2 is done in single-threaded mode, which reduces the benefit
	# of multithreading (there will be little speed improvement above 8 threads)
	Int cpu
	Int memory_gb
	Int preemptible_tries

	Int memory_mb_sort = floor(0.8 * memory_gb * 1024)
	Int diskspace_gb = ceil(6 * size(input_bam, "GB"))

	command <<<
		set -euo pipefail
		
		# sort alignments with ample RAM to avoid disk IO and create final bam file
		sambamba sort -t ${cpu} -m ${memory_mb_sort}M --tmpdir ./tmp \
			-o ${sample_id}.bam ${input_bam}

		mv ${sample_id}.bam.bai ${sample_id}.bai
	>>>

	output {
		File bam = "${sample_id}.bam"
		File bai = "${sample_id}.bai"
	}

	runtime {
		# docker: "djhshih/seqkit:0.1"
		memory: "${memory_gb} GB"
		cpu: "${cpu}"
		disks: "local-disk ${diskspace_gb} HDD"
		preemptible: preemptible_tries
	}

}


workflow fastq_align_paired {
	String sample_id

	call fastq_bwa_mem_paired {
		input:
			sample_id = sample_id
	}

	call bam_sort_coord {
		input:
			sample_id = sample_id,
			input_bam = fastq_bwa_mem_paired.bam
	}

	output {
		File bam = bam_sort_coord.bam
		File bai = bam_sort_coord.bai
	}
}
