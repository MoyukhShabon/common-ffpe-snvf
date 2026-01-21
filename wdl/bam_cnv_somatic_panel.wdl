version 1.0

# Workflow for creating a GATK CNV Panel of Normals given a list of normal samples. Supports both WGS and WES.
#
# Notes:
#
# - The intervals argument is required for both WGS and WES workflows and accepts formats compatible with the
#   GATK -L argument (see https://gatkforums.broadinstitute.org/gatk/discussion/11009/intervals-and-interval-lists).
#   These intervals will be padded on both sides by the amount specified by padding (default 250)
#   and split into bins of length specified by bin_length (default 1000; specify 0 to skip binning,
#   e.g., for WES).  For WGS, the intervals should simply cover the autosomal chromosomes (sex chromosomes may be
#   included, but care should be taken to 1) avoid creating panels of mixed sex, and 2) denoise case samples only
#   with panels containing only individuals of the same sex as the case samples).
#
# - Intervals can be blacklisted from coverage collection and all downstream steps by using the blacklist_intervals
#   argument, which accepts formats compatible with the GATK -XL argument
#   (see https://gatkforums.broadinstitute.org/gatk/discussion/11009/intervals-and-interval-lists).
#   This may be useful for excluding centromeric regions, etc. from analysis.  Alternatively, these regions may
#   be manually filtered from the final callset.
#
#  A reasonable blacklist for excluded intervals (-XL) can be found at:
#   hg19: gs://gatk-best-practices/somatic-b37/CNV_and_centromere_blacklist.hg19.list
#   hg38: gs://gatk-best-practices/somatic-hg38/CNV_and_centromere_blacklist.hg38liftover.list (untested)
#
# - Example invocation:
#
#       java -jar cromwell.jar run cnv_somatic_panel_workflow.wdl -i my_parameters.json
#
#############

task intervals_preprocess_intervals {
    input {
      File? intervals
      File? blacklist_intervals
      File ref_fasta
      File ref_fasta_fai
      File ref_fasta_dict
      Int? padding
      Int? bin_length
      File? gatk4_jar_override

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts
    }

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000
    Int command_mem_mb = machine_mem_mb - 500

    # Determine output filename
    String filename = select_first([intervals, "wgs"])
    String base_filename = basename(filename, ".interval_list")

    command <<<
        set -eu
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx~{command_mem_mb}m -XX:+UseSerialGC" PreprocessIntervals \
            ~{"-L " + intervals} \
            ~{"-XL " + blacklist_intervals} \
            --reference ~{ref_fasta} \
            --padding ~{default="250" padding} \
            --bin-length ~{default="1000" bin_length} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ~{base_filename}.preprocessed.interval_list
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File preprocessed_intervals = "~{base_filename}.preprocessed.interval_list"
    }
}

task intervals_annotate_intervals {
    input {
      File intervals
      File ref_fasta
      File ref_fasta_fai
      File ref_fasta_dict
      File? mappability_track_bed
      File? mappability_track_bed_idx
      File? segmental_duplication_track_bed
      File? segmental_duplication_track_bed_idx
      Int? feature_query_lookahead
      File? gatk4_jar_override

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts
    }

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000
    Int command_mem_mb = machine_mem_mb - 500
    
    # Determine output filename
    String base_filename = basename(intervals, ".interval_list")

    command <<<
        set -eu
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx~{command_mem_mb}m -XX:+UseSerialGC" AnnotateIntervals \
            -L ~{intervals} \
            --reference ~{ref_fasta} \
            ~{"--mappability-track " + mappability_track_bed} \
            ~{"--segmental-duplication-track " + segmental_duplication_track_bed} \
            --feature-query-lookahead ~{default=1000000 feature_query_lookahead} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ~{base_filename}.annotated.tsv
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(ref_fasta, "GB")) + 50]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File annotated_intervals = "~{base_filename}.annotated.tsv"
    }
}

task bam_collect_counts {
    input{
      File intervals
      File bam
      File bam_idx
      File ref_fasta
      File ref_fasta_fai
      File ref_fasta_dict
      Array[String]? disabled_read_filters
      Boolean? enable_indexing
      String? format
      File? gatk4_jar_override
      String? gcs_project_for_requester_pays

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts
    }

    parameter_meta {
      bam: {
        localization_optional: true
      }
      bam_idx: {
        localization_optional: true
      }
    }

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = machine_mem_mb - 1000

    Boolean enable_indexing_ = select_first([enable_indexing, false])

    Array[String] disabled_read_filters_arr = if defined(disabled_read_filters) then prefix("--disable-read-filter ", select_first([disabled_read_filters])) else []

    # Sample name is derived from the bam filename
    String base_filename = basename(bam, ".bam")
    String format_ = select_first([format, "HDF5"])
    String hdf5_or_tsv_or_null_format =
        if format_ == "HDF5" then "HDF5" else
        (if format_ == "TSV" then "TSV" else
        (if format_ == "TSV_GZ" then "TSV" else "null")) # until we can write TSV_GZ in CollectReadCounts, we write TSV and use bgzip
    String counts_filename_extension =
        if format_ == "HDF5" then "counts.hdf5" else
        (if format_ == "TSV" then "counts.tsv" else
        (if format_ == "TSV_GZ" then "counts.tsv.gz" else "null"))
    String counts_index_filename_extension =
        if format_ == "HDF5" then "null" else
        (if format_ == "TSV" then "counts.tsv.idx" else
        (if format_ == "TSV_GZ" then "counts.tsv.gz.tbi" else "null"))
    Boolean do_block_compression =
        if format_ == "HDF5" then false else
        (if format_ == "TSV" then false else
        (if format_ == "TSV_GZ" then true else false))
    String counts_filename = "~{base_filename}.~{counts_filename_extension}"
    String counts_filename_for_collect_read_counts = basename(counts_filename, ".gz")

    command <<<
        set -eu
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        case ~{format_} in
            HDF5 | TSV | TSV_GZ)
                ;;
            *)
                echo "ERROR: Unknown format specified. Format must be one of HDF5, TSV, or TSV_GZ."
                exit 1
                ;;
        esac

        if [ ~{format_} = "HDF5" ] && [ ~{enable_indexing_} = "true" ]; then
            echo "ERROR: Incompatible WDL parameters. Cannot have format = HDF5 and enable_indexing = true."
            exit 1
        fi

        if [ ~{hdf5_or_tsv_or_null_format} = "null" ]; then
            echo "ERROR: Should never reach here."
            exit 1
        fi

        gatk --java-options "-Xmx~{command_mem_mb}m -XX:+UseSerialGC" CollectReadCounts \
            -L ~{intervals} \
            --input ~{bam} \
            --reference ~{ref_fasta} \
            --format ~{default="HDF5" hdf5_or_tsv_or_null_format} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ~{counts_filename_for_collect_read_counts} \
            ~{"--gcs-project-for-requester-pays " + gcs_project_for_requester_pays} \
            ~{sep=' ' disabled_read_filters_arr}

        if [ ~{do_block_compression} = "true" ]; then
            bgzip ~{counts_filename_for_collect_read_counts}
        fi

        if [ ~{enable_indexing_} = "true" ]; then
            gatk --java-options "-Xmx~{command_mem_mb}m -XX:+UseSerialGC" IndexFeatureFile \
                -I ~{counts_filename}
        fi
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(bam, "GB")) + 50]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        String entity_id = base_filename
        File counts = counts_filename
    }
}

task counts_create_read_count_panel_of_normals {
    input {
      String pon_entity_id
      Array[File] read_count_files
      Float? minimum_interval_median_percentile
      Float? maximum_zeros_in_sample_percentage
      Float? maximum_zeros_in_interval_percentage
      Float? extreme_sample_median_percentile
      Boolean? do_impute_zeros
      Float? extreme_outlier_truncation_percentile
      Int? number_of_eigensamples
      Int? maximum_chunk_size
      File? annotated_intervals   #do not perform explicit GC correction by default
      File? gatk4_jar_override

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts
    }

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = machine_mem_mb - 500

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx~{command_mem_mb}m -XX:+UseSerialGC" CreateReadCountPanelOfNormals \
            --input ~{sep=" --input " read_count_files} \
            --minimum-interval-median-percentile ~{default="10.0" minimum_interval_median_percentile} \
            --maximum-zeros-in-sample-percentage ~{default="5.0" maximum_zeros_in_sample_percentage} \
            --maximum-zeros-in-interval-percentage ~{default="5.0" maximum_zeros_in_interval_percentage} \
            --extreme-sample-median-percentile ~{default="2.5" extreme_sample_median_percentile} \
            --do-impute-zeros ~{default="true" do_impute_zeros} \
            --extreme-outlier-truncation-percentile ~{default="0.1" extreme_outlier_truncation_percentile} \
            --number-of-eigensamples ~{default="20" number_of_eigensamples} \
            --maximum-chunk-size ~{default="16777216" maximum_chunk_size} \
            ~{"--annotated-intervals " + annotated_intervals} \
            --output ~{pon_entity_id}.pon.hdf5
    >>>

    runtime {
        docker: "~{gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File read_count_pon = "~{pon_entity_id}.pon.hdf5"
    }
}




workflow bam_cnv_somatic_panel {

    input {
      ##################################
      #### required basic arguments ####
      ##################################
      File intervals
      File? blacklist_intervals
      Array[String] normal_bams
      Array[String] normal_bais
      String pon_entity_id
      File ref_fasta_dict
      File ref_fasta_fai
      File ref_fasta
      String gatk_docker

      ##################################
      #### optional basic arguments ####
      ##################################
      # If true, intervals_annotate_intervals will be run to create GC annotations and explicit GC correction
      # will be performed by the PoN generated by counts_create_read_count_panel_of_normals before PCA is performed on subsequent cases
      Boolean? do_explicit_gc_correction
      File? gatk4_jar_override
      Int? preemptible_attempts

      # Required if BAM/CRAM is in a requester pays bucket
      String? gcs_project_for_requester_pays

      ####################################################
      #### optional arguments for intervals_preprocess_intervals ####
      ####################################################
      Int? padding
      Int? bin_length
      Int? mem_gb_for_preprocess_intervals

      ##################################################
      #### optional arguments for intervals_annotate_intervals ####
      ##################################################
      File? mappability_track_bed
      File? mappability_track_bed_idx
      File? segmental_duplication_track_bed
      File? segmental_duplication_track_bed_idx
      Int? feature_query_lookahead
      Int? mem_gb_for_annotate_intervals

      ##############################################
      #### optional arguments forb count_collect_counts ####
      ##############################################
      String? collect_counts_format
      Int? mem_gb_for_collect_counts

      ##############################################################
      #### optional arguments for counts_create_read_count_panel_of_normals ####
      ##############################################################
      Float? minimum_interval_median_percentile
      Float? maximum_zeros_in_sample_percentage
      Float? maximum_zeros_in_interval_percentage
      Float? extreme_sample_median_percentile
      Boolean? do_impute_zeros
      Float? extreme_outlier_truncation_percentile
      Int? number_of_eigensamples
      Int? maximum_chunk_size
      Int? mem_gb_for_create_read_count_pon
    }

    Array[Pair[String, String]] normal_bams_and_bais = zip(normal_bams, normal_bais)

    call intervals_preprocess_intervals {
        input:
            intervals = intervals,
            blacklist_intervals = blacklist_intervals,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            padding = padding,
            bin_length = bin_length,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_preprocess_intervals,
            preemptible_attempts = preemptible_attempts
    }

    if (select_first([do_explicit_gc_correction, true])) {
        call intervals_annotate_intervals {
            input:
                intervals = intervals_preprocess_intervals.preprocessed_intervals,
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                ref_fasta_dict = ref_fasta_dict,
                mappability_track_bed = mappability_track_bed,
                mappability_track_bed_idx = mappability_track_bed_idx,
                segmental_duplication_track_bed = segmental_duplication_track_bed,
                segmental_duplication_track_bed_idx = segmental_duplication_track_bed_idx,
                feature_query_lookahead = feature_query_lookahead,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                mem_gb = mem_gb_for_annotate_intervals,
                preemptible_attempts = preemptible_attempts
        }
    }

    scatter (normal_bam_and_bai in normal_bams_and_bais) {
        call bam_collect_counts {
            input:
                intervals = intervals_preprocess_intervals.preprocessed_intervals,
                bam = normal_bam_and_bai.left,
                bam_idx = normal_bam_and_bai.right,
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                ref_fasta_dict = ref_fasta_dict,
                format = collect_counts_format,
                enable_indexing = false,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                mem_gb = mem_gb_for_collect_counts,
                preemptible_attempts = preemptible_attempts,
                gcs_project_for_requester_pays = gcs_project_for_requester_pays
        }
    }

    call counts_create_read_count_panel_of_normals {
        input:
            pon_entity_id = pon_entity_id,
            read_count_files = bam_collect_counts.counts,
            minimum_interval_median_percentile = minimum_interval_median_percentile,
            maximum_zeros_in_sample_percentage = maximum_zeros_in_sample_percentage,
            maximum_zeros_in_interval_percentage = maximum_zeros_in_interval_percentage,
            extreme_sample_median_percentile = extreme_sample_median_percentile,
            do_impute_zeros = do_impute_zeros,
            extreme_outlier_truncation_percentile = extreme_outlier_truncation_percentile,
            number_of_eigensamples = number_of_eigensamples,
            maximum_chunk_size = maximum_chunk_size,
            annotated_intervals = intervals_annotate_intervals.annotated_intervals,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_create_read_count_pon,
            preemptible_attempts = preemptible_attempts
    }

    output {
        File preprocessed_intervals = intervals_preprocess_intervals.preprocessed_intervals
        Array[String] read_counts_entity_ids = bam_collect_counts.entity_id
        Array[File] read_counts = bam_collect_counts.counts
        File read_count_pon = counts_create_read_count_panel_of_normals.read_count_pon
    }
}

