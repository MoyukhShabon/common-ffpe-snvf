#!/usr/env/bin bash

rsync -chaPv ../../../gdc/vcf-pass-filter-eval-all-snv/roc-prc-auc/precrec/combined_auc_table.tsv tcga.matched.f-ffpe_f-ff.all-snvs.tsv
rsync -chaPv ../../../gdc/vcf-pass-filter-eval/mutect2/roc-prc-auc/combined_auc_table.tsv tcga.mutect2.matched.f-ffpe_f-ff.ct-snvs.tsv
rsync -chaPv ../../../gdc/vcf-pass-filter-eval/muse/roc-prc-auc/combined_auc_table.tsv tcga.muse.matched.f-ffpe_f-ff.ct-snvs.tsv
rsync -chaPv ../../../gdc/vcf-pass-filter-eval/varscan2/roc-prc-auc/combined_auc_table.tsv tcga.varscan2.matched.f-ffpe_f-ff.ct-snvs.tsv
rsync -chaPv ../../../gdc/vcf-pass-filter-eval/somaticsniper/roc-prc-auc/combined_auc_table.tsv tcga.somaticsniper.matched.f-ffpe_f-ff.ct-snvs.tsv
rsync -chaPv ../../../gdc/vcf-ffpe-eval/plots/aucs.tsv tcga.matched.uf-ffpe_uf-ff.ct-snvs.tsv
rsync -chaPv ../../../gdc/mutect2_tumor-only/evaluations/aucs.tsv tcga.tumor-only.f-ffpe_f-ff.ct-snvs.tsv