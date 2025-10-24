#!/usr/env/bin bash

rsync -chaPv ../../../eval/mutect2-matched-normal_pass-orientation-filtered.vs.filtered-ff/FFG/roc-prc-auc/precrec/combined_auc_table.tsv seqc2.ffg.matched.f-ffpe_f-ff.tsv
rsync -chaPv ../../../eval/mutect2-matched-normal_pass-orientation-filtered.vs.filtered-ff/FFX/roc-prc-auc/precrec/combined_auc_table.tsv seqc2.ffx.matched.f-ffpe_f-ff.tsv
rsync -chaPv ../../../eval/mutect2-matched-normal_pass-orientation-filtered.vs.unfiltered-ff/FFG/roc-prc-auc/precrec/combined_auc_table.tsv seqc2.ffg.matched.f-ffpe_uf-ff.tsv
rsync -chaPv ../../../eval/mutect2-matched-normal_pass-orientation-filtered.vs.unfiltered-ff/FFX/roc-prc-auc/precrec/combined_auc_table.tsv seqc2.ffx.matched.f-ffpe_uf-ff.tsv
