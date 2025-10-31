

include { CALC_UNALIGNED_STATS } from './modules/calc_unaligned_stats.nf'
include { UNALIGNED_SUMMARY } from './modules/unaligned_summary.nf'

all_align_channel = channel.fromPath( params.aligned_dir + "/*")
                           .map{ it -> tuple(it.baseName.replace("all_alignments_", ""), it)}
unalign_channel = channel.fromPath( params.unaligned_dir + "/*")
                         .map{ it -> tuple(it.baseName, it)}
                           

ch_samples = unalign_channel
             .join(all_align_channel)
             .take(2)

ch_samples.view()     

workflow {
    ch_results = CALC_UNALIGNED_STATS(ch_samples)
    
    UNALIGNED_SUMMARY( ch_results.collect())
}
