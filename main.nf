#!/usr/bin/env nextflow

include { CircularPerBasePlot } from './modules/circular_cov_plot.nf'



workflow {
  
    
  
  gbk_ch = channel.fromPath("${params.path}" + "/GBK/*.gbk")
           .map{ it -> tuple( it.baseName, it)}
           
  cov_ch = channel.fromPath("${params.path}" + "/per_base_data/*.csv")
           .map{ it -> tuple( it.baseName.replace("per_base_data_",""), it)}
  
  gbk_cov_ch = gbk_ch.join( cov_ch)
  
  CircularPerBasePlot( gbk_cov_ch )
  
  
}

