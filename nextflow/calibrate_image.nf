#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.calid = 'no_calid'

params.doaoflagger = false

params.modeldir = "/astro/mwavcs/nswainston/code/MWA-fast-image-transients/models"
params.calmodel = 'model-PicA-88comp_withalpha.txt'
params.data_dir = '/astro/mwavcs/nswainston/code/MWA-fast-image-transients/test_data'

params.out_dir = 'results'


params.help = false
if ( params.help ) {
    help = """calibrate_image.nf: A pipeline that will download MWA data and calibrate it.
             |Required argurments:
             |  --calid     MWA Observation ID you want to calibrate [no default]
             |
             |Optional arguments:
             |  --download_dir
             |              The directory of already downloaded data (with the ASVO).
             |  -w          The Nextflow work directory. Delete the directory once the processs
             |              is finished [default: ${workDir}]""".stripMargin()
    println(help)
    exit(0)
}


process calibrate {
    label 'cpu_all'
    label 'mwa_reduce'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    val description
    // Used done channels to force aoflagger to be serial
    val done

    output:
    path "*_solutions_${description}.bin"
    val "done"

    """
    solutions=${params.calid}_${params.calmodel.split(".txt")[0]}_solutions_${description}.bin

    # calibrate
    calibrate -absmem ${task.memory.toGiga()} -j ${task.cpus} -m ${params.modeldir}/${params.calmodel} -minuv 20 -maxuv 2700 ${params.data_dir}/${params.calid}.ms \${solutions}

    # apply calibration
    applysolutions ${params.data_dir}/${params.calid}.ms \${solutions}
    """
}

process aoflagger {
    label 'cpu_all'
    label 'wsclean'

    input:
    each done
    val obsid

    output:
    val obsid
    """
    # run aoflagger
    # default is to work on the corrected data column
    aoflagger -j ${task.cpus} ${params.data_dir}/${obsid}.ms
    """
}

process aocal_plot {
    label 'cpu_one'
    label 'mwa_reduce'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path solution_bin

    output:
    path "*png"

    """
    aocal_plot.py --amp_max=2 ${solution_bin}
    """
}


workflow aoflagger_flow {
    take:
        done
    main:
        aoflagger(
            done,
            Channel.value(params.calid)
        )
        calibrate(
            Channel.value("aoflagged"),
            aoflagger.out
        )
    emit:
        calibrate.out[0]
}

workflow {
    calibrate(
        Channel.value("initial"),
        Channel.value("Done")
    )
    if ( params.doaoflagger ) {
        // Used done channels to force aoflagger to be serial
        aoflagger_flow( calibrate.out[1] )
        solution_bin = aoflagger_flow.out
        plot_solutions = aoflagger_flow.out.mix( calibrate.out[0] )
    }
    else {
        solution_bin = calibrate.out[0]
        plot_solutions = calibrate.out[0]
    }
    aocal_plot( plot_solutions )
}