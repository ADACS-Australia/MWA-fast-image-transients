#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.obsids = 'no_obsid'
params.calid = 'no_calid'

params.data_dir = '/astro/mwavcs/nswainston/code/MWA-fast-image-transients/test_data'
params.cal_bin = ""

params.out_dir = 'results'

// imaging options
params.image_dur = 120.0
params.imsize = 4096
params.pixscale = "32asec"
params.mgain = 1
params.beamsize = "-elliptical-beam"
params.clean = false


params.help = false
if ( params.help ) {
    help = """mwa_fast_image.nf: A pipeline that will download MWA data, calibrate and image it.
             |Required argurments:
             |  --obsids    Observation IDs you want to process seperated by a comma.
             |              Eg. 1111111111,2222222222,3333333333 [no default]
             |  --calid     Calibration ID you want to apply [no default]
             |  --cal_bin   Calibration solution .bin file [no default]
             |
             |Imaging arguments:
             |  --image_dur The duration of each image we will make in seconds [default: ${params.image_dur}]
             |  --imsize    Image size will be imsize x imsize pixels [default: ${params.imsize}]
             |  --pixscale  Image pixel scale [default: ${params.pixscale}]
             |  --mgain     Mgain value in wsclean [default: ${params.mgain}]
             |  --beamsize  Circular beam size in arcsecond [default: is no circular beam]
             |  --clean     Clean image  [default: false]
             |
             |Optional arguments:
             |  --download_dir
             |              The directory of already downloaded data (with the ASVO).
             |  -w          The Nextflow work directory. Delete the directory once the processs
             |              is finished [default: ${workDir}]""".stripMargin()
    println(help)
    exit(0)
}


if ( params.beamsize != "-elliptical-beam" ) {
    // Not using default so convert the input to a command
    params.beamsize = "-beam-size ${params.beamsize}"
}

if ( params.clean ) {
    // Create clean command
    clean_command = "-join-polarizations -niter 10000000 -auto-mask 3 -auto-threshold 1"
}
else {
    clean_command = ""
}

if ( params.obsids == 'no_obsid' ) {
    println("Please input the observation IDs you would like to process with --obsids.")
    exit(0)
}
if ( params.calid == 'no_obsid' ) {
    println("Please input the calibration observation ID you would like to process with --calid.")
    exit(0)
}
if ( params.cal_bin == '' ) {
    println("Please input the calibration solution binary with --cal_bin.")
    exit(0)
}
obsids = Channel.from(params.obsids.toString().split(","))

process apply_cal {
    label 'cpu_one'
    label 'mwa_reduce'

    input:
    val obsid
    each solution_bin

    output:
    val obsid

    """
    applysolutions ${params.data_dir}/${obsid}.ms ${solution_bin}
    """
}


process wsclean {
    label 'cpu_all'
    label 'wsclean'
    publishDir "${params.out_dir}", mode: "copy", pattern: "*dirty.fits"
    // only publish image.fit if clean otherwise they are redundant
    publishDir "${params.out_dir}", mode: "copy", pattern: "*image.fits", enabled: params.clean
    time '3 h'

    input:
    val obsid

    output:
    tuple val(obsid), path("*fits")
    """
    # get metadata from measurment set
    nsamples=\$( taql "NELEMENTS([select distinct TIME from ${params.data_dir}/${obsid}.ms])" )
    time_res=\$( taql "select distinct INTERVAL from ${params.data_dir}/${obsid}.ms" | tail -n 1 )
    image_sample=\$(echo "${params.image_dur} / \$time_res" | bc)

    echo "nsamples: \$nsamples  time_res: \$time_res  image_sample: \$image_sample"

    for i in `seq 0 \${image_sample} \${nsamples}`; do
        j=\$((\$i+\${image_sample}))
        wsclean -name ${obsid}_${params.image_dur}s_s\${i}_e\${j}_t \
            -size ${params.imsize} ${params.imsize} \
            -abs-mem ${task.memory.toGiga()} \
            -weight briggs -1 -mfs-weighting \
            -scale ${params.pixscale} \
            -pol xx,yy -minuv-l 30 \
            ${params.beamsize} -make-psf \
            ${clean_command} \
            -interval \${i} \${j} \
            -mgain ${params.mgain} \
            ${params.data_dir}/${obsid}.ms
    done
    """
}


process make_primary_beam {
    label 'cpu_all'
    label 'mwa_reduce'
    publishDir "${params.out_dir}", mode: "copy"

    input:
    tuple val(obsid), path(fits)

    output:
    tuple val(obsid), path("${obsid}_beam-MFS-*.fits")

    """
    # Dirty hack to link the h5 file until I work out the environmental variable I need
    mkdir -p mwapy/data
    ln -s /pawsey/mwa/mwa_full_embedded_element_pattern.h5 mwapy/data/mwa_full_embedded_element_pattern.h5

    echo "## Making primary beam"
    beam -2016 -proto *.fits -ms ${params.data_dir}/${obsid}.ms -name ${obsid}_beam-MFS
    """
}

process make_stokes {
    label 'cpu_all'
    label 'mwa_reduce'
    publishDir "${params.out_dir}", mode: "copy"

    input:
    tuple val(label), file(primary_beam_fits), file(image_fits)

    output:
    path "*pbcorr-*.fits"

    """
    echo "## Creating stokes iquv images"
    pbcorrect ${label}_t image.fits beam-MFS ${label}-pbcorr
    """
}

include { aoflagger } from './calibrate_image'

workflow {
    aoflagger(
        Channel.value("Done"),
        obsids,
    )
    apply_cal(
        aoflagger.out,
        Channel.fromPath( params.cal_bin )
    )
    wsclean( apply_cal.out  )
    make_primary_beam( wsclean.out.map{ it -> [it[0], it[1][0] ]} )
    make_stokes(
        // Use cross to give a primay beam with every image
        make_primary_beam.out.cross(
        // Group by their images
        wsclean.out.map{ it -> it[1] }.flatten().map{ it -> [ it.baseName.split("_t-")[0], it ] }.groupTuple()
        // Go back to the label just being obsid then it's ready to be crossed
        .map{ it -> [ it[0].split("_")[0], [ it[0], it[1] ] ] }
        // Reformate it to the tuple required
        ).map{ it -> [it[1][1][0], it[0][1], it[1][1][1]] }
    )
}