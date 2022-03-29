#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.obsids = 'no_obsid'
params.calid = 'no_calid'

params.data_dir = '/astro/mwavcs/nswainston/code/MWA-fast-image-transients/test_data'
params.cal_bin = ""

params.out_dir = 'results'

// imaging options
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
obsids = Channel.from(params.obsids.toString().split(",")).view()

process apply_cal {
    label 'cpu_one'
    label 'mwa_reduce'

    input:
    val obsid
    path solution_bin

    output:
    val obsid

    """
    applysolutions ${params.data_dir}/${obsid}.ms ${solution_bin}
    """
}


// process wsclean_split {
//     """
//     echo "##WSCLEAN 28s"
//     for i in `seq 0 3`;
//     do
//         j=$((8+56*i))
//         k=$((j+56))

//     echo "##WSCLEAN for 5s images"
//     if [[ ${dataintervals} -eq 240 ]]
//     then
//         # 240 -> flag first 4s and last 6s
//         #intervals='-interval 8 229'
//         start=8
//         end=218 #(note that interval 218 to 228 will be imaged. Only 0.5s interval 229 is not imaged)
//     else
//         # 224 -> process as below
//         #intervals=''
//         end=0
//         end=224
//     fi
//     for i in `seq ${start} 10 ${end}`;
//     do
//         j=$((i+10))

//     echo "##WSCLEAN for 0.5s images"
//     if [[ ${dataintervals} -eq 240 ]]
//     then
//         # 240 -> flag first 4s and last 6s
//         #intervals='-interval 8 229'
//         start=8
//         end=228 #(note that interval 228 to 229 is imageable, which is done by the code below)
//     else
//         # 224 -> process as below
//         #intervals=''
//         end=0
//         end=224
//     fi
//     for i in `seq ${start} ${end}`;
//     do
//         j=$((i+1))
//     """
// }


process wsclean {
    label 'cpu_all'
    label 'wsclean'
    publishDir "${params.out_dir}", mode: "copy", pattern: "*dirty.fits"
    // only publish image.fit if clean otherwise they are redundant
    publishDir "${params.out_dir}", mode: "copy", pattern: "*image.fits", enabled: params.clean

    input:
    val label
    val obsid

    output:
    tuple val(obsid), path("*fits")
    """
    wsclean -name ${obsid}_${label} \
        -size ${params.imsize} ${params.imsize} \
        -abs-mem ${task.memory.toGiga()} \
        -weight briggs -1 -mfs-weighting \
        -scale ${params.pixscale} \
        -pol xx,yy -minuv-l 30 \
        ${params.beamsize} -make-psf \
        ${clean_command} \
        -mgain ${params.mgain} \
        ${params.data_dir}/${obsid}.ms
    """
}


process make_primary_beam {
    label 'cpu_all'
    label 'mwa_reduce'

    input:
    tuple val(obsid), path(fits)

    output:
    tuple val(obsid), path("${obsid}_beam-MFS-*.fits")

    """
    # Dirty hack to link the h5 file until I work out the environmental variable I need
    mkdir -p mwapy/data
    ln -s /pawsey/mwa/mwa_full_embedded_element_pattern.h5 mwapy/data/mwa_full_embedded_element_pattern.h5

    echo "## Making primary beam"
    beam -2016 -proto *XX-image.fits -ms ${params.data_dir}/${obsid}.ms -name ${obsid}_beam-MFS
    """
}

process make_stokes {
    label 'cpu_all'
    label 'mwa_reduce'
    publishDir "${params.out_dir}", mode: "copy"

    input:
    val label
    tuple val(obsid), path(primary_beam_fits), path(image_fits)

    output:
    path "*pbcorr-*.fits"

    """
    echo "## Creating stokes iquv images"
    pbcorrect ${obsid}_${label} image.fits beam-MFS ${obsid}_${label}-pbcorr
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
    wsclean(
        Channel.value("2m"),
        apply_cal.out
    )
    make_primary_beam( wsclean.out )
    make_stokes(
        Channel.value("2m"),
        make_primary_beam.out.cross( wsclean.out ).map{ it -> [it[0][0], it[0][1], it[1][1]] },
    )
}