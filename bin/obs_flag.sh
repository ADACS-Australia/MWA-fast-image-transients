#! /bin/bash

usage()
{
echo "obs_flag.sh [-d dep] [-q queue] [-t] obsnum
  -d dep      : job number for dependency (afterok)
  -q queue    : job queue, default=gpuq
  -t          : test. Don't submit job, just make the batch file
                and then return the submission command
  obsnum      : the obsid to process" 1>&2;
exit 1;
}

#initialize as empty
dep=
queue='-p gpuq'
tst=


# parse args and set options
while getopts ':td:q:' OPTION
do
    case "$OPTION" in
        d)
            dep=${OPTARG}
            ;;
	q)
	    queue="-p ${OPTARG}"
	    ;;
        t)
            tst=1
            ;;
        ? | : | h)
            usage
            ;;
  esac
done

# set the obsid to be the first non option
shift  "$(($OPTIND -1))"
obsnum=$1

# if obsid is empty then just print help

if [[ -z ${obsnum} ]]
then
    usage
fi

depend=
if [[ ! -z ${dep} ]]
then
    depend="--dependency=afterok:${dep}"
fi


base=/astro/mwasci/phancock/RadioSchool/MWA-fast-image-transients/

script="${base}queue/flag_${obsnum}.sh"
cat ${base}/bin/flag.tmpl | sed -e "s:OBSNUM:${obsnum}:g" \
                                -e "s:BASEDIR:${base}:g" > ${script}

output="${base}queue/logs/flag_${obsnum}.o%A"
error="${base}queue/logs/flag_${obsnum}.e%A"

sub="sbatch --begin=now+15 --output=${output} --error=${error} ${depend} ${queue} ${script}"
if [[ ! -z ${tst} ]]
then
    echo "script is ${script}"
    echo "submit via:"
    echo "${sub}"
    exit 0
fi

# submit job
jobid=($(${sub}))
jobid=${jobid[3]}

# rename the err/output files as we now know the jobid
error=`echo ${error} | sed "s/%A/${jobid}/"`
output=`echo ${output} | sed "s/%A/${jobid}/"`

echo "Submitted ${script} as ${jobid}"
