#!/bin/sh

EXIT_SUCCESS=0
EXIT_FAILURE=1
EXIT_NOTIBGS=255

ICEBERGS_MAJOR_VERSION=0
ICEBERGS_MINOR_VERSION=1

echoerr() {
    echo "$@" 1>&2
}

help_usage() {
    echo "usage: ${myCommand} [options] <in_file> [<in_file> ...] <out_file>"
    echo ""
    echo ""
    echo "  -c              Check if file is an icebergs restart file but do"
    echo "                  not combine.  Exit status will be 0 if the file"
    echo "                  is an icebergs restart file, and 255 if not"
    echo "  -h              Display this help and exit"
    echo "  -v              Be verbose in output"
    echo "  -V              Output version information and exit"
}

help() {
    help_usage
    exit $EXIT_SUCCESS
}

get_NumFilesInSet() {
    # Check for the global attribute "NumFilesInSet". If that attribute is
    # found, the function will return the number of files in this set.  If not
    # found, the funcation will return -1. Get the output from ncdump
    local output=$(ncdump -h $1)
    local rtn=-1 # Return value if NumFilesInSet not in file
    numFilesInSet=$(echo $output | sed -n 's/.*:NumFilesInSet *= *\([0-9][0-9]*\) *;.*/\1/p')
    if [ ! x$numFilesInSet = "x" ]
    then
        rtn=$numFilesInSet
    fi
    echo $rtn
}

checkValidFileVersion() {
    # Check for the global attributes "file_format_{major,minor}_version", and
    # return a decimal separated version number.
    local output=$(ncdump -h $1)
    local ncdump_status=$?
    local rtn=$EXIT_NOTIBGS # Value if unable to get the version number
    local major_version=$(echo $output | sed -n 's/.*:file_format_major_version *= *\([0-9][0-9]*\) *;.*/\1/p')
    local minor_version=$(echo $output | sed -n 's/.*:file_format_minor_version *= *\([0-9][0-9]*\) *;.*/\1/p')
    if [ $ncdump_status -eq 0 ]
    then
        if [ ( $major_version -gt $ICEBERGS_MAJOR_VERSION ) -o \
             ( $major_version -eq $ICEBERGS_MAJOR_VERSION -a \
               $minor_version -ge $ICEBERGS_MINOR_VERSION ) ]
        then
            rtn=$EXIT_SUCCESS
        else
            rtn=$EXIT_NOTIBGS
        fi
    else
        rtn=$EXIT_NOTIBGS
    fi
    echo $rtn
}

checkIfIcebergs() {
    # Get the output from ncdump
    local output=$(ncdump -h $1)
    local ncdump_status=$?
    # A file is an icebergs restart file if the file contains the global
    # attributes file_format_major_version, file_format_minor_version, and
    # and NumFilesInSet
    local validFileVersion=$(checkValidFileVersion $1)
    local numFilesInSet=$(get_NumFilesInSet $1)
    if [ $ncdump_status -ne 0 -a \
         $validFileVersion -eq 0 -a \
         $numFilesInSet -gt 0 ]
    then
        rtn=$EXIT_SUCCESS
    else
        rtn=$EXIT_NOTIBGS
    fi
    return $rtn
}

# Name of the command
myCommand=$( basename $0 )
# Other global variables
check=0
verbose=0
version="0.0.1"

while getopts :chvV OPT; do
    case "$OPT" in
        c)
            check=1
            ;;
        h)
            help
            ;;
        v)
            (( verbose++ ))
            ;;
        V)
            echo "${myCommand}: Version ${version}"
            exit $EXIT_SUCCESS
            ;;
        \?)
            echoerr "${myCommand}: unknown option: -${OPTARG}"
            echoerr ""
            help_usage >&2
            exit $EXIT_FAILURE
            ;;
    esac
done
shift $((OPTIND-1))

if [ $check -gt 0 -a $# -ge 1 ]
then
    checkIfIcebergs $1
    exit $?
elif [ $# -lt 2 ]; then
    # Check that at least 2 arguments have been passed in
    echoerr "ERROR: Not enough arguments given"
    help_usage
    exit $EXIT_FAILURE
fi

# Verify the ncdump executable is in PATH and is executable
ncdumpPath=$( which ncdump 2>&1 )
if [ $? -ne 0 ]; then
    echo "ERROR: Required command 'ncdump' is not in PATH."
    exit $EXIT_FAILURE
else
    if [ ! -x $ncdumpPath ]; then
        echoerr "ERROR: Required command 'ncdump' is not executable ($ncdumpPath)."
        exit 1
    fi
fi

# Verify the ncrcat executable is in PATH and is executable (`ncrcat --version` exit status of 0)
which ncrcat > /dev/null 2>&1
if [ $? -ne 0 ]; then
    echoerr "ERROR: NCO command 'ncrcat' is not in PATH."
    exit 1
else
    # Check that ncrcat is executable
    ncrcat --version > /dev/null 2>&1
    if [ $? -ne 0 ]; then
        echoerr "ERROR: NCO command 'ncrcat' is not executable."
        exit 1
    fi
fi

# All of the input files (<iceberg_res_file> [<iceberg_res_file> ...] need to exist.
# <iceberg_out_file> must NOT exist, and should fail if the file exists (unless an -o option given?).

icebergFiles=""  # list of files with icebergs
while [ $# -gt 0 ]; do
    curFile=$1; shift
    if [ $# = 0 ]; then
        # This is the last file listed, and is the output file
        # Exit if this file exists
        if [[ -e $curFile ]]; then
            echoerr "ERROR: out file '$curFile' exists.  Refusing to overwrite file."
            exit 1
        else
            if [ $verbose -gt 1 ]; then
                echoerr "Out file: $curFile"
            fi
            outFile=$curFile
        fi
    else
        # The in files must exist
        if [[ ! -e $curFile ]]; then
            echoerr "ERROR: in file '$curFile' does not exist."
            exit 1
        else
            ncdumpOut=$( ncdump -h $curFile 2>&1 )
            status=$?
            if [ $status -ne 0 ]; then
                echoerr "WARNING: skipping in file '$curFile' as it is NOT a NetCDF formatted file."
            else
                # Check each iceberg_res_file to see if it has icebergs.
                if [ $( echo "$ncdumpOut" | grep 'UNLIMITED' | awk '{gsub(/\(/," ");print $6}' ) -gt 0 ]; then
                    if [ $verbose -gt 1 ]; then
                        echoerr "Using input file $curFile"
                    fi
                    icebergFiles="$icebergFiles $curFile"
                fi
            fi
        fi
    fi
done

# Collect the group of files that do have icebergs, and pass that to ncrcat
if [ "X$icebergFiles" = "X" ]; then
    echoerr "No files to record concatenate."
else
    ncrcat_cmd="ncrcat $icebergFiles $outFile"
    if [ $verbose -gt 0 ]; then
        echo $ncrcat_cmd
    fi
    eval $ncrcat_cmd
fi
