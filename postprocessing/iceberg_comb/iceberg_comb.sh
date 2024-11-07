#!/bin/sh

script_name=$(basename "$0")

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

find_command() {
    # Check if the command is usable at the given path.  If not, check
    # if a similar command exists in the user's PATH.
    #
    # Arguments:
    #    $1 : Standard command name
    #    $2 : Full path to command
    local rtn=""
    command -v $2 > /dev/null 2>&1
    if [ $? -eq 0 ]
    then
        rtn=$2
    else
        # check for the command in PATH
        command -v $1 > /dev/null 2>&1
        if [ $? -eq 0 ]
        then
            # found something in PATH
            rtn=$(command -v $1)
        fi
    fi
    echo "$rtn"
}

get_NumFilesInSet() {
    # Check for the global attribute "NumFilesInSet". If that attribute is
    # found, the function will return the number of files in this set.  If not
    # found, the funcation will return -1. Get the output from ncdump
    local output=$($ncdump -h $1)
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
    local output=$($ncdump -h $1)
    local ncdump_status=$?
    local rtn=$EXIT_NOTIBGS # Value if unable to get the version number
    local major_version=$(echo $output | sed -n 's/.*:file_format_major_version *= *\([0-9][0-9]*\) *;.*/\1/p')
    local minor_version=$(echo $output | sed -n 's/.*:file_format_minor_version *= *\([0-9][0-9]*\) *;.*/\1/p')
    if [ $ncdump_status -eq 0 ]
    then
        if [  $major_version -gt $ICEBERGS_MAJOR_VERSION -o \
              $major_version -eq $ICEBERGS_MAJOR_VERSION -a \
              $minor_version -ge $ICEBERGS_MINOR_VERSION ]
        then
            rtn=$EXIT_SUCCESS
        else
            rtn=$EXIT_NOTIBGS
        fi
    else
        rtn=$EXIT_NOTIBGS
    fi
    return $rtn
}

checkIfIcebergs() {
    # Default return value
    local rtn=$EXIT_NOTIBGS
    # Use ncdump to determine if the file is a netcdf file
    $ncdump -k $1 > /dev/null 2>&1
    if [ $? -eq $EXIT_SUCCESS ]
    then
        # A file is an icebergs restart file if the file contains the global
        # attributes file_format_major_version, file_format_minor_version, and
        # and NumFilesInSet
        checkValidFileVersion $1
        local validFileVersion=$?
        local numFilesInSet=$(get_NumFilesInSet $1)
        if [ $validFileVersion -eq 0 -a \
             $numFilesInSet -gt 0 ]
        then
            rtn=$EXIT_SUCCESS
        else
            rtn=$EXIT_NOTIBGS
        fi
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

# Check that required commands are still available, if not available
# at the given path, fallback to the something in the path if
# available.
ncdump=$(find_command "ncdump" "${NCDUMP:-/Users/seth.underwood/opt/homebrew/bin/ncdump}")
if [ -z $ncdump ]
then
    echoerr "$script_name: Unable to find command 'ncdump'"
    exit $EXIT_FAILURE
fi
if [ $verbose -gt 0 ]
then
    echoerr "debug1: Using ncdump command '$ncdump'"
fi

ncrcat=$(find_command "ncrcat" "${NCRCAT:-/Users/seth.underwood/opt/homebrew/bin/ncrcat}")
if [ -z $ncrcat ]
then
    echoerr "$script_name: Unable to find command 'ncrcat'"
    exit $EXIT_FAILURE
fi
if [ $verbose -gt 0 ]
then
    echoerr "debug1: Using ncrcat command '$ncrcat'"
fi

ncatted=$(find_command "ncatted" ${NCATTED:-/Users/seth.underwood/opt/homebrew/bin/ncatted})
if [ -z $ncatted ]
then
    echoerr "$script_name: Unable to find command 'ncatted'"
    exit $EXIT_FAILURE
fi
if [ $verbose -gt 0 ]
then
    echoerr "debug1: Using ncrcat command '$ncrcat'"
fi

# The check option overrides the normal procedure.  The user
# should only check one file at a time -- this simplifies the
# reporting
if [ $check -gt 0 ]
then
    if [ $# -gt 1 ]
    then
        echoerr "$script_name: warning: Only checking the first file '$1'"
    fi
    if [ ! -e $1 ]
    then
        echoerr "$script_name: $1: No such file or directory"
        exit $EXIT_FAILURE
    fi
    checkIfIcebergs $1
    exit $?
fi

# The user should have supplied at least three items, at least
# two files to combine, and the output file.
if [ $# -lt 3 ]
then
    echoerr "$script_name: Not enough arguments given"
    exit $EXIT_FAILURE
fi

# All but the last argument should be valid icebergs restart files.
# The last argument file should not exist, and should be written
# to a directory the user can write to.
all_files="$@"
out_file=${all_files##* }
in_files=${all_files% *}
if [ $verbose -gt 1 ]
then
    echoerr "debug2: Received $(echo "$in_files" | wc -w) input files"
fi

file_check=$EXIT_SUCCESS
for f in $in_files
do
    # If any of the files do not exist, or
    if [ ! -e $f ]
    then
        echoerr "$script_name: $f: No such file or directory"
        file_check=$EXIT_FAILURE
    else
        checkIfIcebergs $f
        if [ $? -ne $EXIT_SUCCESS ]
        then
            echoerr "$script_name: $f: Not an icebergs restart file"
            file_check=$EXIT_FAILURE
        fi
    fi
done
if [ $file_check -ne $EXIT_SUCCESS ]
then
    echoerr "$script_name: Unable to combine iceberg files"
    exit $EXIT_FAILURE
fi
# Check if the user supplied all the expected number of input files
# The number of files should match the global attribute NumFilesInSet
# We need to also check that all the files have the same value for
# NumFilesInSet.
# Since we are here, we already know that the files are valid
# iceberg restart files
file_check=$EXIT_SUCCESS
num_setFiles_0=$(get_NumFilesInSet $1)
if [ $verbose -gt 1 ]
then
    echoerr "debug2: File '$1' expects $num_setFiles_0 files in set"
fi
for f in ${in_files#* }
do
    num_setFiles_n=$(get_NumFilesInSet $f)
    if [ $num_setFiles_0 -ne $num_setFiles_n ]
    then
        file_check=$EXIT_FAILURE
    fi
done
if [ $file_check -ne $EXIT_SUCCESS ]
then
    echoerr "$script_name: Incorrect or inconsistent number of files in set"
    exit $EXIT_FAILURE
fi
# Verify that all the set files are given
if [ $num_setFiles_0 -ne $(echo "$in_files" | wc -w) ]
then
    echoerr "$script_name: Expected $num_setFiles_0 files, but got $(echo "$in_files" | wc -w)"
    exit $EXIT_FAILURE
fi

# Check if the output file exists.
if [ -e $out_file ]
then
    echoerr "$script_name: $f: Output file exists"
    exit $EXIT_FAILURE
fi


# Run ncrcat on the files
ncrcat_cmd="$ncrcat $@"
if [ $verbose -gt 0 ]
then
    echoerr "debug1: Running '$ncrcat_cmd'"
fi
eval $ncrcat_cmd
if [ $? -ne 0 ]
then
    echoerr "$script_name: Error during executing '$ncrcat_cmd'"
    exit $EXIT_FAILURE
fi

ncatted_cmd="$ncatted -a NumFilesInSet,global,d,, $out_file"
if [ $verbose -gt 0 ]
then
    echoerr "debug1: Running '$ncatted_cmd'"
fi
eval $ncatted_cmd
if [ $? -ne 0 ]
then
    echoerr "$script_name: Error during executing '$ncatted_cmd'"
    exit $EXIT_FAILURE
fi

exit $EXIT_SUCCESS
