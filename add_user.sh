#!/bin/bash

# Default values
DEFAULT_PASSWORD="LiuLab666!"
DEFAULT_SHELL="/bin/bash"
DEFAULT_GROUP=""
DEFAULT_HOME_PREFIX="/home"
DEFAULT_EXTRA_GROUPS=""
DEFAULT_EXPIRE_DATE=""
IS_ROOT=false

# Usage information
usage() {
  echo "Usage: $0 -u <username> [-p <password>] [-d <home_directory>] [-s <shell>] [-g <group>] [-G <additional_groups>] [-e <expire_date>] [-r]"
  exit 1
}

# Parse command line arguments
while getopts ":u:p:d:s:g:G:e:r" opt; do
  case $opt in
    u) USERNAME="$OPTARG"
    ;;
    p) PASSWORD="$OPTARG"
    ;;
    d) HOMEDIR="$OPTARG"
    ;;
    s) SHELL="$OPTARG"
    ;;
    g) GROUP="$OPTARG"
    ;;
    G) EXTRA_GROUPS="$OPTARG"
    ;;
    e) EXPIRE_DATE="$OPTARG"
    ;;
    r) IS_ROOT=true
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
        usage
    ;;
    :) echo "Option -$OPTARG requires an argument." >&2
        usage
    ;;
  esac
done

# Check if the username has been provided
if [ -z "$USERNAME" ]; then
  echo "Username is required."
  usage
fi

# Set default values
PASSWORD=${PASSWORD:-$DEFAULT_PASSWORD}
SHELL=${SHELL:-$DEFAULT_SHELL}
GROUP=${GROUP:-$DEFAULT_GROUP}
HOMEDIR=${HOMEDIR:-"$DEFAULT_HOME_PREFIX/$USERNAME"}
EXTRA_GROUPS=${EXTRA_GROUPS:-$DEFAULT_EXTRA_GROUPS}
EXPIRE_DATE=${EXPIRE_DATE:-$DEFAULT_EXPIRE_DATE}

# Create user command
CMD="sudo useradd -m -d $HOMEDIR -s $SHELL"

# If a user group has been provided, add it to the command
[ -n "$GROUP" ] && CMD+=" -g $GROUP"

# If additional groups have been provided, add them to the command
[ -n "$EXTRA_GROUPS" ] && CMD+=" -G $EXTRA_GROUPS"

# If an expiry date has been provided, add it to the command
[ -n "$EXPIRE_DATE" ] && CMD+=" -e $EXPIRE_DATE"

# If the user should have root privileges, add them to the sudo group
$IS_ROOT && CMD+=" -G sudo"

# Add username
CMD+=" $USERNAME"

# Execute the command to create the user
eval "$CMD"

# Set the password
echo "$USERNAME:$PASSWORD" | sudo chpasswd

# Output the result
echo "User $USERNAME created successfully with the provided information."
