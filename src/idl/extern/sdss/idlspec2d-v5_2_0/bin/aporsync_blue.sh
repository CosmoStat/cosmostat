#! /bin/sh

#------------------------------------------------------------------------------
# This routine is called by the cron daemon, as set up with "sos_start".
# It syncs files from $RAWDATA_HOST to the local machine (i.e., sos.apo.nmsu.edu).
#
# Note that the directory location on the host machine is hardwired
# as /data/spectro.
#
# We need the executable code "rsync" in the default path
# (e.g., as /usr/bin/rsync).
#
# S. Burles, APO, 4 May 2000
#------------------------------------------------------------------------------

echo "APORSYNC_BLUE: Launched at "`date -u` UID=$UID PPID=$PPID

if [ -z "$RAWDATA_HOST" ] ; then
   echo "Abort: RAWDATA_HOST not set!"
   exit
fi

if [ -z "$RAWDATA_DIR" ] ; then
   echo "Abort: RAWDATA_DIR not set!"
   exit
fi

# This syncs /astrolog/[5-9]???? from the host machine to the local machine,
# exluding the red and guider files.
# Many files might be passed to startapo.sh at once.
rsync -ar --rsh="ssh -c blowfish" \
      --log-format="/data/spectro/%f" \
      --exclude="*-r*" \
      --exclude="*guider*" \
      "${RAWDATA_HOST}:/data/spectro/[5-9]????" $RAWDATA_DIR | startapo.sh 

echo "APORSYNC_BLUE: Finished at "`date -u` UID=$UID PPID=$PPID

