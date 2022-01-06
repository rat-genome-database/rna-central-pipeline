#!/usr/bin/env bash
# shell script to run RNACentral pipeline
. /etc/profile

APPNAME="rna-central-pipeline"
APPDIR=/home/rgddata/pipelines/$APPNAME
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

EMAIL_LIST=mtutaj@mcw.edu
if [ "$SERVER" == "REED" ]; then
  EMAIL_LIST=mtutaj@mcw.edu
fi

cd $APPDIR
java -Dspring.config=$APPDIR/../properties/default_db2.xml \
    -Dlog4j.configurationFile=file://$APPDIR/properties/log4j2.xml \
    -jar lib/${APPNAME}.jar "$@" 2>&1 > $APPDIR/run.log

mailx -s "[$SERVER] RNACentral pipeline run" $EMAIL_LIST < $APPDIR/logs/summary.log

