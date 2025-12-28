#!/bin/bash

OUTPUT="nf_final_results.csv"

echo "process,pod_name,node,start_time,end_time,duration_sec" > $OUTPUT

# Query K8s for all Nextflow pods that have finished (Succeeded or Failed)
# Use jq to parse the dates and calculate duration on the fly
kubectl get pods -l nextflow.io/app=nextflow,role!=launcher -o json | jq -r '
  .items[] |
  select(.status.phase == "Succeeded" or .status.phase == "Failed") |
  [
    .metadata.labels["nextflow.io/processName"],
    .metadata.name,
    .spec.nodeName,
    .status.startTime,
    .status.containerStatuses[0].state.terminated.finishedAt,
    ((.status.containerStatuses[0].state.terminated.finishedAt | fromdateiso8601) -
     (.status.startTime | fromdateiso8601))
  ] | @csv' >> $OUTPUT

echo "Done! Data exported to $OUTPUT"
