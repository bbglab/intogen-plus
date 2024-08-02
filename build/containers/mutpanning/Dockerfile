FROM debian:buster-slim

ENV LC_ALL=C.UTF-8

RUN apt-get update && \
    apt-get install -y --no-install-recommends openjdk-11-jre python3 procps && \
    apt-get clean && \
    apt-get remove -y python3-pip && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /mutpanning
# COPY MutPanning.jar /mutpanning/
RUN chmod -R a+rx /mutpanning/

# argument 0: root file, where all the other files can be found
# argument 1: maf file (standard value: root file/MutationsComplete.maf)
# argument 2: sample annotation file (standard value: root file/SamplesComplete.txt)
# argument 3: path to Hg19 folder (standard value root file/Hg19/)
# java -cp MutPanning.jar MutPanning <path_to_outputs> <mutation_file.maf> <samples_file.txt> Hg19/
ENTRYPOINT [ "/usr/bin/java", "-cp", "/mutpanning/MutPanning.jar", "MutPanning" ]