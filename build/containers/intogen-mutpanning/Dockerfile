FROM debian:buster-slim AS build

ENV LC_ALL=C.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends --no-install-suggests \
        openjdk-11-jdk-headless=11.0.23+9-1~deb10u1 \
        git=1:2.20.1-2+deb10u9 && \
    rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/genepattern/MutPanning.git --branch v2 && \
    cd MutPanning && \
    javac -cp bin/commons-math3-3.6.1.jar:bin/jdistlib-0.4.5-bin.jar \
          src/AffinityCount_Cosmic.java src/AffinityCount.java src/AlignHG19.java \
          src/CBASE_Solutions.java src/ClusteringEntity.java src/ClusteringPanCancer.java \
          src/ComputeMutationRateClusters_Entities.java src/ComputeSignificance.java \
          src/CountDestructiveMutations.java src/Filter_Step1.java src/Filter_Step2.java \
          src/Filter_Step3.java src/MutPanning.java src/ReformatCBASE.java -d bin && \
    cd bin && \
    printf "Manifest-Version: 1.0\nMain-Class: MutPanning\nClass-Path: commons-math3-3.6.1.jar jdistlib-0.4.5-bin.jar\n" >MANIFEST.MF && \
    jar cfm MutPanning.jar MANIFEST.MF ./*.class && \
    rm ./*.class MANIFEST.MF

FROM debian:buster-slim AS final

ENV LC_ALL=C.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends --no-install-suggests \
        openjdk-11-jre-headless=11.0.23+9-1~deb10u1 \
        python3.7-minimal=3.7.3-2+deb10u7 \
        procps=2:3.3.15-2 && \
    apt-get remove -y python3-pip && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /mutpanning
COPY --from=build /MutPanning/bin/*.jar /mutpanning/

# argument 0: root file, where all the other files can be found
# argument 1: maf file (standard value: root file/MutationsComplete.maf)
# argument 2: sample annotation file (standard value: root file/SamplesComplete.txt)
# argument 3: path to Hg19 folder (standard value root file/Hg19/)
# java -cp MutPanning.jar MutPanning <path_to_outputs> <mutation_file.maf> <samples_file.txt> Hg19/
ENTRYPOINT [ "/usr/bin/java", "-cp", "/mutpanning/MutPanning.jar", "MutPanning" ]
