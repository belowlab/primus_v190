
FROM ubuntu:20.04
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    perl \
    plink1.9 \
    r-cran-devtools
WORKDIR /usr/src/
COPY . .

#RUN cpanm --notest -l $PERL_PATH \ Devel::NYTProf
RUN ln -s /bin/plink1.9 /bin/plink
CMD ["mkdir", '-p', "/usr/src/perl"]
ENV PERL_PATH=/usr/src/perl
ENV PERL5LIB=$PERL_PATH:$PERL_PATH/lib/perl5:/usr/src/lib/perl_modules/:/usr/src/lib/perl_modules/PRIMUS/:$PERL5LIB
#RUN cpanm --notest -l $PERL_PATH Devel::NYTProf

WORKDIR /usr/src/bin
ENTRYPOINT ["perl", "./run_PRIMUS.pl"]

# example run: 
#docker build -t primus1 .
#docker run primus1 --file ../example_data/eur_test --genome -o ./test --keep_inter_files --no_PR