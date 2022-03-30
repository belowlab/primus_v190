
FROM ubuntu:20.04
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y perl 
WORKDIR /usr/src
COPY . .

#RUN cpanm --notest -l $PERL_PATH \ Devel::NYTProf

CMD ["mkdir", '-p', "/usr/src/perl"]
ENV PERL_PATH=/usr/src/perl
ENV PERL5LIB=$PERL_PATH:$PERL_PATH/lib/perl5:/usr/src/lib/perl_modules/:/usr/src/lib/perl_modules/PRIMUS/:$PERL5LIB
#RUN cpanm --notest -l $PERL_PATH Devel::NYTProf

ENTRYPOINT ["perl", "./bin/run_PRIMUS.pl"]

# example run: 
#docker build -t primusperl .
#docker run primusperl --file ./example_data/eur_test --genome -o ./output/git3 --keep_inter_files