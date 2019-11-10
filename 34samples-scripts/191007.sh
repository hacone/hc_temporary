#!/bin/bash

export TMPDIR=/data/hacone/ncbitmp/sra/

for f in $( ls . | grep .txt$ ); do

        DIR=${f%%.txt}
        mkdir -p $DIR ; cp $f $DIR 

        cd $DIR

                echo "--- in $DIR"
                rm .tmp.to-fetch && touch .tmp.to-fetch
                for j in $( cat *.txt ); do 
                        if [[ ! -e $TMPDIR/$j.sra ]]; then
                                echo $j
                                echo $j >> .tmp.to-fetch
                        fi
                done

                echo $( wc -l .tmp.to-fetch )" runs to fetch"
                echo "cat .tmp.to-fetch | xargs -n 1 -P 12 prefetch -p 1 -v"
                # TODO
                cat .tmp.to-fetch | xargs -n 1 -P 12 prefetch -p 1 --max-size 500G &

                rm .tmp.to-ext && touch .tmp.to-ext
                for j in $( cat *.txt ); do 
                        if [[ ! -e $TMPDIR/$j.fastq.gz ]]; then
                                echo $j
                                echo $j >> .tmp.to-ext
                        fi
                done

                echo $(wc -l .tmp.to-ext)" runs to extract fastq"
                echo "cat .tmp.to-ext | xargs -n 1 -P 12 -I% fastq-dump --skip-technical --gzip -O $TMPDIR $TMPDIR/%.sra"
                # TODO
                cat .tmp.to-ext | xargs -n 1 -P 12 -I% fastq-dump --skip-technical --gzip -O $TMPDIR $TMPDIR/%.sra &

                for j in $( cat *.txt ); do 
                        if [[ -e $TMPDIR/$j.fastq.gz ]]; then
                                if [[ ! -e $TMPDIR/$j.fastq.stats ]]; then
                                        echo "seqkit stats -a $TMPDIR/$j.fastq.gz > $TMPDIR/$j.fastq.stats &"
                                        seqkit stats -a $TMPDIR/$j.fastq.gz > $TMPDIR/$j.fastq.stats &
                                fi
                        fi
                done

        cd ..
done
