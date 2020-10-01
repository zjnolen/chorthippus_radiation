! usr/bin/perl -w
use strict;
use warnings;
die "perl $0 <position file> <angsd theta output>" unless (@ARGV==2);
open IN, "$ARGV[0]" or die $!;
open FS, "gzip -dc $ARGV[1]|" || die $!;
my $col=$ARGV[2];
my (@alpha, @beta);
my $finish=0;
<FS>;
chomp(my $fl=<FS>);
my @check=split /\t/,$fl;
my $endp=$check[1];
die "the second line don't have asked info" unless ($check[0]=~/chr1/ and @check > 4);
seek (FS, 0, 0);
TNT:while (my $po=<IN>){
        next TNT if ($po=~/^derivedFileNanme/);
        chomp ($po);
        my @a=split /\t/,$po;
        my $gene=(split /\:/,$a[0])[-1];
        my @b=split /\:/,$a[1];
        my $chr=$b[0];
        my $st=(split /\-/,$b[1])[0];
        my $en=(split /\-/,$b[1])[1];
        next TNT if ($en < $endp and $finish==0);
        my $len=$en-$st+1;
#       print "$gene\:$st\-$en\t"
        if(@alpha){
                my $ta; my $tb;
                warn "no match for two array\n" unless (@alpha == @beta);
                my $elen=$#alpha+$finish;
                for my $ia (0..($#alpha-1+$finish)){
                        my $na=shift @alpha;
                        my $ena=exp($na);
                        $ta+=$ena;
                }
                for my $ib (0..($#beta-1+$finish)){
                        my $nb=shift @beta;
                        my $enb=exp($nb);
                        $tb+=$enb;
                }
                print "$elen\t$ta\t$tb\n"
        }
        last TNT if ($finish);
        print "$gene\:$st\-$en\t$a[-1]\t$len\t";
#       for my $i (0..($#alpha-1)){
#               shift @alpha;
#               shift @beta;
#       }
#       last TNT if ($finish);
        TTT:while(my $fst=<FS>){
                next TTT if ($fst=~/^#/);
                if (eof(FS)){
                        $finish=1;
                }
                chomp($fst);
                my @c=split /\t/,$fst;
                if ($c[1]>$en){
                        push @alpha,$c[2];
                        push @beta,$c[3];
                        $endp=$c[1];
                        next TNT;
                }else{
                        push @alpha,$c[2];
                        push @beta,$c[3];
                }
        }
}
close IN;
close FS;

