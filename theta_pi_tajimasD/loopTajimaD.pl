#! /usr/bin/perl -w
use strict;
use warnings;
die "perl $0 <position file> <angsd theta output> <number of chromosome>" unless (@ARGV==3);
my $nc=$ARGV[2];
die "chr number should be positive" unless ($nc >0);
my ($a1,$a2,$b1,$b2,$c1,$c2,$e1,$e2);
for my $i (1..$nc-1){
        $a1+=(1/$i);
        $a2+=(1/($i*$i));
}
$b1=($nc+1)/(3*($nc-1));
$b2=(2*($nc*$nc+$nc+3))/(9*$nc*($nc-1));
$c1=$b1-1/$a1;
$c2=$b2-($nc+2)/($a1*$nc)+$a2/($a1*$a1);
$e1=$c1/$a1;
$e2=$c2/($a1*$a1+$a2);
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
                my $seg=$ta*$a1;
                my $top=$tb-$ta;
                my $bot=sqrt($e1*$seg+$e2*$seg*($seg-1));
                my $tajimaD;
                if ($seg==0){
                        $tajimaD=0;
                }else{
                        $tajimaD=sprintf "%.4f", $top/$bot;
                }
                print "$elen\t$tajimaD\n"
        }
        last TNT if ($finish);
        print "$gene\:$st\-$en\t$a[-1]\t$len\t";
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

