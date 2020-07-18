#!usr/bin/perl
#use strict;

my $path = $ARGV[0];                ## Argument: /path/to/file.xmfa ## MAUVE FORMAT ##
my $ref = $ARGV[1];                 ## Number of order of the sequence used as reference ##
my @path = split /\//, $path;
my $file = pop @path;
$path = join "\/", @path;

print "\nxmfa2fas\.pl\tversion 1.0\tLSB Oct2013\n\n";
print "Program called as\: $path/$file $ref\n\n";
if ($ref){
    print "Reordering as\: sequence $ref\n\n";
}else{
    print "Reordering as\: default \(order in XMFA file\)\n\n";
}

my @names;                          ## Get names of genomes ##
my $numb = 0;                       ## Get number of genomes ##
my %posits;
my %gennames;
my %posct;
my $check = 0;

my ($name, $seq, $num, $lenfile, $sig);
my %rev;
open XMFA, "$path/$file";
while(<XMFA>){
    my $l = $_;
    $lenfile .= $l;
    if ($l =~ /^#Sequence\d+File[\t|\s+](.*)/){
        my $n = $1;
        chomp $n;
        if ($n =~ /.*\\(.*)\r/){
            $n = $1;
        }
        $n =~ s/\s/\_/g;
        $n =~ s/\\/\_/;
        $n =~ s/\r//;
        $n =~ s/\.\.\///g;
        $n =~ s/\..*//;
        $numb++;
        $posits{$numb} = $n;
        $gennames{$numb} = $n;
        $posct{$numb} = 0;
    }elsif ($l =~ /^#BackboneFile/){
        $check = 1;   
    }elsif ($l =~ /^> (.*)\:(\d+)\-\d+ ([\+|\-]) .*\n/ && $check == 1){
        if ($seq){
            $seq =~ s/\=//g;
            if ($posct{$num} == 0){
                $posits{$num} = {$name => $seq};                    # Create a Hash of Hashes!!! (HoH) #
            }else{
                $posits{$num}{$name} = $seq;                        # Add information to a HoH #
            }
            $posct{$num}++;                                         # Control of the number of breaks #
        }
        undef $name, undef $seq, undef $sig;
        $name = $2;
        $num = $1;
        $sig = $3;
    }elsif ($l !~ /^> ((.*)\:\d+\-\d+ [\+|\-] .*)\n/ && $check == 1) {
        if ($name || $name == 0){
            chomp, $_ =~ s/\r//;
            if ($sig eq '+'){
                $seq .= $_;
            }elsif ($sig eq '-'){
                $seq .= $_;
                my $r = "$num\:$name";
                if ($num == $ref){
                    $rev{$r} = $r;
                }
            }
        }
    }
}
close XMFA;
$seq =~ s/\=//g;                                                    ## Get last sequence ##
if ($posct{$numb} == 0){
    $posits{$num} = {$name => $seq};                    
}else{
    $posits{$num}{$name} = $seq;
}
$posct{$num}++;
undef $name, undef $seq, undef $sig;

my @lenfile = split /\=/, $lenfile;                                 ## Split the whole file by segments ##
pop @lenfile;
my %bone;
my $seg = 0;
my @lenfile2 = @lenfile;
do{
    $seg++;
    my $w = $seg;
    my $x = shift @lenfile2;
    my @x = split /\>/, $x;
    shift @x;
    foreach (@x){
        if (/(\d+\:\d+)\-/){
            $bone{$w} .= "\|$1\|";                                    ## Get num of genome + ini position ##
        }
    }
}until scalar @lenfile2 == 0;

my $c = scalar keys %bone;

print "\nNb_segments\: $c\n";

foreach (keys %rev){                                               ## TRANSLATE OTHER SEQUENCES FROM SEGMENT WITH REVERSE REFERENCE ##
    my $r = $_;
    foreach (keys %bone){
        my $b = $_;
        if ($bone{$b} =~ /$r\|/){
            my @b = split /\|/, $bone{$b};
            foreach (@b){
                if (/(\d+)\:(\d+)/){
                    my $g = $1;
                    my $p = $2;
                    my $s = $posits{$g}->{$p};
                    my $rc = reverse($s);
                    $rc =~ tr/ACGTacgt/TGCAtgca/;
                    $posits{$g}->{$p} = $rc;
                }
            }
        }
    }
}

my %seqs;

if ($ref){                                                          ## Get sequence $ARGV[1] as reference for ordering ##
    my @order = keys $posits{$ref};
    @order = sort {$a <=> $b} @order;
    my $incl = "";
    do{
        my $ord = shift @order;
        foreach (keys %bone){                                       ## Get segment number and other genomes from that segment ##
            my $k = $_;
            if ($bone{$k} =~ /\|$ref\:$ord\|/){
                $incl .= "\|$k\|";                                  ## Save the included segment ##
                my @conts = split /\|/, $bone{$k};
                @conts = grep { !/^$/ } @conts;
                my ($ini, $gen);
                my $ch = "";
                my $len;
                foreach(@conts){                                                        
                        if (/^(\d+)\:(\d+)/){
                            $gen = $1;
                            
                            $ini = $2;
                        }
                        $ch .= "\_$gen";
                        foreach (keys %{$posits{$gen}}){
                            if ($_ == $ini){
                                if(!$seqs{$gen}){
                                    $seqs{$gen} = $posits{$gen}->{$ini};
                                    $len = $posits{$gen}->{$ini};
                                }else{
                                    $seqs{$gen} .= $posits{$gen}->{$ini};
                                    $len = $posits{$gen}->{$ini};
                                }
                                last;
                            }
                        }
                }
                my @left;
                for (my $z = 1; $z<=$numb; $z++){                                   ## Get absent genomes in current segment and add '-' ##
                    if ($ch !~ /\_$z/){
                        push @left, $z;
                    }
                }
                my $ln = length $len;
                foreach (@left){
                    my $lf = $_;
                    for (my $x = 1; $x <= $ln; $x++){   
                        $seqs{$lf} .= '-';
                    }
                }
            last;
            }
        }
    }until scalar @order == 0;
    my %left;
    my $count_seg = 0;
    for (my $j=1; $j<=$seg; $j++){             ## Add extra segments at the end ##
        if ($incl !~ /\|$j\|/){                ## If they had not already been included...##
            $count_seg++;
            my $conts = $bone{$j};
            my @conts = split /\|/, $conts;
            @conts = grep { !/^$/ } @conts;
            my ($ini, $gen);
            my $ch = "";
            my $len;
            foreach(@conts){
                if (/^(\d+)\:(\d+)/){
                    $gen = $1;
                    $ini = $2;
                }
                $ch .= "\_$gen";
                my $nam = "$gen\:$ini\_seg$j";
                $left{$nam} = $posits{$gen}->{$ini};
                foreach (keys %{$posits{$gen}}){
                    if ($_ == $ini){
                        if(!$seqs{$gen}){
                            $seqs{$gen} = $posits{$gen}->{$ini};
                            $len = $posits{$gen}->{$ini};
                        }else{
                            $seqs{$gen} .= $posits{$gen}->{$ini};
                            $len = $posits{$gen}->{$ini};
                        }
                        last;
                    }
                }
            }
            my @left;
            for (my $z = 1; $z<=$numb; $z++){                                   ## Get absent genomes in current segment and add '-' ##
                if ($ch !~ /\_$z/){
                    push @left, $z;
                }
            }
            my $ln = length $len;
            foreach (@left){
                my $lf = $_;
                for (my $x = 1; $x <= $ln; $x++){   
                $seqs{$lf} .= '-';
                }
            }
        }
    }
    print "\nNb_left_segments\: $count_seg\n";
    
}else{                                                              ## Without reference, segment order as XMFA ##
    for (my $j=1; $j<=$seg; $j++){
        my $conts = $bone{$j};
        my @conts = split /\|/, $conts;
        @conts = grep { !/^$/ } @conts;
        my ($ini, $gen);
        my $ch = "";
        my $len;
        foreach(@conts){
            if (/^(\d+)\:(\d+)/){
                $gen = $1;
                $ini = $2;
            }
            $ch .= "\_$gen";
            foreach (keys %{$posits{$gen}}){
                if ($_ == $ini){
                    if(!$seqs{$gen}){
                        $seqs{$gen} = $posits{$gen}->{$ini};
                        $len = $posits{$gen}->{$ini};
                    }else{
                        $seqs{$gen} .= $posits{$gen}->{$ini};
                        $len = $posits{$gen}->{$ini};
                    }
                    last;
                }
            }
        }
        my @left;
        for (my $z = 1; $z<=$numb; $z++){                                   ## Get absent genomes in current segment and add '-' ##
            if ($ch !~ /\_$z/){
                push @left, $z;
            }
        }
        my $ln = length $len;
        foreach (@left){
            my $lf = $_;
            for (my $x = 1; $x <= $ln; $x++){   
            $seqs{$lf} .= '-';
            }
        }
    }
}

$file =~ s/\..*//;
open OUT, ">$path/$file\_$num\_ref$ref\.fas";
foreach (keys %gennames){
    my $k = $gennames{$_};
    $k =~ s/\s//g;
    print OUT "\>$k\n$seqs{$_}\n";
}
close OUT;


exit;