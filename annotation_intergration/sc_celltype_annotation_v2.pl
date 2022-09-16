use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw/max min sum maxstr minstr shuffle/;
my $ver="1.1";

my %opts;
GetOptions(\%opts,"i=s", "d=s", "o=s", "h" );

if(!defined($opts{i}) || !defined($opts{d}) || defined($opts{h})){
	print <<"	Usage End.";
	Description:
                
		Version: $ver   

	Usage:

		-i           annotation  file                              <infile>     must be given
		-d           database dir file                             <infile>     must be given
		-o           output file                                   <outfile>    option given
		-h           Help document
	
	Usage End.

	exit;
}

my $annotation = $opts{i} ;
my $db = $opts{d} ;
my $outfile = defined $opts{o} ? $opts{o} : "final_anotation.txt";

open O, ">$outfile";
#my $db = shift;
my %anno = &read_database($db);
=head
for my $k (sort keys %anno){
	print "$k\t$anno{$k}\n";
}
=cut

my $in = "$annotation";
open IN, "$in";
print O "cluster\tSingleR\tscCATCH\tclustermole\tcelltype\n";
open O1, ">cluster2cellType.txt";
print O1 "cluster\tcelltype\n";
while(<IN>){
	chomp;
	next if $. == 1;
	my ($cluster, $singler, $sccatch, $clustermole)=split/\t/, $_;
	my @anno_info;
	my @a;
	for(split/\|/, $clustermole){
		push @a, "$_($anno{$_})";
		push @anno_info, $anno{$_};
	}
	my @b;
	my @sccatch_info;
	#if ($sccatch =~/\,/){
	#$sccatch =~ s/\s//g;
	for(split/\,/, $sccatch){
		$_ =~ s/^\s+|\s+$//g;
		my $count = $_ =~ tr/ / /;
		my $t = (split/\s+/, $_)[1];
		$_ = (split/\s+/, $_, 2)[1] if $t =~/^[A-Z]/;;
		$anno{$_} = "undefined" if (!exists $anno{$_});
		push @b, "$_($anno{$_})";
		push @anno_info, $anno{$_};
	}
	#}
	$anno{$singler} = "undefined" if (!exists $anno{$singler});
	$anno{$sccatch} = "undefined" if (!exists $anno{$sccatch});
	
	push @anno_info, $anno{$singler};
	push @anno_info, $anno{$sccatch};
	
	#print O "$cluster\t$singler($anno{$singler})\t$sccatch($anno{$sccatch})\t";
	my @clustermole_anno = join "|", @a;
	#print O join "|", @a;
	my @sccatch_anno = join "|", @b;
	
	print O "$cluster\t$singler($anno{$singler})\t", @sccatch_anno, "\t", @clustermole_anno, "\t";
	#print O "\t";
	
	print O getMax(@anno_info);
	print O "\n";

	print O1 "$cluster\t";
	print O1 getMax(@anno_info);
	print O1 "\n";
}

sub read_database(){
	#my $dir = shift;
	my ($dir) = @_;
	#print "$dir\n";
	my @file = <$dir/*/*txt>;
	my %f;
	for my $file (@file){
		open IN, "$file";
		while(<IN>){
			chomp;
			my ($raw_type, $new_type) = (split/\t/, $_)[0, 1];
			$f{$raw_type} = $new_type;
			#print "$file\t$raw_type\t$new_type\n";
		}
		close IN;
	}
	return %f;
}

sub getMax {
    my (%temp, @ret);
    my $max = 0;
    foreach ( @_ ) {
        my $count = $temp{$_}++;
        $max = $count if $count > $max; 
    }
    foreach ( keys %temp ) {
        push @ret, $_ if $temp{$_} == $max + 1
    }    
    return (sort @ret)[0];
}

##
sub getmax{
        my ($input)=@_;
        my %occ;
        for(@$input){
                $occ{$_} +=1;
        }
        my $maxocc = (sort {$b <=> $a} values %occ)[0];
        while(my($k, $v) = each %occ){
                if($v == $maxocc ){
                        print "$k ";
                }
        }
}
