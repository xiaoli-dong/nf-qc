#!/usr/bin/env perl
use strict;
use File::Basename;
use Data::Dumper;

my $AUTHOR = 'Xiaoli Dong <xdong@ucalgary.ca>';
my $VERSION = "0.0.1";
my $quiet = 0;
my $EXE = basename($0);

my(@Options, $seqtype);
setOptions();

@ARGV or die "Please provide a directory name containing the fastq format sequence files!\n";
my $seqdir = $ARGV[0];

opendir(my $dh, $seqdir) || die "Can't open $seqdir: $!";

my %seqfile_prefixnames = ();
my %stats = ();

while (readdir $dh) {

    next if /^\.\.?$/;
    if(/(\S+?)(_1|_2|_R1|_R2|\.R1|\.R2)/){
	push(@{$seqfile_prefixnames{$1}}, $_);
    }
}
closedir $dh;

#print Dumper(\%seqfile_prefixnames);



for my $prefixname (keys %seqfile_prefixnames){

    my @seqfile_names = @{$seqfile_prefixnames{$prefixname}};
    
    my @seq_files = map {"$seqdir/$_"} @seqfile_names;

    fq($prefixname, $seqtype, \@seq_files);
}

my @header = ("Seqfile", "seqtype", "Reads", "Yield", "GeeCee", "MinLen", "AvgLen", "MaxLen", "AvgQual", "ErrQual", "Ambiguous");

print join(",", @header), "\n";
for my $key (keys %stats){
    print join(",", @{$stats{$key}}), "\n";
    
}


sub fq{
    my($prefixname, $seqtype, $seq_files) = @_;
    
    #min_len: 35; max_len: 151; avg_len: 147.63; 6 distinct quality values
    #POS     #bases  %A      %C      %G      %T      %N      avgQ    errQ    %low    %high
    #ALL     134985892       21.7    28.2    28.6    21.5    0.0     32.1    22.9    11.3 >
    #1       914379  19.4    23.9    43.8    13.0    0.0     30.9    25.7    4.3     95.7
    #2       914379  18.3    28.6    18.0    35.1    0.0     31.1    26.3    3.6     96.4
    #....
    #150     785526  22.8    26.1    31.3    19.8    0.0     27.3    19.8    23.8    76.2
    #151     530460  31.4    0.0     44.7    23.9    0.0     23.9    17.9    38.3    61.7

    my $filename_str = join(" ", @$seq_files);
    msg("processed data from $prefixname");
    
    my %stat;
    my $cmd = "cat $filename_str  | seqtk fqchk -q0 -";
    msg("running command: $cmd");
    open my $IN, '-|', $cmd or err("could not run command: $cmd");

    while (<$IN>) {
	if (m/^min_len/) {
	    s/\s//g;
	    for my $pair (split m';') {
		my($k,$v) = split m':', $pair;
		$stat{$k} = $v if $v;
	    }
	}
	elsif (m/^ALL/) {
	    my @x = split ' ';
	    $stat{total_bp} = $x[1];
	    $stat{gee_cee} = $x[3] + $x[4];
	    $stat{avg_qual} = $x[7];
	    $stat{err_qual} = $x[8];
	    $stat{ambig_bp_pc} = $x[6];
	}
	elsif (m/^1\s+(\d+)\b/) {
	    $stat{num_reads} = $1;
	}
    }
    msg("processed", $stat{num_reads}, "reads from $prefixname dataset.");
    my @out = ();
        
    push(@out, $prefixname);
    push(@out, $seqtype);

    if($seqtype eq "paired"){
	push(@out, int($stat{num_reads}/2));
    }
    else{
	push(@out, $stat{num_reads});
    }
    push(@out,  $stat{total_bp});
    push(@out,  $stat{gee_cee});
    push(@out,  $stat{min_len});
    push(@out,  int($stat{avg_len}));
    push(@out,  $stat{max_len});
    push(@out, $stat{avg_qual});
    push(@out, $stat{err_qual});
    push(@out, $stat{ambig_bp_pc});
    #print join(",", @out), "\n";
    
    $stats{$prefixname} = \@out;
    
}

sub msg { print STDERR "[$EXE] @_\n" unless $quiet; }
sub err { $quiet=0; msg("ERROR:", @_); exit(-1); }


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
    use Getopt::Long;
    
    @Options = (
        {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
        {OPT=>"quiet!",  VAR=>\$quiet, DEFAULT=>0, DESC=>"Quiet mode: no progress output"},
	{OPT=>"seqtype", VAR=>\$seqtype, DEFAULT=>'paired', DESC=>"paired|single, default is \"paired\""},
        );

    (!@ARGV) && (usage());
    
    &GetOptions(map {$_->{OPT}, $_->{VAR}} grep { ref } @Options) || usage();
    
    # Now setup default values.
    foreach (@Options) {
	if (ref $_ && defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
	    ${$_->{VAR}} = $_->{DEFAULT};
	}
    }
}

sub usage {
    print STDERR
	"Name:\n  ", $EXE, " $VERSION by $AUTHOR\n",
	"\nSynopsis:\n  Generate basic stats for the sequence files contained in the directory\n",
	"\nUsage:\n  $EXE [options] <directory path>\n";
    
    foreach (@Options) {
	if (ref) {
	    my $def = defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
	    $def = ($def ? ' (default OFF)' : '(default ON)') if $_->{OPT} =~ m/!$/;
	    my $opt = $_->{OPT};
	    $opt =~ s/!$//;
	    $opt =~ s/=s$/ [X]/;
	    $opt =~ s/=i$/ [N]/;
	    $opt =~ s/=f$/ [n.n]/;
	    printf STDERR "  --%-15s %s%s\n", $opt, $_->{DESC}, $def;
	}
	else {
	    print STDERR "\n$_\n";
	}
    }
    exit(1);
}

