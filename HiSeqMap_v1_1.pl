#!/usr/bin/perl -w
use strict ;

use IO::File ;

###################################################################
#
# author: Friedhelm Pfeiffer, MPI of Biochemistry, fpf@biochem.mpg.de
# - version 1.1 (25-SEP-2017) (code unaltered; documentation updated)
# - version 1.0 (15-OCT-2015) (original code as used in the associated paper)
#
# this tool was written to do:
# - read a genome sequence (reference genome)
# - read HiSeq read sequences
# - compare kmers in the genome with those in the reads
# - thus, detect problems in the reference genome
#
# associated publication:
# - Pfeiffer et al,
#   "The complete and fully assembled genome sequence of Aeromonas salmonicida
#    subsp. pectolytica and its comparative analysis with other
#    Aeromonas species: investigation of the mobilome in environmental
#    and pathogenic strains"
#   BMC Genomics (2017) submitted
#
# WARNING NOTES
# - this script/source_code is provided "AS IS"
#   - it was developed with exclusively internal usage in mind
#   - it was not developed for public release but public release was requested;
#     the script has not been improved upon public release
# - documentation is only marginal
#   - if a member of the open source community wants to convert this script
#     into a reasonable piece of software, best in Python instead of Perl,
#      the author offers support (explanations, documentation, co-scripting)
#     - send a request to fpf@biochem.mpg.de
# - this program is "Spaghetti code"!
# - the "fastq file parser" accepts only one specific implementation of
#   that format (the one encountered when analyzing Aeromonas pectinolytica
#   strain 34mel Illumina reads)
#   - WARNING: for reads with quality values consisting exclusively of
#     base characters (A,C,G,T,N), the parser will misunderstand the
#     quality line as a sequence line and will fail!
# - the program is "brute force" with EXCESSIVE memory requirements
#   - even on a 16 GByte main memory computer, the program may fail with an
#     "out of memory" error
#   - due to EXCESSIVE memory usage, other jobs on the computer may fail
# - the program is SLOW (... be patient ... be very patient ...)
#
# - algorithm
#   - the program requires three input parameters:
#     - the genome sequence file (fasta format)
#     - the read sequence file (a special version of the fastq format)
#     - the output file name core
#       - output files with extensions as defined in @outextlist are generated
#   - note: the length of the kmer is fixed in the program
#     (variable $fix_kmerlen)
#   - phase 1:
#     - the genome sequence is read (fasta format)
#       - data are stored in a multidimensional hash
#       - first level key: code (of the replicon, contig)
#       - second level keys: {seq} {title}
#   - phase 2: the Illumina sequencing reads are read
#              (special version of fastq format)
#       - data are stored in a multidimensional hash
#       - first level: paired-read direction (1 or 2)
#       - second level key: "code" (excluding the direction indicator)
#         - codes consist of non-space characters with a terminal "1" or "2"
#           as direction indicator
#       - second level value: sequence
#   - phase 3: check reads
#     - call checkReads()
#     - each "code" (excluding the direction character) must occur twice,
#       once for each direction
#     - the two sequences must have THE SAME LENGTH
#       - this is NOT COMPATIBLE with quality trimming of reads
#   - phase 4 (be patient!): build kmer catalog for genome
#     - call kmerCatalogGenome()
#       - returns two hash references, one for %gkmerhash, one for %revkmers
#     - one-base stepping through the genome
#     - at each step extract the kmer (sequence fragment of kmer length)
#     - three cases
#       - (a) the reverse complement is already in the catalog
#       - (b) neither the fragment nor its reverse complement is in the catalog
#       - (c) the fragment is in the catalog
#     - (case b):
#       - generate an entry for %gkmerhash
#       - subkeys {isDup} "N" {cnt} "1" {readcnt} "0"
#     - (case c): set {isDup} "Y" and increase cnt for kmer
#     - (case a): set {isDup} "Y" and increase cnt for reverse complement
#     - for additional details see the routine-associated documentation
#   - phase 5 (be VERY patient!): map read kmers to genome
#     - call mapReads2Genome()
#     - one-base stepping through each of the sequencing reads
#     - at each step extract the kmer (sequence fragment of kmer length)
#     - five cases
#       - (a) the kmer occurs in the genome
#             (a key in %gkmerhash)
#       - (b) the reverse complement of the kmer occurs in the genome
#             (a key in %revkmers)
#       - (c) the kmer is not in the genome but was seen before
#             (a key in %tmpnewkmers)
#       - (d) the kmer is not in the genome but its reverse complement
#             was seen before (a key in %tmpnewrevkmers)
#       - (e) the kmer is novel and was not seen before
#     - cases a,b are fragments which occur in the genome
#       - genome coverage is affected
#     - cases c,d,e are fragments which do not occur in the genome
#       - these are "novel kmers"
#     - case (e)
#       - generate a key in %tmpnewkmers and initialize to "0" (SHOULD BE "1"!)
#       - compute reverse complement and add to %tmpnewrevkmers
#     - case (c) increase {cnt} for kmer in %tmpnewkmers
#     - case (d) for reverse complement of kmer increase {cnt} in %tmpnewkmers
#     - case (a) increase {readcnt} for kmer in %gkmerhash
#     - case (d) for reverse complement increase {readcnt} in %gkmerhash
#     - note: all kmers that traverse an undefined base ("N") are skipped
#     - note:
#       reads which do not share even a single fragment of kmer length
#       with the genome are written to the "unmatched" output file
#     - for additional details see the routine-associated documentation
#   - phase 6: print novel kmers
#     - call printNovelKmers()
#       - novel kmers are sorted from high to low frequency,
#         alphabetically in case of an identical frequency
#       - quantile data are printed to evaluate if a frequency is relevant
#       - frequent novel kmers point to potential genome sequence errors
#   - phase 7 (be patient!): compute matching statistics
#     - call computeMatchStat()
#       - kmers matching to the genome are sorted from high to low frequency,
#         alphabetically in case of an identical frequency
#       - computes which coverage is exceeded by 75%, 50%, or 25% of the
#         kmers that match to the genome (%quantilhash)
#       - all matching kmers are printed to the "matched" file, sorted
#         from high to low frequency, alphabetically in case of an
#         identical frequency
#         - it is indicated if the kmer is duplicated ("D") or not (".")
#       reads which do not share even a single fragment of kmer length
#       with the genome are written to the "unmatched" output file
#   - phase 8: compute coverage drops and coverage slopes
#     - call detectCoverageDrops()
#     - the algorithm is explained in the routine-associated documentation
#       - in case of zero-division, the max-value "999" is assigned
#     - coverage drops are written to the "covdrop", coverage slopes to the
#       "covslope" file, with sorting from high to low values
#       - the coverage drop is not only printed for the position, but also
#         for adjacent positions (+/- 1 kmer length); high coverage drops
#         which are printed as adjacent to a position are skipped from
#         further printing (same for coverage slopes)
#     - large coverage drops or coverage slopes may point to genome
#       sequence errors
#     - note:
#       at each genome position, the kmers in the vicinity (from -1.5*length
#       to +1.5-times length) are analyzed; to avoid recomputation of kmer
#       sequences, these are stored in the routine-internal hash %tmpkmers
#   - miscellaneous
#     - routine print2log writes parameters (e.g. names of input files)
#       to the .log file
#     - the routine make_complement() computes the reverse complement
#     - the routine fillDuplHash() is no longer used
#
###################################################################

unless (defined $ARGV[2]) {
  print
    "you must specify:\n",
    "- the genome sequence file (fasta format)\n",
    "- the read sequences file (fasta format)\n",
    "- the output file name core\n" ;
  exit ;
}

my $genomfnam = $ARGV[0] ;
my $readfnam = $ARGV[1] ;
my $outfnamcore = $ARGV[2] ;

unless (-e "$genomfnam") {
  print "genome squence file $genomfnam does not exist\n" ;
  exit ;
}
unless (-e "$readfnam") {
  print "read squences file $readfnam does not exist\n" ;
  exit ;
}

my $outfnam ;
my %fho ;
my @outextlist =
  ("log","stat","unmatched","novelkmers","matched","covdrop","covslope") ;
my $filext ;

foreach $filext (@outextlist) {
  $outfnam = "$outfnamcore.$filext" ;
  if (-e "$outfnam") {
    print "file $outfnam already exists\n" ;
    exit ;
  }
  $fho{$filext} = new IO::File (">$outfnam") ;
  if (!defined $fho{$filext}) {
    print "could not open $outfnam\n" ;
    exit ;
  }
}

# reducing $fix_kmerlen from 49 to 45 (hoping to escape memory overflow)
#   did not work
my $fix_kmerlen = 49 ;

my ($hashref1,$hashref2) ;

my %genomhash ;
my %readhash ;
my %gkmerhash ;
my %revkmers ;
my %duphash ;
my %novelkmers ;
my %quantilhash ;

print "phase 1: reading genome\n" ;

if (!defined $fix_kmerlen) {
  print "kmer is nor defined\n" ;
  exit ;
} elsif ($fix_kmerlen =~ /^\d+$/) {
  if ($fix_kmerlen%2 == 1) {
# kmer is ok, continue
  } else {
    print "kmer must be an uneven number: $fix_kmerlen\n" ;
    exit ;
  }
} else {
  print "kmer must be an integer: fix_kmerlen\n" ;
  exit ;
}

print2log(1,"") ;

$hashref1 = parseFastaFile($genomfnam) ;
if (defined $hashref1) {
  %genomhash = %{$hashref1} ;
} else {
  print "could not retrieve genome sequence data\n" ;
  exit ;
}

print "phase 2: reading reads (... be patient ...)\n" ;

$hashref1 = parseFastqFile($readfnam) ;
if (defined $hashref1) {
  %readhash = %{$hashref1} ;
} else {
  print "could not retrieve sequencing reads data\n" ;
  exit ;
}

print "phase 3: checking reads\n" ;

checkReads() ;

print "phase 4: making kmer catalog for genome (... be patient ...)\n" ;

($hashref1,$hashref2) = kmerCatalogGenome() ;
if (defined $hashref1) {
  %gkmerhash = %{$hashref1} ;
} else {
  print "could not retrieve sequencing reads data\n" ;
  exit ;
}
if (defined $hashref2) {
  %revkmers = %{$hashref2} ;
} else {
  print "could not get the reverse kmer hash\n" ;
  exit ;
}

#print "phase 5: filling duplicate hash\n" ;
#fillDuplHash() ;

print "phase 5: kmer mapping of reads to genome (... be very patient ...)\n" ;

$hashref1 = mapReads2Genome() ;
if (defined $hashref1) {
  %novelkmers = %{$hashref1} ;
} else {
  print "could not get the set of novel kmers\n" ;
  exit ;
}

# MEMORY ISSUE: %readhash no longer required
%readhash = () ;

# FILE ISSUE: close "unmatched" and "stat" file ;
$fho{unmatched}->close() ;
$fho{stat}->close() ;

print "phase 6: printing unassigned kmers\n" ;

printNovelKmers() ;

# MEMORY ISSUE: %novelkmers are no longer required
%novelkmers = () ;

# FILE ISSUE: close "novelkmers" file ;
$fho{novelkmers}->close() ;

print "phase 7: computing match statistics (... be patient ...)\n" ;

$hashref1 = computeMatchStat() ;
if (defined $hashref1) {
  %quantilhash = %{$hashref1} ;
} else {
  print "could not get the set of novel kmers\n" ;
  exit ;
}

# FILE ISSUE: close "matched" file ;
$fho{matched}->close() ;

print "phase 8: finding coverage drops and coverage slopes\n" ;

detectCoverageDrops() ;

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
sub print2log {
  my $serial = shift() ;
  my $data = shift() ;

  if (!defined $serial) {
    print "[print2log] serial number must be supplied\n" ;
  }
  if (!defined $data) {
    print "[print2log] data must be supplied\n" ;
  }

# when called first, the input files and kmer length are recorded
  if ($serial == 1) {
    print {$fho{log}} "genome file: $genomfnam\n" ;
    print {$fho{log}} "reads file: $readfnam\n" ;
    print {$fho{log}} "kmer: $fix_kmerlen\n" ;
  } elsif ($serial == 2) {
    print {$fho{log}} "genome size: $data\n" ;
  } elsif ($serial == 3) {
    print {$fho{log}} "number of reads: $data\n" ;
  } elsif ($serial == 4) {
    print {$fho{log}} "genompos dup: $data->{dup}\n" ;
    print {$fho{log}} "genompos adjacent: $data->{adj}\n" ;
    print {$fho{log}} "genompos unique: $data->{uniq}\n" ;
  } elsif ($serial == 5) {
    print {$fho{log}} "total reads: $data->{allreads}\n" ;
    print {$fho{log}} "unmatched reads: $data->{nomatch}\n" ;
    print {$fho{log}} "novel kmers from reads: $data->{novelkmers}\n" ;
  } else {
    print "unexpected serial number $serial\n" ;
    exit ;
  }
}

#----------------------------------------------------------------------------
# definition
# - halflength is ($fix_kmerlen-1)/2, half of the kmer length
# - for discussion on the contig extension by 4*$halflen see below
#
# coverage drops
# - for individual sequence errors (e.g. point mutations, one-base indels)
#   the coverage by Illumina reads must drop
#   - when analyzed by a sliding window, the drop should start halflength before
#     the error and should terminate halflength after the error
#   - accordingly, an "inner" and an "outer" region is used to compute the
#     average coverage
#     - the "inner" region covers halflength upstream plus halflength
#        downstreamof the position under analysis
#     - the "outer" region covers the adjacent halflength upstream/downstream
#   - the ratio of average(outer)/average(inner) is the coverage drop indicator
#   - coverage drop ranges
#     - coverage drop indicators above 900 are set to 900
#     - if there is no matching inner read (div0), coverage drop is set to 999
# - for larger deletions, the coverage will be high on one side while it will
#   be low on the other side
#   - low coverage will start halflength upstream of the deletion point
#   - accrodingly, a "left" and a "right" region is used to compute the
#     average coverage
#     - the "left" region is one kmer, starting halflength upstream of the
#       of the position under analysis
#     - for reasons of symmetry, the "right" region also starts halflength
#       downstream of the position under analysis
#   - the ratio larger/smaller (left/right or right/left) is computed as
#     the coverage slope indicator
#   - coverage slope ranges
#     - coverage slope indicators above 900 are set to 900
#     - if there is no matching read (div0), coverage slope is set to 999
#
# - output data
#   - including description of the %isprinted hash
#   - the most severe coerage drop or slope is printed first
#     - however, the coverage in the surrounding is relevant to know and
#       thus a whole series of coverages is printed
#     - the region starts 2*$halflen upstream and ends 2*$halflen downstream
#   - to avoid printing the same region over and over again, isprinted is
#     set for $halflen positions upstream and downstream of the current
#     position, so that re-printing is suppressed
#     - for positions outside $halflen, isprinted is not set as values in
#       this region are not connected to the current position (kmers outside
#       $halflen will not overlap with the current position)
#   - this makes sure that the most significant position in the area
#     is presented in the output file)
#
# extending the contig by 4*$halflen
# - the covslope analysis averages 2*$halflen positions on either side but
#   starting at a distance of 1*$halflen (to skip kmers which cannot
#   match due to the sequence problem)
# - the outmost positions thus are 3*halflen away but that is the center
#   of the kmer which extends by another $halflen on each side
# - thus, the contig must be extended by 4*$halflen to be able to generate
#   all of the kmers required
#----------------------------------------------------------------------------
sub detectCoverageDrops {
  my $fix_min4print = 4 ;

  my $cnt1 ;
  my $pos ;
  my $code ;
  my $kmer ;
  my $halflen ;
  my $wrkpos ;
  my $tmpseq ;
  my $contiglen ;
  my ($frag5,$frag3) ;
  my ($innerscore,$outerscore) ;
  my ($leftscore,$rightscore) ;
  my ($innerdups,$outerdups) ;
  my ($leftdups,$rightdups) ;
  my $seqshift ;

  my %tmpkmers ;
  my %tmpstat ;
  my %isprinted ;

  $tmpstat{covdrop} = {} ;
  $tmpstat{covslope} = {} ;
  $tmpstat{zerodrop} = {} ;
  $tmpstat{zeroslope} = {} ;

  $halflen = ($fix_kmerlen-1) / 2 ;
# see above for an explanation of 4*$halflen
  $seqshift = 4*$halflen ;

  foreach $code (keys %genomhash) {
    $tmpseq = $genomhash{$code}{seq} ;
    $contiglen = length($tmpseq) ;

    $frag5 = substr($tmpseq,0,$seqshift) ;
    $frag3 = substr($tmpseq,-$seqshift,$seqshift) ;
    $tmpseq = $frag3.$tmpseq.$frag5 ;

my $tmppos ;
my $iii ;

    for ($pos = 1 ; $pos <= $contiglen ; $pos++) {

      if ($pos%500000 == 1) {
        print "having computed coverage for $pos positions for $code\n" ;
      }

$tmppos = $pos + $seqshift ;
      $innerscore = 0 ;
      $outerscore = 0 ;

      $leftscore = 0 ;
      $rightscore = 0 ;

      $innerdups = 0 ;
      $outerdups = 0 ;

      $leftdups = 0 ;
      $rightdups = 0 ;

$iii = 0 ;
      for ($cnt1 = -$halflen ; $cnt1 < 0 ; $cnt1++) {
$iii++ ;
	$wrkpos = $pos+$seqshift+$cnt1 ;
        if (!defined $tmpkmers{$wrkpos}) {
          $kmer = substr($tmpseq,($wrkpos-1)-$halflen,$fix_kmerlen) ;

          if (!defined $gkmerhash{$kmer}) {
            if (defined $revkmers{$kmer}) {
              $kmer = $revkmers{$kmer} ;
	    } else {
              print
                "PROGBUG: the kmer must be defined as %gkmerhash or ",
                "%revkmers: $kmer $wrkpos\n" ;

############# while testing, just use the fist 49 bp of the genome
print "TEST1I1 $iii $tmppos $pos $wrkpos $kmer\n" ;
#              $kmer = "ATGGTGTGGACATCCTTGGGGGCAGTGGTGTATCCACTGTGGGTATCCG" ;
              exit ;
	    }
	  }
          $tmpkmers{$wrkpos} = $kmer ;
	} else {
          $kmer = $tmpkmers{$wrkpos} ;
	}
        $innerscore += $gkmerhash{$kmer}{readcnt} ;
        if ($gkmerhash{$kmer}{isDup} eq "Y") {
          $innerdups++ ;
	}
      }

$iii = 0 ;
      for ($cnt1 = 1 ; $cnt1 <= $halflen ; $cnt1++) {
$iii++ ;
	$wrkpos = $pos+$seqshift+$cnt1 ;
        if (!defined $tmpkmers{$wrkpos}) {
          $kmer = substr($tmpseq,($wrkpos-1)-$halflen,$fix_kmerlen) ;

          if (!defined $gkmerhash{$kmer}) {
            if (defined $revkmers{$kmer}) {
              $kmer = $revkmers{$kmer} ;
	    } else {
              print
                "PROGBUG: the kmer must be defined as %gkmerhash or ",
                "%revkmers: $kmer\n" ;

############# while testing, just use the fist 49 bp of the genome
print "TEST1I2 $iii $tmppos $pos $wrkpos $kmer\n" ;
#              $kmer = "ATGGTGTGGACATCCTTGGGGGCAGTGGTGTATCCACTGTGGGTATCCG" ;
              exit ;
	    }
	  }
          $tmpkmers{$wrkpos} = $kmer ;
	} else {
          $kmer = $tmpkmers{$wrkpos} ;
	}
        $innerscore += $gkmerhash{$kmer}{readcnt} ;
        if ($gkmerhash{$kmer}{isDup} eq "Y") {
          $innerdups++ ;
	}
      }

$iii = 0 ;
      for ($cnt1 = -2*$halflen ; $cnt1 < -$halflen ; $cnt1++) {
$iii++ ;
	$wrkpos = $pos+$seqshift+$cnt1 ;
        if (!defined $tmpkmers{$wrkpos}) {
          $kmer = substr($tmpseq,($wrkpos-1)-$halflen,$fix_kmerlen) ;

          if (!defined $gkmerhash{$kmer}) {
            if (defined $revkmers{$kmer}) {
              $kmer = $revkmers{$kmer} ;
	    } else {
              print
                "PROGBUG: the kmer must be defined as %gkmerhash or ",
                "%revkmers: $kmer\n" ;

############# while testing, just use the fist 49 bp of the genome
print "TEST1O1 $iii $tmppos $pos $wrkpos $kmer\n" ;
#              $kmer = "ATGGTGTGGACATCCTTGGGGGCAGTGGTGTATCCACTGTGGGTATCCG" ;
              exit ;
	    }
	  }
          $tmpkmers{$wrkpos} = $kmer ;
        } else {
          $kmer = $tmpkmers{$wrkpos} ;
	}
        $outerscore += $gkmerhash{$kmer}{readcnt} ;
        if ($gkmerhash{$kmer}{isDup} eq "Y") {
          $outerdups++ ;
	}
      }

$iii = 0 ;
      for ($cnt1 = $halflen+1 ; $cnt1 <= 2*$halflen ; $cnt1++) {
$iii++ ;
	$wrkpos = $pos+$seqshift+$cnt1 ;
        if (!defined $tmpkmers{$wrkpos}) {
          $kmer = substr($tmpseq,($wrkpos-1)-$halflen,$fix_kmerlen) ;

          if (!defined $gkmerhash{$kmer}) {
            if (defined $revkmers{$kmer}) {
              $kmer = $revkmers{$kmer} ;
	    } else {
              print
                "PROGBUG: the kmer must be defined as %gkmerhash or ",
                "%revkmers: $kmer\n" ;

############# while testing, just use the fist 49 bp of the genome
print "TEST1O2 $iii $tmppos $pos $wrkpos $kmer\n" ;
#              $kmer = "ATGGTGTGGACATCCTTGGGGGCAGTGGTGTATCCACTGTGGGTATCCG" ;
              exit ;
	    }
	  }
          $tmpkmers{$wrkpos} = $kmer ;
	} else {
          $kmer = $tmpkmers{$wrkpos} ;
	}
        $outerscore += $gkmerhash{$kmer}{readcnt} ;
        if ($gkmerhash{$kmer}{isDup} eq "Y") {
          $outerdups++ ;
	}
      }
      if ($innerscore == 0) {
	if ($outerscore == 0) {
          $tmpstat{zerodrop}{$pos}++ ;
          $tmpstat{covdrop}{$pos} = -1 ;
        } else {
          $tmpstat{covdrop}{$pos} = 999 ;
        }
      } else {
        $tmpstat{covdrop}{$pos} = $outerscore/$innerscore ;
        $tmpstat{covdrop}{$pos} = sprintf "%.0f",$tmpstat{covdrop}{$pos} ;
        if ($tmpstat{covdrop}{$pos} > 900) {
          $tmpstat{covdrop}{$pos} = 900 ;
	}
      }

$iii = 0 ;
      for ($cnt1 = -3*$halflen ; $cnt1 < -$halflen ; $cnt1++) {
$iii++ ;
	$wrkpos = $pos+$seqshift+$cnt1 ;
        if (!defined $tmpkmers{$wrkpos}) {
          $kmer = substr($tmpseq,($wrkpos-1)-$halflen,$fix_kmerlen) ;

          if (!defined $gkmerhash{$kmer}) {
            if (defined $revkmers{$kmer}) {
              $kmer = $revkmers{$kmer} ;
	    } else {
              print
                "PROGBUG: the kmer must be defined as %gkmerhash or ",
                "%revkmers: $kmer\n" ;

############# while testing, just use the fist 49 bp of the genome
print "TEST2L $iii $tmppos $pos $wrkpos $kmer\n" ;
#              $kmer = "ATGGTGTGGACATCCTTGGGGGCAGTGGTGTATCCACTGTGGGTATCCG" ;
              exit ;
	    }
	  }
          $tmpkmers{$wrkpos} = $kmer ;
	} else {
          $kmer = $tmpkmers{$wrkpos} ;
	}
        $leftscore += $gkmerhash{$kmer}{readcnt} ;
        if ($gkmerhash{$kmer}{isDup} eq "Y") {
          $leftdups++ ;
	}
      }

$iii = 0 ;
      for ($cnt1 = $halflen+1 ; $cnt1 <= 3*$halflen ; $cnt1++) {
$iii++ ;
	$wrkpos = $pos+$seqshift+$cnt1 ;
        if (!defined $tmpkmers{$wrkpos}) {
          $kmer = substr($tmpseq,($wrkpos-1)-$halflen,$fix_kmerlen) ;

          if (!defined $gkmerhash{$kmer}) {
            if (defined $revkmers{$kmer}) {
              $kmer = $revkmers{$kmer} ;
	    } else {
              print
                "PROGBUG: the kmer must be defined as %gkmerhash or ",
                "%revkmers: $kmer\n" ;
############# while testing, just use the fist 49 bp of the genome
print "TEST2R $iii $tmppos $pos $wrkpos $kmer\n" ;
#              $kmer = "ATGGTGTGGACATCCTTGGGGGCAGTGGTGTATCCACTGTGGGTATCCG" ;
              exit ;
	    }
	  }
          $tmpkmers{$wrkpos} = $kmer ;
	} else {
          $kmer = $tmpkmers{$wrkpos} ;
	}
        $rightscore += $gkmerhash{$kmer}{readcnt} ;
        if ($gkmerhash{$kmer}{isDup} eq "Y") {
          $rightdups++ ;
	}
      }
      if ($leftscore >= $rightscore) {
	if ($rightscore == 0) {
          if ($leftscore == 0) {
            $tmpstat{zeroslope}{$pos}++ ;
            $tmpstat{covslope}{$pos} = -1 ;
          } else {
            $tmpstat{covslope}{$pos} = 999 ;
          }
        } else {
          $tmpstat{covslope}{$pos} = $leftscore/$rightscore ;
          $tmpstat{covslope}{$pos} = sprintf "%.0f",$tmpstat{covslope}{$pos} ;
          if ($tmpstat{covslope}{$pos} > 900) {
            $tmpstat{covslope}{$pos} = 900 ;
	  }
	}
      } else {
	if ($leftscore == 0) {
          $tmpstat{covslope}{$pos} = 999 ;
        } else {
          $tmpstat{covslope}{$pos} = $rightscore/$leftscore ;
          $tmpstat{covslope}{$pos} = sprintf "%.0f",$tmpstat{covslope}{$pos} ;
          if ($tmpstat{covslope}{$pos} > 900) {
            $tmpstat{covslope}{$pos} = 900 ;
	  }
	}
      }

#if ($pos >= 4) { exit ; }
    }
  }

  print "done with data collection, writing results to file\n" ;

  %isprinted = () ;

# first print positions in regions not having any matching kmer
  foreach $pos (sort {$a <=> $b} keys %{$tmpstat{covdrop}}) {
    if (defined $tmpstat{zerodrop}{$pos}) {
      print {$fho{covdrop}}
        "$pos: no matching kmer for $fix_kmerlen on either side\n" ;
      $isprinted{$pos}++ ;
    }
  }

# then print the positions with a high coverage drop
  foreach $pos (sort {$tmpstat{covdrop}{$b} <=> $tmpstat{covdrop}{$a} or
                      $a <=> $b}
                keys %{$tmpstat{covdrop}}) {

# skip positions which have already been in focus with a higher covdrop
    if (defined $isprinted{$pos}) {
      next ;
    }

# stop printing when covdrop decreases below $fix_min4print
    if ($tmpstat{covdrop}{$pos} < $fix_min4print) {
      last ;
    }

    print {$fho{covdrop}} "\ncovdrop for $pos ($tmpstat{covdrop}{$pos})\n" ;

    for ($cnt1 = -2*$halflen; $cnt1 <= 2*$halflen; $cnt1++) {
      $wrkpos = $pos+$cnt1 ;
      if ($wrkpos < 1) {
	print {$fho{covdrop}}
          "below pos_1: $wrkpos $pos $cnt1\n" ;
      } elsif ($wrkpos > $contiglen) {
	print {$fho{covdrop}}
          "beyond contiglen $contiglen: $wrkpos $pos $cnt1\n" ;
      } else {
        print {$fho{covdrop}}
          "$wrkpos ($pos $cnt1): $tmpstat{covdrop}{$wrkpos}\n" ;
        if ($cnt1 >= -$halflen and $cnt1 <= $halflen) {
          $isprinted{$wrkpos}++ ;
	}
      }
    }
  }

  %isprinted = () ;

# first print positions in regions not having any matching kmer
  foreach $pos (sort {$a <=> $b} keys %{$tmpstat{covslope}}) {
    if (defined $tmpstat{zeroslope}{$pos}) {
      print {$fho{covslope}}
        "$pos: no matching kmer for more than $fix_kmerlen on either side\n" ;
      $isprinted{$pos}++ ;
    }
  }

# then print the positions with a high coverage slope
  foreach $pos (sort {$tmpstat{covslope}{$b} <=> $tmpstat{covslope}{$a} or
                      $a <=> $b}
                keys %{$tmpstat{covslope}}) {

# skip positions which have already been in focus with a higher covslope
    if (defined $isprinted{$pos}) {
      next ;
    }

# stop printing when covslope decreases below $fix_min4print
    if ($tmpstat{covslope}{$pos} < $fix_min4print) {
      last ;
    }

    print {$fho{covslope}} "\ncovslope for $pos ($tmpstat{covslope}{$pos})\n" ;

    for ($cnt1 = -2*$halflen; $cnt1 <= 2*$halflen; $cnt1++) {
      $wrkpos = $pos+$cnt1 ;
      if ($wrkpos < 1) {
	print {$fho{covslope}}
          "below pos_1: $wrkpos $pos $cnt1\n" ;
      } elsif ($wrkpos > $contiglen) {
	print {$fho{covslope}}
          "beyond contiglen $contiglen: $wrkpos $pos $cnt1\n" ;
      } else {
        print {$fho{covslope}}
          "$wrkpos ($pos $cnt1): $tmpstat{covslope}{$wrkpos}\n" ;
        if ($cnt1 >= -$halflen and $cnt1 <= $halflen) {
          $isprinted{$wrkpos}++ ;
	}
      }
    }
  }
}

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
sub printNovelKmers {
  my $cnt ;
  my $kmer ;

  my %tmphash ;

  print "reshuffling novelkmers by frequency\n" ;

  foreach $kmer (keys %novelkmers) {
    $cnt = $novelkmers{$kmer} ;

    if (!defined $tmphash{$cnt}) {
      $tmphash{$cnt} = {} ;
    }
    delete $novelkmers{$kmer} ;
    $tmphash{$cnt}{$kmer}++ ;
  }

  print "now printing novelkmers\n" ;

  foreach $cnt (sort {$b <=> $a} keys %tmphash) {
    foreach $kmer (sort keys %{$tmphash{$cnt}}) {
      print {$fho{novelkmers}} "$cnt $kmer\n" ;
    }
  }
}

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
sub computeMatchStat {
  my $cnt1 ;
  my $kmer ;
  my $totalcnt ;
  my $dupchar ;

  my %tmpstat ;

  $cnt1 = 0 ;
  $totalcnt = keys %gkmerhash ;
  foreach $kmer (sort {$gkmerhash{$b}{readcnt} <=> $gkmerhash{$a}{readcnt} or
                       $a cmp $b}
                 keys %gkmerhash) {
    if ($gkmerhash{$kmer}{isDup} eq "Y") {
# ok, continue
    } elsif ($gkmerhash{$kmer}{isDup} eq "N") {
      $cnt1++ ;

      if (!defined $tmpstat{med25} and $cnt1/$totalcnt >= (1-0.25)) {
        $tmpstat{med25} = $gkmerhash{$kmer}{readcnt} ;
        print {$fho{matched}} "!! quantil_25: $tmpstat{med25}\n" ;
        print {$fho{log}} "!! quantil_25: $tmpstat{med25}\n" ;
      } elsif (!defined $tmpstat{med50} and $cnt1/$totalcnt >= (1-0.50)) {
        $tmpstat{med50} = $gkmerhash{$kmer}{readcnt} ;
        print {$fho{matched}} "!! quantil_50: $tmpstat{med50}\n" ;
        print {$fho{log}} "!! quantil_50: $tmpstat{med50}\n" ;
      } elsif (!defined $tmpstat{med75} and $cnt1/$totalcnt >= (1-0.75)) {
        $tmpstat{med75} = $gkmerhash{$kmer}{readcnt} ;
        print {$fho{matched}} "!! quantil_75: $tmpstat{med75}\n" ;
        print {$fho{log}} "!! quantil_75: $tmpstat{med75}\n" ;
      }
    } else {
      print "invalid isDup value for $kmer: $gkmerhash{$kmer}{isDup}\n" ;
      exit ;
    }

    if ($gkmerhash{$kmer}{isDup} eq "Y") {
      $dupchar = "D" ;
    } elsif ($gkmerhash{$kmer}{isDup} eq "N") {
      $dupchar = "." ;
    } else {
      print "invalid isDup value for $kmer: $gkmerhash{$kmer}{isDup}\n" ;
      exit ;
    }

    print {$fho{matched}} "$gkmerhash{$kmer}{readcnt} $dupchar $kmer\n" ;
  }

  return \%tmpstat ;
}

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
sub fillDuplHash {
  my $cnt1 ;
  my $pos ;
  my $code ;
  my $contiglen ;
  my $tmpseq ;
  my $kmer ;
  my $wrkkmer ;
  my ($frag5, $frag3) ;
  my $halflen ;
  my $tmppos ;

  my %tmpstat ;

  $halflen = ($fix_kmerlen-1) / 2 ;

  foreach $code (keys %genomhash) {
    if (defined $duphash{$code}) {
      print "duplicate contig code $code\n" ;
      exit ;
    }
    $duphash{$code} = {} ;

    $tmpseq = $genomhash{$code}{seq} ;
    $contiglen = length($tmpseq) ;

    $frag5 = substr($tmpseq,0,$halflen) ;
    $frag3 = substr($tmpseq,-$halflen,$halflen) ;
    $tmpseq = $frag3.$tmpseq.$frag5 ;

    for ($pos = 1 ; $pos <= $contiglen ; $pos++) {

    if ($pos%1000000 == 1) {
      print "having checked $pos positions for $code\n" ;
    }

# $pos-1: positions are 1-based, substr works 0-based
    $kmer = substr($tmpseq,($pos-1)-$halflen,$fix_kmerlen) ;

    if (!defined $gkmerhash{$kmer} and !defined $revkmers{$kmer}) {
      print "all kmers must be defined by now but $kmer is not; ABORTING\n" ;
      exit ;
    }

    if (defined $gkmerhash{$kmer}) {
      $wrkkmer = $kmer ;
    } elsif (defined $revkmers{$kmer}) {
      $wrkkmer = $revkmers{$kmer} ;
    }

    if ($gkmerhash{$wrkkmer}{isDup} eq "N") {
      if (!defined $duphash{$code}{$pos}) {
        $duphash{$code}{$pos} = "-" ;
      }
    } elsif ($gkmerhash{$wrkkmer}{isDup} eq "Y") {
	for ($cnt1 = -$halflen ; $cnt1 <= $halflen ; $cnt1++) {
          $tmppos = $pos+$cnt1 ;
          $duphash{$code}{$tmppos} = "dup" ;
        }
	for ($cnt1 = -($fix_kmerlen-1) ; $cnt1 <= -($halflen-1) ; $cnt1++) {
          $tmppos = $pos+$cnt1 ;
          if ($duphash{$code}{$tmppos} eq "-") {
            $duphash{$code}{$tmppos} = "adj" ;
	  }
        }
	for ($cnt1 = $halflen+1 ; $cnt1 <= ($fix_kmerlen-1) ; $cnt1++) {
          $tmppos = $pos+$cnt1 ;
          if (!defined $duphash{$code}{$tmppos}) {
            $duphash{$code}{$tmppos} = "adj" ;
          } elsif ($duphash{$code}{$tmppos} eq "-") {
            $duphash{$code}{$tmppos} = "adj" ;
	  }
        }
      } else {
        print "invalid isDup value for $kmer: $gkmerhash{$wrkkmer}{isDup}\n" ;
        exit ;
      }
    }
  }

  $tmpstat{dummy} = 0 ;
  $tmpstat{adj} = 0 ;
  $tmpstat{dup} = 0 ;
  foreach $code (keys %genomhash) {
    $tmpseq = $genomhash{$code}{seq} ;
    $contiglen = length($tmpseq) ;

    for ($pos = 1 ; $pos <= $contiglen ; $pos++) {
      if (!defined $duphash{$code}{$pos}) {
          print
            "duphash not defined for $code $pos: ",
            "$duphash{$code}{$pos}\n" ;
          exit ;
      } elsif ($duphash{$code}{$pos} eq "-") {
        $tmpstat{uniq}++ ;
      } elsif ($duphash{$code}{$pos} eq "adj") {
        $tmpstat{adj}++ ;
      } elsif ($duphash{$code}{$pos} eq "dup") {
        $tmpstat{dup}++ ;
      } else {
          print
            "unexpected term in duphash for $code $pos: ",
            "$duphash{$code}{$pos}\n" ;
          exit ;
      }
    }
  }
  print2log(4,\%tmpstat) ;
}

#----------------------------------------------------------------------------
# some info on algorithm
# - collected data
#   - number of matching kmers (hitcnt)
#   - first matching kmer
#   - last matching kmer
#
# - reads with hitcnt zero
#   - these do not map to the genome, write to .nonmapped
#
# - reads with at least one hit
#   - if kmers hits: increase {readcnt} in the corresponding kmer
#     - this is checked against kmers and revdups!
#   - if kmers do not hit: add to %nohitkmerhash
#----------------------------------------------------------------------------
sub mapReads2Genome {
  my $pos ;
  my $term ;
  my $code ;
  my $endnr ;
  my $seq ;
  my $seqlen ;
  my ($kmer,$revkmer) ;
  my $totalcnt ;
  my $codcnt ;
  my $allcnt ;
  my $cnt4read ;

  my %tmpnewkmers ;
  my %tmpnewrevkmers ;
  my %tmpstat ;

  $allcnt = (keys %{$readhash{1}}) + (keys %{$readhash{2}}) ;
  $codcnt = 0 ;

  $tmpstat{allreads} = $allcnt ;
  $tmpstat{nomatch} = 0 ;

  foreach $endnr (1,2) {
    foreach $code (keys %{$readhash{$endnr}}) {
      $codcnt++ ;
      if ($codcnt%500000 == 1) {
        print "having done $codcnt of $allcnt\n" ;
      }

      $seq = $readhash{$endnr}{$code} ;
      $seqlen = length $seq ;

      $cnt4read = 0 ;
      for ($pos = 1 ; $pos <= $seqlen-$fix_kmerlen ; $pos++) {
        $kmer = substr($seq,$pos-1,$fix_kmerlen) ;

# skip all kmers with an undefined base (N)
        if ($kmer =~ /N/) {
	  next ;
	}
        if (defined $gkmerhash{$kmer} or defined $revkmers{$kmer}) {
          $cnt4read++ ;
	}
      }
      $tmpstat{hitcnt}{$cnt4read}++ ;

      if ($cnt4read == 0) {
        print {$fho{unmatched}} ">$code $code\n$readhash{$endnr}{$code}\n" ;
        $tmpstat{nomatch}++ ;
      } else {
        for ($pos = 1 ; $pos <= $seqlen-$fix_kmerlen ; $pos++) {
          $kmer = substr($seq,$pos-1,$fix_kmerlen) ;

# skip all kmers with an undefined base (N)
          if ($kmer =~ /N/) {
            next ;
          }
          if (defined $gkmerhash{$kmer}) {
            $gkmerhash{$kmer}{readcnt}++ ;
          } elsif (defined $revkmers{$kmer}) {
            $revkmer = $revkmers{$kmer} ;
            $gkmerhash{$revkmer}{readcnt}++ ;
          } elsif (defined $tmpnewkmers{$kmer}) {
            $tmpnewkmers{$kmer}++ ;
          } elsif (defined $tmpnewrevkmers{$kmer}) {
            $revkmer = $tmpnewrevkmers{$kmer} ;
            $tmpnewkmers{$revkmer}++ ;
          } else {
            $tmpnewkmers{$kmer} = 0 ;
            $revkmer = make_complement($kmer) ;
            $tmpnewrevkmers{$revkmer} = $kmer ;
            $tmpstat{novelkmers}++ ;
	  }
	}
      }
    }
    print "having done $codcnt of $allcnt at end of set $endnr\n" ;
  }

# at this point, %tmpnewrevkmers are no longer needed
  %tmpnewrevkmers = () ;

#  print {$fho{stat}} "total number of reads: $tmpstat{allreads}\n" ;
#  print {$fho{stat}} "unmatched reads: $tmpstat{nomatch}\n" ;

#  print {$fho{stat}} "\n" ;
#  print {$fho{stat}} "hitcnt distribution\n" ;
#  $term = "hitcnt" ;
#  foreach $pos (sort {$a <=> $b} keys %{$tmpstat{$term}}) {
#    print {$fho{stat}} "$term $pos $tmpstat{$term}{$pos}\n" ;
#  }

  print2log(5,\%tmpstat) ;

  return \%tmpnewkmers ;
}

#----------------------------------------------------------------------------
# some info on algorithm
# - kmer characteristic and halflen
#   - kmers must be of uneven size
#   - halflen is half of (kmer-1)
#
# - handling of circular genomes
#   - genomes are assumed to be circular. Thus, a special twist is required
#     to collect kmers at the point of ring opening
#   - the 3pr halflen bases are appended upstream of the 5pr end
#   - the 5pr halflen bases are appended at the 3pr end
#
# - reverse complements
#   - genome duplications may be in the same or opposite orientation
#   - opposite orientation duplicates lead to kmers which are
#     reverse complements of each other
#   - such reverse complements are excluded from the kmer list and the
#     data are accumulated with the partner kmer
#   - the duplication flag is set of the partner kmer
#   - as reads can be from the FWD and REV strand, all kmer reverse complements
#     are produced and kept in a separate hash
#----------------------------------------------------------------------------
sub kmerCatalogGenome {
  my $code ;
  my $pos ;
  my $contigcnt ;
  my $tmpseq ;
  my $revkmer ;
  my $contiglen ;
  my $genomlen ;
  my $kmer ;
  my ($frag5, $frag3) ;
  my $halflen ;

  my %tmphash ;
  my %tmpcompl ;

  $halflen = ($fix_kmerlen-1) / 2 ;

  $genomlen = 0 ;
  foreach $code (keys %genomhash) {
    $tmpseq = $genomhash{$code}{seq} ;
    $contiglen = length($tmpseq) ;
    $genomlen += $contiglen ;

    $frag5 = substr($tmpseq,0,$halflen) ;
    $frag3 = substr($tmpseq,-$halflen,$halflen) ;
    $tmpseq = $frag3.$tmpseq.$frag5 ;

    for ($pos = 1 ; $pos <= $contiglen ; $pos++) {
# $pos-1: positions are 1-based, substr works 0-based

      if ($pos%500000 == 1) {
        print "working on $code at $pos (of $contiglen bp)\n" ;
      }

      $kmer = substr($tmpseq,$pos-1,$fix_kmerlen) ;

      if (defined $tmpcompl{$kmer}) {
        $revkmer = $tmpcompl{$kmer} ;
        $tmphash{$revkmer}{isDup} = "Y" ;
        $tmphash{$revkmer}{cnt}++ ;
      } elsif (!defined $tmphash{$kmer}) {
        $tmphash{$kmer} = {} ;
        $tmphash{$kmer}{isDup} = "N" ;
        $tmphash{$kmer}{cnt} = 1 ;
        $tmphash{$kmer}{readcnt} = 0 ;

        $revkmer = make_complement($kmer) ;
        $tmpcompl{$revkmer} = $kmer ;
      } else {
        $tmphash{$kmer}{isDup} = "Y" ;
        $tmphash{$kmer}{cnt}++ ;
      }
    }
  }

  $contigcnt = keys %genomhash ;
  print2log(2,"$genomlen bp in $contigcnt contigs") ;

  return (\%tmphash,\%tmpcompl) ;
}

#----------------------------------------------------------------------------
# a number of checks is performed
# - actually, many of the properties checked are not mandatory for this script
#
# the following checks are performed
# - there must always be read pairs
# - the length of the read pair must be identical
#----------------------------------------------------------------------------
sub checkReads {
  my %fix_endnrs = (
    1 => 2,
    2 => 1
  ) ;

  my $code ;
  my ($endnr,$otherend) ;
  my ($len1,$len2) ;

  foreach $endnr (sort keys %fix_endnrs) {
    $otherend = $fix_endnrs{$endnr} ;
    foreach $code (keys %{$readhash{$endnr}}) {
      if (!defined $readhash{$otherend}{$code}) {
        print "only end $endnr defined for $code\n" ;
        exit ;
      }
    }
  }

  foreach $code (keys %{$readhash{1}}) {
    $len1 = length($readhash{1}{$code}) ;
    $len2 = length($readhash{2}{$code}) ;
    if ($len1 != $len2) {
      print "length discrepancy for $code: $len1 $len2\n" ;
      exit ;
    }
  }

# print number of reads to log file
  $len1 = (keys %{$readhash{1}}) + (keys %{$readhash{2}}) ;
  print2log(3,$len1) ;
}

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
sub parseFastaFile {
  my $wrkfnam = shift() ;

  my $INFIL ;

  my $line ;
  my $code ;
  my $title ;
  my $wrkseq ;

  my %tmphash ;

  if (!defined $wrkfnam) {
    print "[parseFastaFile] the file name must be supplied!\n" ;
    exit ;
  }

  $INFIL = new IO::File ("$wrkfnam") ;
  if (!defined $INFIL) {
    print "could not open $wrkfnam\n" ;
    exit ;
  }

  while ($line = $INFIL->getline) {
    chomp $line ;

    if ($line =~ /^>(\S+)(\s+(.+))?/) {
      $code = $1 ;
      $title = $3 ;

      if (defined $tmphash{$code}) {
        print "duplicate code $code\n" ;
        exit ;
      }
      $tmphash{$code} = {} ;
      $tmphash{$code}{seq} = "" ;

      if (defined $title) {
        $tmphash{$code}{title} = $title ;
      } else {
        $tmphash{$code}{title} = "DUMMY" ;
      }

    } else {
      if (!defined $code) {
        print "this seems not to be a fasta file: $line\n" ;
        exit ;
      }
      $tmphash{$code}{seq} .= $line ;
    }
  }
  $INFIL->close ;

  foreach $code (sort keys %tmphash) {
    $tmphash{$code}{seq} =~ s/\s//g ;

    if ($tmphash{$code}{seq} !~ /^[ACGT]+$/) {
      print "invalid sequence for $code: '$tmphash{$code}{seq}'\n" ;
      exit ;
    }
  }

  return \%tmphash ;
}

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
sub parseFastqFile {
  my $wrkfnam = shift() ;

  my $INFIL ;

  my $line ;
  my $code ;
  my $seq ;
  my $endnr ;
  my $lintyp ;
  my $seqlen ;
  my $codecnt ;

  my %tmphash ;

  $tmphash{1} = {} ;
  $tmphash{2} = {} ;

  if (!defined $wrkfnam) {
    print "[parseFastqFile] the file name must be supplied!\n" ;
    exit ;
  }

  $INFIL = new IO::File ("$wrkfnam") ;
  if (!defined $INFIL) {
    print "could not open $wrkfnam\n" ;
    exit ;
  }

  $codecnt = 0 ;
  $lintyp = 1 ;
  $seqlen = 0 ;
  while ($line = $INFIL->getline) {
    chomp $line ;

    if ($line =~ /^\@(\S+)\/([12])$/) {
      $code = $1 ;
      $endnr = $2 ;

      if ($lintyp != 1) {
        print "line type mismatch for codeline: $lintyp\n" ;
        exit ;
      }
      $lintyp++ ;

      if (!defined $endnr) {
        print "PROGBUG: end number is not defined\n" ;
        exit ;
      } elsif ($endnr == 1 or $endnr == 2) {
# ok, continue
      } else {
        print "invalid end number for $code: $endnr\n" ;
        exit ;
      }

      if (defined $tmphash{$endnr}{$code}) {
        print "duplicate code end $code $endnr\n" ;
        exit ;
      }
      $tmphash{$endnr}{$code} = "" ;

      $codecnt++ ;
      if ($codecnt%500000 == 1) {
        print "having read $codecnt entries\n" ;
      }
    } elsif ($line =~ /^([ACGTN]+)$/) {
      $seq = $1 ;

      if ($lintyp != 2) {
        print "line type mismatch for seqline: $lintyp\n" ;
        exit ;
      }
      $lintyp++ ;

      $seqlen = length $seq ;

      if ($tmphash{$endnr}{$code} eq "") {
	  $tmphash{$endnr}{$code} = $seq ;
      } else {
        print "there is already a sequence for $code\n" ;
        exit ;
      }
    } elsif ($line eq "+") {
      if ($lintyp != 3) {
        print "line type mismatch for plusline: $lintyp\n" ;
        exit ;
      }

      $lintyp++ ;
    } else {
      if ($lintyp != 4) {
        print "line type mismatch for xline: $lintyp '$line'\n" ;
        exit ;
      }
      $lintyp = 1 ;

      if (length($line) ne $seqlen) {
        print "length mismatch in xline, length should be $seqlen\n" ;
        exit ;
      }
    }
  }

  return \%tmphash ;
}

###############################################################################
# make_complement
# - this routine takes a sequencde and returns the complement
###############################################################################
sub make_complement {
  my $routine = "make_complement" ;

  my $seq = shift ;
  my $seqlen = length ($seq) ;
  my $seqchar ;
  my $complement ;
  my $pos ;

  my %compl = (
    "A" => "T",
    "G" => "C",
    "C" => "G",
    "T" => "A",
    "N" => "N",
    "R" => "Y",
    "Y" => "R",
    "S" => "S",
    "W" => "W",
    "B" => "V",
    "H" => "D",
    "D" => "H",
    "V" => "B",
    "X" => "X"
  ) ;

  $complement = "" ;

  for ($pos = $seqlen-1 ; $pos >= 0 ; $pos--) {
    $seqchar = substr $seq,$pos,1 ;
    unless (defined $compl{$seqchar}) {
      die "program error! do not know the complement of $seqchar\n" ;
    }
    $complement .= $compl{$seqchar} ;
  }
  return ($complement) ;
}
