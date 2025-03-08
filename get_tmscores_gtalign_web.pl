#!/usr/bin/env perl
BEGIN {$^W=1}

use strict;

## NOTE: all directories and executables must exist!
my $WRKDIR = "/home/mindaugas/projects/data/gtalign-web-benchmark";
my $TMalign ="/home/mindaugas/install/TM-align/TMalign";
my $TMscore = "/home/mindaugas/install/TM-score--2022-2-27/TMscore_from_alignment";
my $getchainprog = "/data/installed-software/comer-ws-backend/bin/getchain.py";

my $QRYDIR = "$WRKDIR/queries";
my $RFNDIR = "/data-SSD/databases/wwpdb/mmCIF";
##my $RDIR = "$WRKDIR/gtalign_web_speed13_prescore04";
my $RDIR = "$WRKDIR/gtalign_web_speed9_prescore04";
my $WDIR = "$WRKDIR/tmp_wrk_dir";

## NOTE: 1st argument is a query structure filename
my $QIDfle = $ARGV[0];
my $GDT = ($ARGV[1] && ($ARGV[1] eq "gdt"))? 1: 0;

$TMalign = $TMscore if $GDT;

exit 0 unless $QIDfle;

my $QID = $QIDfle; $QID =~ s/.pdb$//;

exit 0 unless $QID;

my $QSTR = substr($QID,0,4);
my $QCHN = substr($QID,5,1);
my $MAXNHITS = 1000;

die "ERROR: Query $QID not found." unless -f "$QRYDIR/$QIDfle";

my $resfile = glob("$RDIR/$QID*.json");
my $cnt = 0;
my ($rfn, $qaln, $raln, $hits) = ("","","",0);
my ($qlen, $rlen, $qstart, $qend, $rstart, $rend) = (0,0,0,0,0,0);
my ($gtm1, $gtm2, $tm1, $tm2, $gdt) = (-1,-1,-1,-1,-1);

open(F, $resfile) || die "ERROR: Failed to open results: $resfile";

while(<F>) {
  if(eof(F) || /{"hit_record":\s+{/)
  {
    $hits = 1 if /{"hit_record":\s+{/;
    if($rfn && $qaln && $raln)
    {
      my $wfile = "$WDIR/${QID}_${rfn}.aln";
      my $tmfile = "$WDIR/${QID}_${rfn}.tm";
      my $rstr = substr($rfn,0,4);
      my $rchn = substr($rfn,5);
      my $rmdl = substr($rfn,1,2);
      my $rfnfile = "$RFNDIR/$rmdl/${rstr}.cif.gz";
      my $qalnlen = length($qaln); my $ralnlen = length($raln);
      my $qnnb = $qstart - 1; my $qnne = $qlen - $qend;
      my $rnnb = $rstart - 1; my $rnne = $rlen - $rend;
      $qaln = "X"x$qnnb . $qaln . "X"x$qnne;
      $raln = "-"x$qnnb . $raln . "-"x$qnne;
      $qaln = "-"x$rnnb . $qaln . "-"x$rnne;
      $raln = "X"x$rnnb . $raln . "X"x$rnne;
      unless(-f $rfnfile) {
        print(STDERR "WARNING: Reference $rfn not found: $rfnfile\n");
      }
      else {
        last if $MAXNHITS <= $cnt++;
        if("$QSTR$QCHN" ne "$rstr$rchn")
        {
          my $trgfile = "$WDIR/${QID}_${rfn}.pdb";
          if(GetTrgFile($rfnfile, $trgfile, $rchn))
          {
            open(W, ">$wfile") || die "ERROR: Failed to open file for writing: $wfile";
            print(W ">$QSTR $QCHN $QID\n$qaln\n>$rstr $rchn $rfn\n$raln\n");
            close(W);
            my $prog = "$TMalign $QRYDIR/$QIDfle $trgfile -I $wfile -het 1 >$tmfile";
            if(system($prog) != 0) {
              print(STDERR "ERROR: Failed to execute: $prog\n");
            } else {
              open(W, $tmfile) || die "ERROR: Failed to open file: $tmfile";
              while(<W>) {
                if($GDT) {
                  $gdt = $1 if /^GDT-TS-score=\s+(\S+)\s+/;
                  last if /^\(":" denotes the residue pairs/;
                } else {
                  $tm1 = $1 if /^TM\-score=\s+(\S+)\s+\(if normalized by length of Chain_1/;
                  $tm2 = $1 if /^TM\-score=\s+(\S+)\s+\(if normalized by length of Chain_2/;
                  last if /^\(":" denotes residue pairs/;
                }
              }
              my $tmqaln = <W>; my $tmraln = <W>; $tmraln = <W>; chomp($tmqaln); chomp($tmraln);
              close(W);
              my $strprn;
              if($GDT) {
                $strprn = sprintf("%6s %10s %6s %9s gtm1= %8s gtm2= %8s gbest= %8.6f  gdtts= %8s\n",
                  $QID,$rfn,$QSTR.$QCHN,$rstr.$rchn,$gtm1,$gtm2,($gtm1<$gtm2)?$gtm2:$gtm1, $gdt);
              } else {
                $strprn = sprintf("%6s %10s %6s %9s gtm1= %8s gtm2= %8s gbest= %8.6f  tm1= %8s tm2= %8s best= %8.6f\n",
                  $QID,$rfn,$QSTR.$QCHN,$rstr.$rchn,$gtm1,$gtm2,($gtm1<$gtm2)?$gtm2:$gtm1, $tm1,$tm2,($tm1<$tm2)?$tm2:$tm1);
              }
              syswrite(STDOUT, $strprn);
              if((length($tmqaln) == length($qaln) && length($tmraln) == length($raln)) ||
                 (substr($tmqaln,$qnnb+$rnnb,$qalnlen) eq substr($qaln,$qnnb+$rnnb,$qalnlen) && substr($tmraln,$qnnb+$rnnb,$qalnlen) eq substr($raln,$qnnb+$rnnb,$qalnlen)))
              {
                unlink($wfile, $tmfile, $trgfile);
              } else {
                print(STDERR "WARNING: Inconsistent alignments: $QID $rfn: $QSTR $QCHN  $rstr $rchn\n");
              }
            }
          } #if(GetTrgFile)
        } #if(ne)
      } #if($rfnfile)
      $gtm1 = $gtm2 = -1; $qaln = $raln = $rfn = ""; $tm1 = $tm2 = $gdt = -1;
      $qlen = $rlen = $qend = $rend = -1; undef($qstart); undef($rstart);
    } #if(data)
  }
  next unless $hits;
  do {$rfn = $1; $rfn =~ s/^\S+\/(\S+)\.cif\.gz\s+Chn:(\S+).*$/${1}_$2/;} if /"reference_description":\s+"([^"]+)"/;
  $qlen = $1 if /"query_length":\s+(\d+)/;
  $rlen = $1 if /"reference_length":\s+(\d+)/;
  $gtm2 = $1 if /"tmscore_refn":\s+([\d\.]+)/;
  $gtm1 = $1 if /"tmscore_query":\s+([\d\.]+)/;
  $qstart = $1 if /"query_from":\s+(\d+)/;
  $qend = $1 if /"query_to":\s+(\d+)/;
  $rstart = $1 if /"refn_from":\s+(\d+)/;
  $rend = $1 if /"refn_to":\s+(\d+)/;
  $qaln = $1 if /"query_aln":\s+"([^"]+)"/;
  $raln = $1 if /"refrn_aln":\s+"([^"]+)"/;
}

close(F);

## -------------------------------------------------------------------

sub GetTrgFile
{
  my ($rfnfile, $trgfile, $rchn) = @_;
  my $tmpfile = "${trgfile}_full.cif";
  my $unpfile = "${trgfile}_chain";
  my $cmd = "gunzip -c $rfnfile >$tmpfile";
  if(system($cmd) != 0) {
      print(STDERR "ERROR: Failed to execute: $cmd\n");
      return 0;
  }
  $cmd = "python3 $getchainprog -i $tmpfile -c $rchn -o $unpfile";
  if(system($cmd) != 0) {
      print(STDERR "ERROR: Failed to execute: $cmd\n");
      return 0;
  }
  unless(open(R, $unpfile)) {
      print(STDERR "ERROR: Failed to open file: $unpfile\n");
      return 0;
  }
  unless(open(W,">$trgfile")) {
      print(STDERR "ERROR: Failed to open file for writing: $trgfile\n");
      close(R);
      return 0;
  }
  my ($prvchn, $prvstr, $prvnum) = ("_","    ",0);
  my ($rec, $n, $ca, $c, $o) = ("",0,0,0,0);
  while(<R>) {
    if(/^(?:ATOM|TER|HETATM)/) {
      if((substr($_,22,4) ne $prvstr) || (substr($_,21,1) ne $prvchn)) {
        $prvchn = substr($_,21,1);
        $prvstr = substr($_,22,4);
        $prvnum++;
        print(W $rec) if($ca); ##if($n && $ca && $c && $o);
        $rec = ""; $n = $ca = $c = $o = 0;
      }
      $ca = 1 if(substr($_,12,4) eq " CA " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A"));
      $n = 1  if(substr($_,12,4) eq " N  " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A"));
      $c = 1  if(substr($_,12,4) eq " C  " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A"));
      $o = 1  if(substr($_,12,4) eq " O  " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A"));
      $rec .= $_;
    }
    else {
      print(W $rec) if($ca && $rec); ##if($n && $ca && $c && $o && $rec);
      $rec = ""; $n = $ca = $c = $o = 0;
      print(W);
    }
    last if /^ENDMDL/;
  }
  print(W $rec) if($ca && $rec); ##if($n && $ca && $c && $o && $rec);
  close(W);
  unlink($tmpfile, $unpfile);
  return 1;
}

