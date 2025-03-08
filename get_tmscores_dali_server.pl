#!/usr/bin/env perl
BEGIN {$^W=1}

use strict;

## NOTE: all directories and executables must exist!
my $WRKDIR = "/home/mindaugas/projects/data/gtalign-web-benchmark";
my $TMalign ="/home/mindaugas/install/TM-align/TMalign";
my $TMscore = "/home/mindaugas/install/TM-score--2022-2-27/TMscore_from_alignment";
my $getchainprog = "/data/installed-software/comer-ws-backend/bin/getchain.py";

my $QRYDIR = "$WRKDIR/queries";
my $RFNDIR = "/data-SSD/databases/wwpdb/PDB";
my $RDIR = "$WRKDIR/Dali-server";
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

my $resfile = "$RDIR/$QID.html";
my $cnt = 0;
my ($nrtt, $qidtt, $rfntt, $zsctt);
my ($nr, $rfn, $zsc) = ("","",-1);
my ($qseq, $rseq) = ("","");
my ($tm1, $tm2, $gdt) = (-1,-1,-1);

open(F, $resfile) || die "ERROR: Failed to open results: $resfile";

while(<F>) {
  chomp;
  if(eof(F) || /\s+No\s+(\S+)\s+Query=(\S+)\s+Sbjct=(\S+)\s+Z\-score=([\d\.]+)/)
  {
    $nrtt = ""; $qidtt = $rfntt = ""; $zsctt = -1;
    do { $nrtt = $1; $qidtt = $2; $rfntt = $3; $zsctt = $4; } unless eof(F);
    if($rfn && $qseq && $rseq)
    {
      my $wfile = "$WDIR/${QID}_${rfn}.aln";
      my $tmfile = "$WDIR/${QID}_${rfn}.tm";
      my $rstr = substr($rfn,0,4);
      my $rchn = substr($rfn,4,1);
      my $rmdl = substr($rfn,1,2);
      my $rfnfile = "$RFNDIR/$rmdl/pdb${rstr}.ent.gz";
      unless(-f $rfnfile) {
        print(STDERR "ERROR: Reference $rfn not found: $rfnfile\n");
      }
      else {
        last if $MAXNHITS <= $cnt++;
        if("$QSTR$QCHN" ne "$rstr$rchn")
        {
          my $trgfile = "$WDIR/${QID}_${rfn}.pdb";
          if(GetTrgFile($rfnfile, $trgfile, $rchn))
          {
            open(W, ">$wfile") || die "ERROR: Failed to open file for writing: $wfile";
            print(W ">$QSTR $QCHN $QID\n$qseq\n>$rstr $rchn $rfn\n$rseq\n");
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
                $strprn = sprintf("No.%-5s %6s %6s %6s %6s Z= %5s  gdtts= %8s\n",
                  $nr,$QID,$rfn,$QSTR.$QCHN,$rstr.$rchn,$zsc,$gdt);
              } else {
                $strprn = sprintf("No.%-5s %6s %6s %6s %6s Z= %5s  tm1= %8s tm2= %8s best= %8.6f\n",
                  $nr,$QID,$rfn,$QSTR.$QCHN,$rstr.$rchn,$zsc,$tm1,$tm2,($tm1<$tm2)?$tm2:$tm1);
              }
              syswrite(STDOUT, $strprn);
              if(length($tmqaln) == length($qseq) && length($tmraln) == length($rseq)) {
                unlink($wfile, $tmfile, $trgfile);
              } else {
                print(STDERR "WARNING: Inconsistent alignments: $QID $rfn: $QSTR $QCHN  $rstr $rchn\n");
              }
            }
          } #if(GetTrgFile)
        } #if(ne)
      } #if($rfnfile)
    } #if(data)
    $nr = $nrtt; $rfn = $rfntt; $zsc = $zsctt; $qseq = $rseq = ""; $tm1 = $tm2 = $gdt = -1;
  }
  $qseq .= $1 if /^Query\s+(\S+)\s+/;
  $rseq .= $1 if /^Sbjct\s+(\S+)\s+/;
}

close(F);

## -------------------------------------------------------------------

sub GetTrgFile
{
  my ($rfnfile, $trgfile, $rchn) = @_;
  my $tmpfile = "${trgfile}_full";
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
        print(W $rec) if($n && $ca && $c && $o);
        $rec = ""; $n = $ca = $c = $o = 0;
      }
      if(substr($_,12,4) eq " CA " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A")) {
        $ca = 1;
      }
      $n = 1  if(substr($_,12,4) eq " N  " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A"));
      $c = 1  if(substr($_,12,4) eq " C  " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A"));
      $o = 1  if(substr($_,12,4) eq " O  " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A"));
      $rec .= $_;
    }
    else {
      print(W $rec) if($n && $ca && $c && $o && $rec);
      $rec = ""; $n = $ca = $c = $o = 0;
      print(W);
    }
    last if /^ENDMDL/;
  }
  print(W $rec) if($n && $ca && $c && $o && $rec);
  close(W);
  unlink($tmpfile, $unpfile);
  return 1;
}

