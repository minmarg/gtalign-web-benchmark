export WRKDIR=/home/mindaugas/projects/data

## get singleton structure chains
mkdir single.org
awk '{print $1}' ${WRKDIR}/gtalign-benchmark-afspdb/singletons.lst | xargs -i sh -c 'n=$(echo {}|sed -re "s/_.+\$//"); c=$(echo {}|sed -re "s/^.+_//"); gunzip -c ../wwpdb-CRISPR-Cas/$n.pdb.gz >$n.pdb; python3 /data/installed-software/comer-ws-backend/bin/getchain.py -i $n.pdb -c $c -o single.org/${n}_${c}.pdb; rm $n.pdb'

## run gtalign for singletons
bash -c '(time /home/mindaugas/local/gtalign/bin/gtalign -v --qrs=${WRKDIR}/gtalign-web-benchmark/single.org --rfs=/data-SSD/databases/wwpdb/mmCIF -o gtalign_16_speed9_prescore04_addss_s05__single --hetatm --dev-queries-per-chunk=2 --dev-queries-total-length-per-chunk=1500 --dev-min-length=3 --dev-max-length=4000 --speed=9 --pre-score=0.4 -s 0.5 --add-search-by-ss --nhits=20000 --nalns=20000 --referenced --ter=0 --split=2 --sort=2 --dev-N=,2 -c cachedir) >gtalign_16_speed9_prescore04_addss_s05__single.log 2>&1 &'
rm cachedir/*

## get lengths
l -S gtalign_16_speed9_prescore04_addss_s05__single|perl -e 'open(F,"/home/mindaugas/projects/data/gtalign-benchmark-afspdb/singletons.lst")or die"Error: singletons.lst";while(<F>){chomp;@a=split(/\s+/);$h{$a[0]}=$a[1]}close(F); while(<>){chomp;@a=split(/\s+/);$o=$a[7];$o=$1 if $o=~/^(.+)__.+$/;$l="";$l=$h{$o} if exists $h{$o}; printf("%-80s %10s\n",$_,$l)}' >gtalign_16_speed9_prescore04_addss_s05__single.lst

## randomly select 100 different length representatives
awk '{if($9>=100)print}' gtalign_16_speed9_prescore04_addss_s05__single.lst|tac|perl -e 'while(<>){@a=split(/\s+/);$l=$a[8]; do{$l1++;print if $l1<=10} if($l>=100 && $l<200); do{$l2++;print if $l2<=30} if($l>=200 && $l<300); do{$l3++;print if $l3<=30} if($l>=300 && $l<400); do{$l4++;print if $l4<=30} if($l>=400);}' >gtalign_16_speed9_prescore04_addss_s05__single100.lst

## copy the structures
cut -b60-66 gtalign_16_speed9_prescore04_addss_s05__single100.lst|xargs -i cp single.org/{}.pdb queries.org/

## prepare structures so that all methods interpret them equivalently: no residues w/o " CA [ |A]" atoms;
## leave only residues that have at least N, CA, C, and O atoms (valid for Dali);
## also, remove HETATM residues which makes foldseek and TM-align read the same sequence of residues
mkdir queries 
ls -1 queries.org | xargs -i -P 40 perl -e '
  open(F,">queries/{}") or die "ERROR: Failed to open file for writing.";
  $prvchn = "_"; $prvstr="    "; $prvnum=0;
  while(<>) {
    if(/^(?:ATOM|TER)/) {
      if((substr($_,22,4) ne $prvstr) || (substr($_,21,1) ne $prvchn)) {
        $prvchn = substr($_,21,1);
        $prvstr = substr($_,22,4);
        $prvnum++;
        print(F $rec) if($n && $ca && $c && $o);
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
      print(F $rec) if($n && $ca && $c && $o && $rec);
      $rec = ""; $n = $ca = $c = $o = 0;
      print(F) unless /^HETATM/;
    }
    last if /^ENDMDL/;
  }
  print(F $rec) if($n && $ca && $c && $o && $rec);
  close(F)' queries.org/{}

## -------------------------------------------------------------------------

## run Dali server
## NOTE: 2h36_X, 3t6a_C, 3t6a_D, 6bl5_A, and 7wwv_A were tested before run (Dsrv doesn't run repeatedly);
## NOTE: their average length is 206 (108, 308, 308, 129, 179);
## NOTE: estimate cumulative time by the time of 5swc_C (208): 87 (sec) * 5 = 435, and this to the total time
ls -1 queries | perl -e '
  use POSIX;
  $odir = "Dali-server";
  $http = "http://ekhidna2.biocenter.helsinki.fi/cgi-bin/sans/dump.cgi";
  while(<>) {
    chomp;
    $qname = $_; $qname =~ s/.pdb$//; $Q{$qname}{CH} = $1 if $qname =~ /_(\S)$/;
    $submit = `curl --form file1=\@queries/$qname.pdb --form method=search --form submit=Submit "$http"`;
    unless($submit =~ /^Redirected\s+to\s+(http:\S+)/m) {
      $dt = strftime "%Y-%m-%d %H:%M:%S", gmtime time;
      print(STDERR "$dt E: Query failed: $qname: [$submit]\n");
      sleep(1);
      next;
    }
    $Q{$qname}{HP} = $1;
    $dt = strftime "%Y-%m-%d %H:%M:%S", gmtime time;
    print("$dt I: $qname submitted: $Q{$qname}{HP}\n");
    sleep(1);
  }
  $N = scalar(keys %Q);
  $dt = strftime "%Y-%m-%d %H:%M:%S", gmtime time;
  print("$dt I: Queries submitted: $N\n");
  while(1) {
    $n = 0;
    foreach $q(keys %Q) {$n++ if exists $Q{$q}{TE}}
    last if $N <= $n;
    foreach $q(keys %Q) {
      $addr = $Q{$q}{HP};
      $ch = $Q{$q}{CH};
      next if(exists $Q{$q}{TE});
      if(exists $Q{$q}{TB}) {
        $res = `curl --head --silent --fail $addr/s001${ch}.html`;
        if($res =~ /Content-Type:\s+text\/html/m) {
          $Q{$q}{TE} = time();
          $dt = strftime "%Y-%m-%d %H:%M:%S", gmtime time;
          print("$dt I: $q finished: $Q{$q}{TE}\n");
          $get = `wget -O $odir/$q.html $addr/s001${ch}.html`;
        }
        sleep(1);
        next;
      }
      $log = `wget -O $odir/$q.log $addr`;
      unless(open(F, "$odir/$q.log")) {
        $dt = strftime "%Y-%m-%d %H:%M:%S", gmtime time;
        print(STDERR "$dt E: Status failed: $q: [$log]\n");
        sleep(1);
        next;
      }
      while(<F>) {
        if(/Status:\s+Running/){
          $Q{$q}{TB} = time();
          $dt = strftime "%Y-%m-%d %H:%M:%S", gmtime time;
          print("$dt I: $q running: $Q{$q}{TB}\n");
          last
      } }
      close(F);
      sleep(1);
    }
  }
' >dali_server.out 2>dali_server.err &

## calculate runtimes
## NOTE: only for inspection because the total runtime is 52489 sec, considerably 
## NOTE: larger than the runtime obtained after all jobs have finished (28643 sec)
## NOTE: due to parallel processing
perl -e 'while(<>){$h{$1}=$2 if /I:\s+(\S+)\s+running:\s+(\d+)/; if(/I:\s+(\S+)\s+finished:\s+(\d+)/){$t+=$2-$h{$1};printf("%10s %10d\n",$1,$2-$h{$1})}} printf("Total: %10d\n",$t)' dali_server.out >dali_server.runtimes

## calculate TM-scores
ls -1 queries | xargs -i -P 10 sh -c 'n=$(basename {} .pdb); ./get_tmscores_dali_server.pl {} >Dali-server-processed.out/${n}.out 2>Dali-server-processed.err/${n}.err'

## calculate GDT_TS scores
ls -1 queries | xargs -i -P 20 sh -c 'n=$(basename {} .pdb); ./get_tmscores_dali_server.pl {} gdt >Dali-server-processed.gdtts.out/${n}.out 2>Dali-server-processed.gdtts.err/${n}.err'

## -------------------------------------------------------------------------

## run Foldseek server

export MODE=3diaa
export MODE=tmalign
ls -1 queries | perl -e '
  use POSIX;
  $mode = "$ENV{MODE}";
  $odir = "Foldseek-server-$mode";
  $http = "https://search.foldseek.com/api/ticket";
  while(<>) {
    chomp;
    $qname = $_; $qname =~ s/.pdb$//;
    next if -f "$odir/${qname}.json";
    $submit = `curl -X POST -F q=\@queries/$qname.pdb -F "mode=$mode" -F "database[]=pdb100" "$http"`;
    unless($submit =~ /{"id":"([^"]+)"/m) {
      $dt = strftime "%Y-%m-%d %H:%M:%S", gmtime time;
      print(STDERR "$dt E: Query failed: $qname: [$submit]\n");
      last;
    }
    $Q{$qname}{HP} = $1;
    $Q{$qname}{CH} = $1 if $qname =~ /_(\S)$/;
    $dt = strftime "%Y-%m-%d %H:%M:%S", gmtime time;
    print("$dt I: $qname submitted: $Q{$qname}{HP}\n");
  }
  sleep(1);
  $N = scalar(keys %Q);
  $dt = strftime "%Y-%m-%d %H:%M:%S", gmtime time;
  print("$dt I: Queries submitted: $N\n");
  while(1) {
    $n = 0;
    foreach $q(keys %Q) {$n++ if exists $Q{$q}{TE}}
    last if $N <= $n;
    foreach $q(keys %Q) {
      $addr = $Q{$q}{HP};
      $ch = $Q{$q}{CH};
      next if(exists $Q{$q}{TE});
      $res = `curl -X GET https://search.foldseek.com/api/ticket/$addr`;
      if($res =~ /"status":"COMPLETE"/m) {
        $Q{$q}{TE} = time();
        $dt = strftime "%Y-%m-%d %H:%M:%S", gmtime time;
        print("$dt I: $q finished: $Q{$q}{TE}\n");
        $get = `curl -X GET https://search.foldseek.com/api/result/$addr/0 >$odir/$q.json`;
        next;
      }
      if($res =~ /"status":"RUNNING"/m) {
        $Q{$q}{TB} = time();
        $dt = strftime "%Y-%m-%d %H:%M:%S", gmtime time;
        print("$dt I: $q running: $Q{$q}{TB}\n");
      }
    }
    sleep(1);
  }
' >>foldseek_server_${MODE}.out 2>>foldseek_server_${MODE}.err &

## calculate TM-scores
ls -1 queries | xargs -i -P 20 sh -c 'n=$(basename {} .pdb); ./get_tmscores_foldseek_server.pl {} >Foldseek-server-3diaa-processed.out/${n}.out 2>Foldseek-server-3diaa-processed.err/${n}.err'

ls -1 queries | xargs -i -P 20 sh -c 'n=$(basename {} .pdb); ./get_tmscores_foldseek_server.pl {} >Foldseek-server-tmalign-processed.out/${n}.out 2>Foldseek-server-tmalign-processed.err/${n}.err'

## calculate GDT_TS scores
ls -1 queries | xargs -i -P 20 sh -c 'n=$(basename {} .pdb); ./get_tmscores_foldseek_server.pl {} gdt >Foldseek-server-3diaa-processed.gdtts.out/${n}.out 2>Foldseek-server-3diaa-processed.gdtts.err/${n}.err'

ls -1 queries | xargs -i -P 20 sh -c 'n=$(basename {} .pdb); ./get_tmscores_foldseek_server.pl {} gdt >Foldseek-server-tmalign-processed.gdtts.out/${n}.out 2>Foldseek-server-tmalign-processed.gdtts.err/${n}.err'

## -------------------------------------------------------------------------

## process GTalign-web results

## calculate GDT_TS scores
## (GTalign TM-scores can be slightly lower than those produced by TM-align given a GTalign alignment due to relaxed superposition refinement)
ls -1 queries | xargs -i -P 20 sh -c 'n=$(basename {} .pdb); ./get_tmscores_gtalign_web.pl {} gdt >gtalign_web_speed13_prescore04-processed.gdtts.out/${n}.out 2>gtalign_web_speed13_prescore04-processed.gdtts.err/${n}.err'

