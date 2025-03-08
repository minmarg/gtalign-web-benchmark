#!/bin/bash
name=$(basename $0 .sh)
## PARAMS
#SORT=legend
MAXN=1000
## FILES
WIDTH=9.5; HEIGHT=3.2
WIDTH=7.08; HEIGHT=2.4
if [ "$SORT" == "legend" ]; then WIDTH=2.4; fi

GTAW_speed9_pre04="gtalign_web_speed9_prescore04-processed.gdtts.out"
GTAW_speed13_pre04="gtalign_web_speed13_prescore04-processed.gdtts.out"
DALI_TM="Dali-server-processed.out"
DALI_GDT="Dali-server-processed.gdtts.out"
FS3d_TM="Foldseek-server-3diaa-processed.out"
FS3d_GDT="Foldseek-server-3diaa-processed.gdtts.out"
FStm_TM="Foldseek-server-tmalign-processed.out"
FStm_GDT="Foldseek-server-tmalign-processed.gdtts.out"

## FUNCTIONS
function GetXYSum() {
  local dir=$1
  local cnt=$2
  local fieldqid=$3
  local fieldrfn=$4
  local fieldv1=$5
  local fieldv2=$6
  local fieldk=$7
  local srt=$8
  NUM=$cnt  FLDQID=$fieldqid FLDRFN=$fieldrfn   FLDV1=$fieldv1 FLDV2=$fieldv2 FLDK=$fieldk SRT=$srt TYPE=$SORT perl -e '
    while(<>){
      @a=split(/\s+/); push @sk,$a[$ENV{FLDK}];
      $id=$a[$ENV{FLDQID}].$a[$ENV{FLDRFN}];
      push @ids,$id;
      push @sv, ($a[$ENV{FLDV1}]<$a[$ENV{FLDV2}])? $a[$ENV{FLDV2}]: $a[$ENV{FLDV1}];
      last if ($ENV{TYPE} eq "legend") && (10000 < $c++);
    } 
    if($ENV{SRT} eq "A") {
      @sndx=sort{$sv[$a]<=>$sv[$b]} 0..$#sv; 
    } else {
      @sndx=sort{$sv[$b]<=>$sv[$a]} 0..$#sv; 
    }
    $y=$x="c("; 
    $sum=0;
    $NUM=$ENV{NUM};
    for($i=$k=0;$k<$NUM && $i<=$#sndx;$i++){
      ##NOTE: GTalign-web includes different models: use one (always the first by chain extraction):
      next if exists $h{ $ids[$sndx[$i]] };
      $h{ $ids[$sndx[$i]] } = 1;
      $last=($k+1)*0.001 if 0.5<=$sk[$sndx[$i]]; 
      $sum+=$sk[$sndx[$i]]; 
      if($k%400==0 || $NUM<=$k+1 || $#sndx<$i+1){
        do{$y.=",";$x.=","} if $k; 
        $y.=sprintf("%.4f",$sum*0.001); 
        $x.=sprintf("%.3f",($k+1)*0.001);
      }
      $k++;
    } 
    $y.=")"; $x.=")"; 
    print $x, " ", $y, " ", $last, " ", $sum*0.001' $dir/*.out
}

## OUTPUT
output="${name}_${MAXN}${SORT}"

## total number of Dali hits:
NN=65371



## sort by tm-score
GTAW_speed13_pre04_TM_dat=("GTalign-web speed=13" $(GetXYSum "../$GTAW_speed13_pre04" $NN  0 1  5 5 5 D) 34)
GTAW_speed9_pre04_TM_dat=("GTalign-web speed=9" $(GetXYSum "../$GTAW_speed9_pre04" $NN  0 1  5 5 5 D) 48)
DALI_TM_dat=("DALI server" $(GetXYSum "../$DALI_TM" $NN  1 2  8 8 8 D) 477.4)
FS3d_TM_dat=("Foldseek server mode=3Di/AA" $(GetXYSum "../$FS3d_TM" $NN  0 1  11 11 11 D) 11.7)
FStm_TM_dat=("Foldseek server mode=TM-align" $(GetXYSum "../$FStm_TM" $NN  0 1  11 11 11 D) 14.7)


## sort by gdt_ts
GTAW_speed13_pre04_GDT_dat=("GTalign-web speed=13" $(GetXYSum "../$GTAW_speed13_pre04" $NN  0 1  11 11 11 D))
GTAW_speed9_pre04_GDT_dat=("GTalign-web speed=9" $(GetXYSum "../$GTAW_speed9_pre04" $NN  0 1  11 11 11 D))
DALI_GDT_dat=("DALI server" $(GetXYSum "../$DALI_GDT" $NN  1 2  8 8 8 D))
FS3d_GDT_dat=("Foldseek server mode=3Di/AA" $(GetXYSum "../$FS3d_GDT" $NN  0 1  11 11 11 D))
FStm_GDT_dat=("Foldseek server mode=TM-align" $(GetXYSum "../$FStm_GDT" $NN  0 1  11 11 11 D))



Rtext="
library(ggplot2);
library(scales);
library(gridExtra);

mypalette <- c('#006000', '#466f46',  '#ff7000',  '#BC6FF1', '#892CDC')

myltypes <- c('224282F2', 'F282',  'solid',  '224282F2', 'F1')

dat <- list(
  '${GTAW_speed13_pre04_TM_dat[0]}'=cbind(x=${GTAW_speed13_pre04_TM_dat[1]},y=${GTAW_speed13_pre04_TM_dat[2]}), 
  '${GTAW_speed9_pre04_TM_dat[0]}'=cbind(x=${GTAW_speed9_pre04_TM_dat[1]},y=${GTAW_speed9_pre04_TM_dat[2]}), 
  '${DALI_TM_dat[0]}'=cbind(x=${DALI_TM_dat[1]},y=${DALI_TM_dat[2]}), 
  '${FS3d_TM_dat[0]}'=cbind(x=${FS3d_TM_dat[1]},y=${FS3d_TM_dat[2]}), 
  '${FStm_TM_dat[0]}'=cbind(x=${FStm_TM_dat[1]},y=${FStm_TM_dat[2]}));

list.names <- names(dat)
lns <- sapply(dat, nrow)
dat <- as.data.frame(do.call('rbind', dat))
dat\$group <- rep(list.names, lns)

p10wlegend <- ggplot(dat, aes(x = x, y = y, colour = group)) +
  theme_bw() +
  theme(legend.text = element_text(size = 8), legend.key.width = unit(0.04,'npc'), legend.key.height = unit(0.02,'npc')) +
  geom_line(aes(linetype = group), linewidth = 0.3) +
  scale_linetype_manual(name = '', breaks = list.names, values = myltypes) +
  scale_colour_manual(name = '', breaks = list.names, values = mypalette)

p12wlegend <- ggplot(dat, aes(x = x, y = y)) + ##, colour = group)) +
  theme_bw() +
  theme(legend.text = element_text(size = 8), legend.key.width = unit(0.02,'npc'), legend.key.height = unit(0.02,'npc')) +
  geom_point(aes(shape = group), size = 1.8) +
  scale_shape_manual(name = '', breaks = list.names, values = seq(0,length(list.names))) +
  scale_colour_manual(name = '', breaks = list.names, values = mypalette)

vlines <- c(
  ${GTAW_speed13_pre04_TM_dat[3]},
  ${GTAW_speed9_pre04_TM_dat[3]},
  ${DALI_TM_dat[3]},
  ${FS3d_TM_dat[3]},
  ${FStm_TM_dat[3]});

P1 <- ggplot(dat, aes(x = x, y = y, colour = group, linetype = group)) +
  ggtitle('query length-normalized TM-score') +
  theme_bw() +
  theme(legend.position = 'none', plot.title = element_text(size = 8), axis.title = element_text(size = 8), axis.text = element_text(size = 8)) +
  labs(x=expression(paste('# top hits ('%*%10^{3},')')), y=expression(paste('Cumulative TM-score ('%*%10^{3},')'))) +
  geom_line(linewidth = 0.3) +
  scale_linetype_manual(breaks = list.names, values = myltypes) +
  scale_colour_manual(breaks = list.names, values = mypalette) +
  geom_vline(xintercept = vlines, linetype = myltypes, color = mypalette, linewidth = 0.2)


dat <- list(
  '${GTAW_speed13_pre04_TM_dat[0]}'=cbind(x=${GTAW_speed13_pre04_TM_dat[5]},y=${GTAW_speed13_pre04_TM_dat[4]}), 
  '${GTAW_speed9_pre04_TM_dat[0]}'=cbind(x=${GTAW_speed9_pre04_TM_dat[5]},y=${GTAW_speed9_pre04_TM_dat[4]}), 
  '${DALI_TM_dat[0]}'=cbind(x=${DALI_TM_dat[5]},y=${DALI_TM_dat[4]}), 
  '${FS3d_TM_dat[0]}'=cbind(x=${FS3d_TM_dat[5]},y=${FS3d_TM_dat[4]}), 
  '${FStm_TM_dat[0]}'=cbind(x=${FStm_TM_dat[5]},y=${FStm_TM_dat[4]}));

list.names <- names(dat)
dat <- as.data.frame(do.call('rbind', dat))
dat\$group <- rep(list.names, 1)

P2 <- ggplot(dat, aes(x = x, y = y, colour = group)) +
  ggtitle('') +
  theme_bw() +
  theme(legend.position = 'none', panel.grid.minor = element_blank(),
    plot.title = element_text(size = 8), axis.title = element_text(size = 8), axis.text = element_text(size = 8)) +
  labs(x=expression(paste('Runtime (min)')), y='') + ##expression(paste('Cumulative TM-score ('%*%10^{3},')'))) +
  ylim(0,NA) +
  geom_point(aes(shape = group), size = 1.4) + #, stroke = 1) +
  annotation_logticks(sides = 'b', short = unit(0.003,'npc'), mid = unit(0.006,'npc'), long = unit(0.01,'npc')) +
  scale_shape_manual(breaks = list.names, values = seq(0,length(list.names))) +
  scale_colour_manual(breaks = list.names, values = mypalette) +
  scale_x_continuous(trans = log10_trans(),
    limits = c(5,1000),
    breaks = c(10,100,1000), #trans_breaks('log10', function(x) 10^x),
    labels = trans_format('log10', math_format(10^.x)))

sprintf('TM-scores:');
sprintf('%-40s %.4f %.4f','${GTAW_speed13_pre04_TM_dat[0]}',${GTAW_speed13_pre04_TM_dat[4]},${GTAW_speed13_pre04_TM_dat[5]});
sprintf('%-40s %.4f %.4f','${GTAW_speed9_pre04_TM_dat[0]}',${GTAW_speed9_pre04_TM_dat[4]},${GTAW_speed9_pre04_TM_dat[5]});
sprintf('%-40s %.4f %.4f','${DALI_TM_dat[0]}',${DALI_TM_dat[4]},${DALI_TM_dat[5]});
sprintf('%-40s %.4f %.4f','${FS3d_TM_dat[0]}',${FS3d_TM_dat[4]},${FS3d_TM_dat[5]});
sprintf('%-40s %.4f %.4f','${FStm_TM_dat[0]}',${FStm_TM_dat[4]},${FStm_TM_dat[5]});
sprintf('');


dat <- list(
  '${GTAW_speed13_pre04_GDT_dat[0]}'=cbind(x=${GTAW_speed13_pre04_GDT_dat[1]},y=${GTAW_speed13_pre04_GDT_dat[2]}), 
  '${GTAW_speed9_pre04_GDT_dat[0]}'=cbind(x=${GTAW_speed9_pre04_GDT_dat[1]},y=${GTAW_speed9_pre04_GDT_dat[2]}), 
  '${DALI_GDT_dat[0]}'=cbind(x=${DALI_GDT_dat[1]},y=${DALI_GDT_dat[2]}), 
  '${FS3d_GDT_dat[0]}'=cbind(x=${FS3d_GDT_dat[1]},y=${FS3d_GDT_dat[2]}), 
  '${FStm_GDT_dat[0]}'=cbind(x=${FStm_GDT_dat[1]},y=${FStm_GDT_dat[2]}));

list.names <- names(dat)
lns <- sapply(dat, nrow)
dat <- as.data.frame(do.call('rbind', dat))
dat\$group <- rep(list.names, lns)

P3 <- ggplot(dat, aes(x = x, y = y, colour = group, linetype = group)) +
  ggtitle('alignment length-normalized GDT_TS') +
  theme_bw() +
  theme(legend.position = 'none',
    plot.title = element_text(hjust=1,size = 8),
    plot.title.position = 'plot', axis.title = element_text(size = 8), axis.text = element_text(size = 8)) +
  labs(x=expression(paste('# top hits ('%*%10^{3},')')), y=expression(paste('Cumulative GDT_TS ('%*%10^{3},')'))) +
  ##geom_line(linetype = 1, linewidth = 0.4) +
  geom_line(linewidth = 0.3) +
  scale_linetype_manual(breaks = list.names, values = myltypes) +
  scale_colour_manual(breaks = list.names, values = mypalette)
  ##scale_colour_discrete(breaks = list.names)

sprintf('GDT_TS scores:');
sprintf('%-40s %.4f','${GTAW_speed13_pre04_GDT_dat[0]}',${GTAW_speed13_pre04_GDT_dat[4]});
sprintf('%-40s %.4f','${GTAW_speed9_pre04_GDT_dat[0]}',${GTAW_speed9_pre04_GDT_dat[4]});
sprintf('%-40s %.4f','${DALI_GDT_dat[0]}',${DALI_GDT_dat[4]});
sprintf('%-40s %.4f','${FS3d_GDT_dat[0]}',${FS3d_GDT_dat[4]});
sprintf('%-40s %.4f','${FStm_GDT_dat[0]}',${FStm_GDT_dat[4]});
sprintf('');


# extract legend from ggplot (credit: Joachim Schork)
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1\$grobs, function(x) x\$name) == 'guide-box')
  step3 <- step1\$grobs[[step2]]
  return(step3)
}

legend10 <- extract_legend(p10wlegend);
legend12 <- extract_legend(p12wlegend);

#plot <- grid.arrange(arrangeGrob(P1, P2, P3, ncol=3), common_legend, nrow=2, heights = c(10, 2.4));
if('${SORT}'=='legend') {
  plot <- grid.arrange(arrangeGrob(legend10, legend12, nrow=2), ncol=1);
} else {
  plot <- grid.arrange(arrangeGrob(P1, P2, P3, ncol=3), nrow=1);
}

ggsave(filename='${output}.pdf',plot,width=${WIDTH},height=${HEIGHT})
#dev.off();
"
echo "$Rtext" | R --vanilla --slave

