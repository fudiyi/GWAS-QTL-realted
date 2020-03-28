# perl脚本汇总

### 这是一些以前师兄用过的脚本，部分整理如下，若有时间可查阅相关使用方法，目前推荐编写python脚本完成

#### add_P_into_gene.pl

```perl
#!usr/bin/perl
open GENE, 'GWAS_p_results.txt' or die "cannot open gene file.";
while(<GENE>){
	chomp;	
	@gene=split/\t+/;
	if($gene[10]=~/(\AA.*FG)(\d\d\d)/){
		print "$1P$2\t$_\n";}
	elsif($gene[10]=~/(\AG.*\d)/){print"$1\t$_\n";}
	elsif($gene[10]=~/(NA)/){print "$1\t$_\n";}
	}
```

#### CMDscirpt.pl

```perl
#!/usr/bin/perl
#用来生成在windows平台上用R批量作图的CMD文件。生成的.cmd文件要配合R脚本一起使用。
#每次作图时，两个地方需要改：
#1. n的初始数值和while的循环次数要和脚本的名字一致
#2.脚本每次放置的路径需要改
$n=1;
open CMD, '+> zuotu.cmd' or die;
while($n<19){
print CMD "R CMD BATCH \"E:\\daixie\\190304_warm\\Rscript\\$n.R\"\n";
$n++;
}
```

#### find_SNP.pl

```perl
#!usr/bin/perl
open HAPMAP, '/data/GWAS/Yanglab_GWAS/513_genotype.hmp' or die "cannot open hapmap";
my %hash;
while(<HAPMAP>){
	chomp;
	@m=split/\t+/;
	$hash{$m[0]}=$m[1];
}
close HAPMAP;

open FILE, 'all.txt' or die "cannot opne all.txt";
while(<FILE>){
	chomp($_);
	@n=split/\t+/;
	if($hash{$n[2]}){
	print "$hash{$n[2]}\t$_";
	}
	else{
	print"NA\t$_";
	}
}
close FILE;

```

#### find_the_smallest.pl

```perl
#!usr/bin/perl
$filename=1;
while($filename<=1407){
$l=0;
$num=0;
$tmp=0;
@line=();
open FILE, "./results/$filename.txt" or die "cannot open file";
$newname=$filename+1407;
print "$newname.warm.txt\t";
@file=<FILE>;
foreach $x(@file){
        chomp($x);
        @line=split/\s+/,$x;
	if($line[6]>0){
		$P=log($line[6])/log(10)*-1;
		if($P>$tmp){
		$tmp=$P;
		$l=$x;
		}
	$num++;
	}
}
print"$l\n";
$filename++;
close FILE;
}

```

#### get_gene_description.pl

```perl
#!usr/bin/perl
my %hash;
open DES, '/data/software/maize_gene/maize_5b_FGS_anno.tab' or die "cannot open the description file.";
while(<DES>){
	chomp;
	@des=split/\t+/;
	if($des[0]=~/\AA/){
	$hash{$des[0]}=$des[1];
	}
	elsif($des[0]=~/\AG/){
	@aa=split/_/,$des[0];
	$hash{$aa[0]}=$des[1];
	}
}
close DES;
open GENE, 'add_P' or die "cannot open gene file.";
while(<GENE>){
	chomp;
	@gene=split/\t+/;
	if($hash{$gene[0]}){
	print "$hash{$gene[0]}\t$_\n";
	}
	else{
	print "NA\t$_\n";
	}

}
close GENE;

```

#### rename.pl

```perl
#!/usr/bin/perl
#一般放入GWAS结果所在的目录使用，将有空的x_output2.txt统一改名为x.txt方便作图时R脚本运行
foreach my$file (glob "*2.txt"){ #查找当前目录下所有以"2.txt"结尾的文件
my $newFile = $file;
$newFile =~ s/\A(\d*)_output2.txt/$1.txt/; #取出文件名中的数字并替换为纯数字。
if(-e $newFile){ #如果修改后会导致文件重名，则输出警告，不作处理
warn "Can't rename $file to $newFile. The $newFile exists!\n";
}else{
rename $file, $newFile #重命名文件
or
warn "Rename $file to $newFile failed: $!\n"; #如果重命名失败，则输出警告
}
}

```

#### Rscript.pl

```perl
#!usr/bin/perl
#用来生成windows平台上用R作图时的脚本文件。配合.cmd文件可以实现一键作n张图。
#每次作图有两处需要修改：
#前提是GWAS结果文件时以x.txt来命名的（x为数字）
#例如结果文件为1.txt到100.txt共100个，则需要修改的地方是：
#1.while()循环中改为n<101，使n循环100次，输出100个脚本
#2.GWAS结果文件所存放的路径每次都需要修改（setwd后面跟的路径）
$n=1;
while($n<352){
open S, "+> $n.R" or die;
print S "library(ggplot2)
library(grid)
setwd(\"/data/FDY_analysis/mGWAS/analysis/GWAS/results_cold\")
fl=\"$n.txt\"
data=read.table(file=fl,sep=\"\\t\",head =T,stringsAsFactors=F,fill=T)
data\$p=as.numeric(data\$p)
data=data[-1,]
data=data[!data\$Chr==0,]
chr=aggregate(data\$Pos,by=list(data\$Chr),FUN=\"max\")
names(chr)=c(\"Chr\",\"pos\")

p0=ggplot(data=data, aes(Pos/1000000, -log10(p)))+
geom_point(aes(colour = factor(Chr)),size=1)+
facet_grid(.~ Chr, space=\"free_x\",scales=\"free_x\")+
geom_abline(intercept =5, slope = 0,linetype=2, size=0.3)+
labs(title='$n')+
geom_hline(yintercept=-0.03,linetype=1,size=0.2)+
geom_text(data=chr,aes(x=pos/2000000,y=0,label=Chr),size=3,vjust=1.4)+
labs(x='Chromosome', y=expression(-log[10]~~(P[observed])))

p=p0+theme_bw() +theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    panel.grid=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    panel.border=element_rect(size=0.3),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.spacing=unit(0,\"lines\"),
    axis.line=element_line(size=0.3,linetype=0),
    axis.line.y=element_line(size=0.2,linetype=1),
    legend.position=\"none\")+
  scale_color_manual( values=rep(c(\"orangered\",\"darkblue\"),times=5))+
  scale_y_continuous(limits=c(-0.5,15))
ggsave(filename=\"$n.perl.tiff\",p,width=120,height=100,units=\"mm\")

";
$n++;
}

```

