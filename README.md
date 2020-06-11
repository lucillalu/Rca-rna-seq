# Rca-rna-seq

1.测试环境
服务器CPU为Intel Xeon E4280，主频为2.0GHZ，内存为64GB，操作系统版本为Linux Ubuntu 16.04。

2.操作与使用
（1）分别在src/abPOA-master、src/deSALT-master、src/deBGA-master目录下执行make指令；在src/目录下执行make指令编译生成可执行文件 rrs
（2）索引参考基因组命令如下：
./rrs index genome_file genome_index_dir
分割RNA序列命令如下：
./rrs par genome_index_dir fasta_file
比对RNA序列命令如下：
./rrs aln genome_index_dir hq_reads.fasta
得到结果文件aln.sam
（3）评测结果文件，具体方法见evaluation/目录下readme.md
（4）模拟实验数据，具体方法见simulation/目录下readme.md