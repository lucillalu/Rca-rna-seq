# RCA-RNA-seq

模拟器: PBSIM.
基因组: Drosophila melanogaster DM6 (fruitfly).
基因注释: DM6(version 94) from Ensemble.
测试版本：Python3.7/Ubuntu 16.04

使用方法：
python Eval_sim_data.py simulation_folder alignment.sam annotations.gtf group_list All_SS_iso.txt out.csv hq_read_file

 [1] simulation_folder: PBSIM生成的模拟数据集的文件夹相对路径
 
 [2] alignment.sam: 结果文件相对路径
 
 [3] annotations.gtf: 基因注释文件相对路径
 
 [4] group_list: PBSIM产生的不同测序深度的数据集的group的相对路径
 
 [5] All_SS_iso.txt: 含所有单剪接转录本的标识符文件（simulation时生成）
 
 [6] out.csv: 结果文件
 
 [7] hq_read_file：比对过程中产生的高质量RNA序列文件


