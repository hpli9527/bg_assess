# bg_assess
用于对回交后代的背景恢复率进行评估
# 使用说明
```shell
# 首先从变异信息文件中找到CHROM、POS、REF、ALT、QUAL、对应的亲本和回交后代个体基因型信息，并过滤存在缺失基因型、亲本和后代间没有多态性的位点
cut -f1,2,4,5,6,42,13,73 Bnabg.filter.SNPs.txt | grep -vP '^scaffold\d' | grep -v '\./\.' | grep -v '0/0.*0/0.*0/0' | grep -v '1/1.*1/1.*1/1' | grep -v '0/1.*0/1.*0/1' > QS_Q10.txt
# 滑窗计算每个窗口和供体亲本基因型相同的SNP，并初步绘图
Rscript bg.R -v QS_Q01.txt -r s4270 -d QS_T01 -s QS_Q01
# 需要人工判断导入片段的边界并创建导入片段位置文件
Rscript.exe ./plot_bg.R -c ./chromosome_length.txt -m ./S8722_2.bg_estimate.txt -r ./S8722_2.region.txt -o ./S8722_2.bg_plot
# 也可以不提供导入片段位置文件，将不会画出对应信息
Rscript.exe ./plot_bg.R -c ./chromosome_length.txt -m ./S8722_2.bg_estimate.txt -o ./S8722_2.bg_plot
```
## 文件示例
### QS_Q10.txt
|CHROM|POS|REF|ALT|QUAL|s8800|S195G_1|Y648|
|---|---|---|---|---|---|---|---|
|scaffoldA01|276015|C|T|2777.51|0/1:11,32:43:99:1120,0,366|1/1:0,6:6:18:227,18,0|0/1:2,18:20:12:613,0,12|
|scaffoldA01|326049|A|T|5037.56|1/1:0,46:46:99:1\|1:326049_A_T:2003,138,0|1/1:0,8:8:27:1\|1:326049_A_T:402,27,0|0/1:1,35:36:99:0\|1:326049_A_T:1400,0,121|
|scaffoldA01|326064|C|T|4830.56|1/1:0,41:41:99:1\|1:326049_A_T:1845,123,0|1/1:0,8:8:24:1\|1:326049_A_T:360,24,0|0/1:1,30:31:99:0\|1:326049_A_T:1344,0,127|
|...|...|...|...|...|...|...|...|
### chromosome_length.txt
|CHROM|LEN|
|---|---|
|A01|38004428|
|A02|35943954|
|A03|44868710|
|...|...|
### S8722_2.bg_estimate.txt
|chr|start|end|sum|mid|
|---|---|---|---|---|
|A01|410735|1252351|0|831543|
|A01|705375|1252545|0|978960|
|A01|1047073|1252548|0|1149810.5|
|...|...|...|...|...|
### S8722_2.region.txt, 对于导入片段的边界需要人工判断
|CHROM|START|END|COLOR|
|---|---|---|---|
|A10|0|15033101|#66ff33|
|C08|39867737.5|42366716.5|#66ff33|
|C09|20792390|26689228|#66ff33|
|...|...|...|...|
## 结果示例
<img src="https://user-images.githubusercontent.com/35584208/125244766-491d9200-e322-11eb-8e3f-02318beddd5a.png" width = "400" alt="S8722_2" />
<img src="https://user-images.githubusercontent.com/35584208/125244786-4fac0980-e322-11eb-9b5d-9f1c007bef00.png" width = "400" alt="S8722_2 bg_plot1" />
<img src="https://user-images.githubusercontent.com/35584208/125244801-53d82700-e322-11eb-91eb-8a1ba9051495.png" width = "400" alt="S8722_2 bg_plot2" />
