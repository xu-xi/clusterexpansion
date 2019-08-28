# ATAT使用简介
## 下载与安装
1. 从[ATAT主页](https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/)下载压缩包,稳定版或者测试版均可。
2. 解压缩
    ```sh
    $ gunzip atatX_XX.tar.gz
    $ tar -xvf atatX_XX.tar.gz
    ```
3. 修改`Makefile`中的安装路径 
4. 编译
    ```sh
    $ make
    $ make install
    ```
5. 添加`$PATH`
6. 运行`ezvasp`生成配置文件`~/.ezvasp.rc`并编辑以修改vasp运行的命令和势函数所在位置

## 必要输入文件
### lat.in
格点文件，具体格式请参见[ATAT手册](https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/manual.pdf)
```sh
$ makelat Al,Ti fcc
$ makelat Ba:In:O,Vac E21
```
上面两个例子分别产生Al-Ti合金fcc的格点文件和含有O缺陷的BaInO<sub>1-x</sub>的钙钛矿结构的格点文件。

说明：
+ 不同的位点用冒号分隔开，同一位点混占的原子之间用逗号分隔开，空位用`Vac`表示
+ 结构文件的模板在`atatX_XX/data`中,可以在`atatX_XX/data/str`中找到更多并复制到上一级目录使用
+ 初始的晶胞参数根据原子半径产生，原子半径的默认值在`atatX_XX/data/radii.in`文件中
+ 如果所需结构不在结构库中，可以使用`py_conv`转换得到，`py_conv`为[tmckit](https://github.com/ccme-tmc/tmckit)的脚本，其使用请参考其帮助文档

### vasp.wrap
计算参数控制文件，支持INCAR的绝大多数参数

说明：
+ `KPPRA` 为不同大小的晶胞设置统一的k点密度，若`KPPRA=1000`，则当晶胞中只有一个原子时，k点数目为1000，有两个原子时，k点数目为500，即原子数*k点数=1000。目前还不支持通过`KSPACING`设置k点密度
+ `USEPOT`指定所用赝势或PAW势，与`~/.ezvasp.rc`中定义势文件的位置的变量相对应
+ `SUBATOM` 可用于使用vasp的[推荐赝势](https://cms.mpi.univie.ac.at/vasp/vasp/Recommended_PAW_potentials_DFT_calculations_using_vasp_5_2.html)，例如`SUBATOM=s/Ti$/Ti_sv/g`即表示对Ti原子用Ti_sv的势
+ `DOSTATIC` 在做完结构弛豫后再做一个静态计算，这样最后的能量更为准确


## 可选输入文件
### crange.in
例子：
+ 直接限制组分
    ```sh
    1*B<=0.2
    ```
+ 电荷平衡
    ```sh
    -2*O+4*Ce+3*Sm<=0.2
    -2*O+4*Ce+3*Sm>=-0.2
    ```
说明：
+ 只能用不等式限制而不能用等式，并且ATAT依然会产生少数在限制条件之外的结构来检测可能的基态构型
+ 如果条件限制得过于严格，在产生了一定数量的结构后ATAT产生结构的速度将会极慢，原因不明。这种情况建议自己写脚本产生所需的结构。
## 运行

```sh
maps -d &
pollmach runstruct_vasp
```
说明：
+ `maps`不断产生结构和构建团簇展开模型，`-d`表示使用默认参数，直接输入`maps`可查看可用参数
+ `runstruct_vasp`调用`vasp`计算所产生的结构的能量，每个目录包含一个结构及其计算结果，每个目录下的`energy`文件即从vasp的计算中读取的能量
+ `maps.log`中记载了现有的团簇展开模型的状态，包括对基态构型的预测和LOOCV等

什么时候停止计算：
+ 团簇展开模型能够成功预测基态 (查看maps.log)
+ LOOCV足够小 (<0.025 eV suggested by ATAT manual) 
+ 有效团簇相互作用随着团簇距离的增大而衰减

以上的三个判据并不是绝对的，判断一个团簇展开模型是合理的还可以有其它方式。


停止计算：
```
$ touch stop
$ touch stoppoll
```

重启计算: 与正常开始一个计算相同

#### 常用命令
+ `corrdump -2= -3= -4=` 产生团簇
+ `getclus` 查看团簇
    + 第一列为团簇的阶数，第二列为团簇直径，第三列为团簇的多重度
    + `clusters.out`中有团簇格点位置的详细信息
    
+ `clusterexpand -e -cv energy` 得到团簇展开系数
    + 此时产生的ECI文件为`energy.eci`，而`maps`产生的ECI为`eci.out`，后者是形成能的团簇展开结果，前者是直接对vasp能量(即`energy`文件)的团簇展开
    + 如果需要对其它性质做团簇展开，例如带隙，则可以在每个结构的目录下创建一个`bandgap`文件,其内容即为带隙的大小，然后运行`clusterexpand -e -cv -pa bandgap`以得到带隙团簇展开的系数`bandgap.eci`，其中`-pa`表示所展开的量已经是一个平均的量

#### 检查计算结果
+ VASP是否正常结束
    + ATAT会自动检查一些错误(`checkerr_vasp`)
    + 但如果结构优化的离子步数不够大，不会有提示
+ 晶胞的应变是否很大
    + `checkrelax` (<0.1 suggested by ATAT）

#### 结果可视化
+ `mapsrep` 会调用`gnuplot`作图并产生相应的脚本`mapsrep.gnu`


#### 热力学计算
必要的输入文件：
+ `lat.in`
+ `clusters.out`
+ `eci.out`
+ `gs_str.out`

#####  `emc2`使用示例
蒙特卡洛模拟
```sh
emc2 -T0=300 -T1=3000 -dT=100 -cm -x=0 -keV -gs=-1 -er=30 -aq=0 -dx=1e-5 -sigdig=12
```
+ `-cm`表示正则系综
+ `-x`其实是点团簇函数，若c表示格点的占据浓度，则x=2c-1
+ `-sigdig=12`是为了表示有足够的位数输出能量的方差以判断相变
+ `-aq`和`-dx`是自动确定MC平衡步数和取样步数的算法，其原理可参考文献[Ref.1](#ref_1)，也可以用`-eq`和`-n`来直接指定步数，需要做收敛性测试
+ 主要输出文件为`mc.out`，输入`emc2 -h`可查看每一列所表示的内容
+ 输出文件`mcsnapshot.out`为最后一个模拟温度下的MC快照。如果想得到某个温度下的snapshot，可以只设置T0而不设置T1和dT。
    ```sh
    $ cp mcsnapshot.out str.out
    $ runstruct_vasp -nr 
    $ vesta POSCAR
    ```
#### 相图计算
基本原理可参考文献[Ref.1](#ref_1)，Martin Baker的[笔记](https://arxiv.org/abs/1907.10151)提供了更多的例子和细节

例子：
```sh
phb -dT=10 -ltep=1e-3 -er=30 -gs1=0 -gs2=1 -o=phb.out -keV -dx=1e-5
```

#### 构建SQS/SQoS
SQS/SQoS是对无序体系的代表性结构，其中SQS假定原子完全随机的占据，而SQoS是SQS的扩展，可以考虑一定的短程有序性，其基本原理可参考相关文献[Ref.2](#ref_2)，Eric的[网页](http://grandcentral.apam.columbia.edu:5555/tutorials/dft_procedures/sqs/index.html)提供了一个例子和流程

必要文件：
+ rndstr.in
+ clusters.out

`rndstr.in`的一个例子：

```
4.1000 0.000000 0.000000
0.000000 4.1000 0.000000
0.000000 0.000000 4.1000
1.000000 0.000000 0.000000
0.000000 1.000000 0.000000
0.000000 0.000000 1.000000
0.000000 0.500000 0.500000 O=0.6667,N=0.3333
0.500000 0.000000 0.500000 O=0.6667,N=0.3333
0.500000 0.500000 0.000000 O=0.6667,N=0.3333
0.000000 0.000000 0.000000 Ba
0.500000 0.500000 0.500000 Ta
```

运行`mcsqs -n= &`
说明：
+ `-n`设定SQS晶胞的原子数
+ 可通过`-tcf`设定目标团簇函数（可从`mc.out`中提取），即SQoS方法
+ 停止计算：`touch stopsqs`
+ 输出文件为`mcsqs.log`和`bestsqs.out`
+ 将`bestsqs.out`复制为`str.out`，然后运行`runstruct_vasp -nr`可生成vasp的输入文件

## 参考文献
1. <a name="ref_1"></a> A. van de Walle and M. Asta, Self-driven lattice-model Monte Carlo simulations of alloy thermodynamic properties and phase diagrams, [Modell. Simul. Mater. Sci. Eng., 2002, 10, 521](https://iopscience.iop.org/article/10.1088/0965-0393/10/5/304) 
2. <a name="ref_2"></a> J. Liu, M. V. Fernandez-Serra and P. B. Allen, Special quasiordered structures: Role of short-range order in the semiconductor alloy (GaN)<sub>1−x</sub>(ZnO)<sub>x</sub>, [Phys. Rev. B, 2016, 93, 054207](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.93.054207)
