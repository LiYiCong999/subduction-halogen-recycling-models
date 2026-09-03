# Nature Communications投稿前的GitHub–Zenodo代码发布指南

本指南适用于本项目的Fig. 4、Fig. 5和Fig. 6三个Python脚本。建议把三个脚本放在同一个GitHub仓库中，并为论文最终使用的代码版本创建一个GitHub Release，再由Zenodo自动存档并生成DOI。

推荐仓库名称：

    subduction-halogen-recycling-models

推荐仓库描述：

    Python codes and input data for Monte Carlo modelling of halogen release, slab-derived fluids, mélange formation and mantle mixing.

## 第一阶段：上传前整理

### 第1步：确定最终代码版本

1. 保留一份原始代码备份。
2. 确认提交到论文中的Fig. 4、Fig. 5和Fig. 6分别由仓库中的哪一个脚本生成。
3. 不要在创建GitHub Release以后继续直接修改同一个版本。
4. 如果代码仍会修改，先使用开发版本，例如v0.9.0；论文投稿对应的稳定版本再命名为v1.0.0。

本模板使用以下英文文件名：

    fig4_aoc_release_contributions.py
    fig5_scenario1_fluid_mantle_mixing.py
    fig6_scenario2_melange_mantle_mixing.py

三个仓库副本仅修改了文件名，代码内容与上传的三个脚本一致。

### 第2步：整理输入数据

至少应提供以下内容：

- Fig. 4所用AOC、蓝片岩、HP榴辉岩和UHP榴辉岩卤素数据；
- Fig. 5所用全岩卤素和H2O数据、矿物模式、分配系数和阶段定义；
- Fig. 6所用蛇纹岩、变沉积岩、蓝片岩、HP榴辉岩和UHP榴辉岩数据；
- 每条数据的样品编号、岩性、单位和原始来源；
- 文献汇编数据的完整参考文献；
- 数据筛选、单位换算和排除规则。

如果数据可以公开，建议直接放入仓库的data目录。如果部分文献数据不允许重新分发，应在data/README.md中列出原始来源、DOI和获取方式，并提供能够合法共享的最小输入数据。

不要上传：

- 未经许可的第三方版权表格；
- 审稿人身份或审稿意见；
- 尚未公开的合作者数据且未取得同意；
- 密码、访问令牌、邮箱授权码或GitHub Personal Access Token；
- 电脑用户名和私人绝对路径。

### 第3步：核对三个脚本的路径

目前三个脚本都依赖相对路径。推荐在code目录中运行。

Fig. 4要求输入文件位于code目录：

    AOC-BS-HP-UHP Halogens data.csv

Fig. 5要求在code目录下建立：

    data source/

并放入四个CSV文件。Fig. 6要求输入文件位于code目录：

    Forearc serpentinite-Metasediment-BS-HP-UHP halogen_data.csv

运行前必须保证文件名、空格、大小写和扩展名完全一致。

### 第4步：记录软件环境

在真正生成论文最终图件的电脑上运行：

    python --version
    python -m pip show numpy pandas matplotlib scipy

把版本写入README，并生成锁定文件：

Windows PowerShell：

    python -m pip freeze | Select-String "numpy|pandas|matplotlib|scipy" | Set-Content requirements-lock.txt

Windows命令提示符：

    python -m pip freeze | findstr /I "numpy pandas matplotlib scipy" > requirements-lock.txt

macOS或Linux：

    python -m pip freeze | grep -Ei "numpy|pandas|matplotlib|scipy" > requirements-lock.txt

必须使用实际生成论文图件的环境，而不是随意填写最新版本。

### 第5步：在干净环境中复现

建议在一个新的文件夹中测试：

    python -m venv .venv

激活环境并安装依赖后，从code目录依次运行三个脚本。核对：

- 程序没有报错；
- 输出CSV包含有效结果；
- Fig. 4、Fig. 5和Fig. 6与论文中的最终图一致；
- 固定随机种子的脚本重复运行产生相同结果；
- 图中文字、单位、颜色、坐标范围和误差棒定义正确。

Fig. 4没有固定随机种子，因此不能保证逐次运行完全一致。投稿前最好决定是否增加固定种子，并在代码、README和Methods中保持一致。

### 第6步：准备许可证

公开仓库并不自动允许他人合法复用代码。需要与所有代码作者确认许可证。科研脚本常用MIT、BSD-3-Clause或GPL-3.0。

本模板包含MIT许可证文本模板。如果选择MIT：

1. 将LICENSE_TEMPLATE.txt中的年份和版权持有人替换完整；
2. 将文件重命名为LICENSE；
3. 确认所有代码所有者同意。

不要同时保留多个互相冲突的许可证。

### 第7步：准备引用元数据

填写CITATION.cff.template：

- title：代码仓库正式标题；
- authors：所有软件作者，顺序需经作者确认；
- ORCID：使用完整ORCID链接；
- version：首次正式版本建议1.0.0；
- date-released：实际发布日期；
- repository-code：GitHub仓库网址；
- url：Zenodo DOI生成后补入；
- license：与LICENSE一致。

填写完成后将其重命名为CITATION.cff，并放在仓库根目录。

Zenodo会优先读取.zenodo.json；如果仓库同时存在.zenodo.json，则CITATION.cff不会控制Zenodo元数据。因此，除非你明确需要.zenodo.json，否则不要添加它。

## 第二阶段：建立GitHub仓库

### 第8步：登录或注册GitHub

1. 打开https://github.com。
2. 登录个人账号。
3. 建议在GitHub个人资料中添加ORCID链接。
4. 确认账号邮箱已经验证。

### 第9步：新建空仓库

1. 点击右上角加号。
2. 选择New repository。
3. Owner选择你的个人账号或课题组组织。
4. Repository name填写：

       subduction-halogen-recycling-models

5. Description填写前文建议的英文描述。
6. Visibility选择Public。
7. 如果准备上传本模板中的完整目录，不要在网页中预先勾选Add a README、Add .gitignore或Choose a license，以免和本地文件冲突。
8. 点击Create repository。

如果论文采用双盲审稿，先阅读本指南的“双盲审稿注意事项”，不要直接公开带有作者身份的仓库。

## 第三阶段：上传文件

### 方法A：使用GitHub网页上传

该方法适合文件较少、每个文件小于25 MiB的情况。

1. 进入新建的空仓库。
2. 点击uploading an existing file，或者Add file → Upload files。
3. 把本模板根目录中的文件和文件夹拖入上传区域。
4. 检查目录层级。README.md必须位于仓库根目录，而不是多套一层文件夹。
5. Commit message填写：

       Initial release of halogen modelling codes

6. 选择Commit directly to the main branch。
7. 点击Commit changes。
8. 回到仓库首页，确认README已正常显示。

网页上传单个文件不能超过25 MiB；GitHub会阻止超过100 MiB的普通Git文件。大型Monte Carlo结果不应直接塞入代码仓库，可压缩后作为Zenodo数据记录或单独的数据集存档。

### 方法B：使用Git命令行上传

该方法更适合需要持续修改和版本管理的情况。

在本地打开Git Bash、PowerShell或终端，进入仓库模板目录：

    cd path/to/nc-halogen-models

初始化并提交：

    git init
    git add .
    git status
    git commit -m "Initial release of halogen modelling codes"
    git branch -M main

连接远程仓库：

    git remote add origin https://github.com/YOUR-USERNAME/subduction-halogen-recycling-models.git

上传：

    git push -u origin main

如果GitHub要求认证，应使用浏览器登录、GitHub CLI、SSH密钥或Personal Access Token。不要把Token写入脚本、README或任何被提交的文件。

### 第10步：检查GitHub仓库

逐项打开并检查：

- README.md是否正常显示；
- 三个代码文件是否可以打开；
- docs目录是否包含三个脚本的注意事项；
- 数据文件或数据访问说明是否完整；
- requirements-lock.txt是否记录实际版本；
- CITATION.cff是否没有REPLACE、TODO或占位符；
- LICENSE是否已经确认；
- 所有私人路径、用户名和敏感信息是否已经移除；
- 仓库中是否不存在大体积临时文件和缓存文件。

建议在GitHub搜索框中分别搜索：

    D:/
    C:/
    Desktop
    TODO
    REPLACE
    password
    token

确认没有不应公开的内容。

## 第四阶段：连接Zenodo

### 第11步：登录Zenodo

1. 打开https://zenodo.org。
2. 使用GitHub账号登录，或登录Zenodo后在账号设置中连接GitHub。
3. 建议同时关联ORCID。

### 第12步：授权GitHub集成

1. 在Zenodo个人设置中进入GitHub或Applications相关页面。
2. 点击连接或授权GitHub。
3. GitHub会显示Zenodo请求的权限，确认后授权。
4. 返回Zenodo的GitHub仓库列表。
5. 如果新仓库没有出现，点击Sync now或刷新仓库列表。
6. 找到subduction-halogen-recycling-models。
7. 把该仓库右侧开关切换为On或Enabled。

必须在创建GitHub Release之前启用该仓库，否则Zenodo可能不会自动接收该Release。

如果仓库属于GitHub Organization但未显示，需要让组织管理员批准Zenodo应用访问。

## 第五阶段：创建GitHub正式版本

### 第13步：发布前最后检查

在创建Release前再次确认：

- 论文最终使用的代码已经提交到main分支；
- 输入数据和README对应；
- 三个脚本均能从干净环境运行；
- 版本依赖已经记录；
- CITATION.cff已经填写；
- LICENSE已经确认；
- Release中不包含占位符；
- 所有作者同意公开。

### 第14步：创建v1.0.0 Release

1. 进入GitHub仓库首页。
2. 点击右侧Releases。
3. 点击Draft a new release。
4. 点击Choose a tag。
5. 输入v1.0.0。
6. 选择Create new tag: v1.0.0 on publish。
7. Target选择main。
8. Release title填写：

       v1.0.0 – Nature Communications submission release

9. Release notes建议填写：

       This release contains the Python codes used to generate Figs. 4–6 of the associated manuscript. It includes the stage-specific halogen-release model, the slab-fluid–mantle mixing scenario, the mélange–mantle mixing scenario, documentation, and reproducibility information.

10. 不要勾选Set as a pre-release。
11. 点击Publish release。

GitHub Release不是普通Commit。Zenodo监听的是Release事件，因此只有提交代码但不创建Release不会自动获得Zenodo DOI。

## 第六阶段：检查Zenodo记录并获得DOI

### 第15步：等待Zenodo存档

发布GitHub Release后：

1. 返回Zenodo的GitHub集成页面。
2. 等待数分钟。
3. 打开Uploads或My dashboard。
4. 查找与v1.0.0对应的软件记录。

如果没有生成：

- 确认仓库开关在Release发布前已经打开；
- 确认发布的是正式Release而不是Draft；
- 确认Zenodo仍有GitHub访问权限；
- 刷新仓库列表；
- 检查GitHub Release页面是否正常。

不要连续重复发布多个相同标签。

### 第16步：核对Zenodo元数据

打开Zenodo记录并检查：

- Resource type：Software；
- Title：代码正式标题；
- Creators：软件作者姓名和ORCID；
- Description：代码与论文的关系及三个模型的用途；
- Version：v1.0.0或1.0.0；
- Publication date：Release日期；
- License：与GitHub LICENSE一致；
- Keywords：halogens、subduction、Monte Carlo simulation、mantle mixing、geochemistry；
- Related identifiers：GitHub仓库网址、论文预印本DOI或文章DOI；
- Funding：与代码开发直接相关的项目；
- Files：GitHub Release归档文件是否存在。

Zenodo发布后的文件和持久标识符不能直接替换，但元数据可以修改。如果代码文件需要变化，应创建新版本，而不是试图覆盖v1.0.0。

### 第17步：选择论文中引用的DOI

Zenodo通常为具体版本和整个版本集合提供不同标识。论文需要严格复现v1.0.0时，应引用该版本对应的DOI，而不是只引用始终指向最新版本的概念DOI。

建议把版本号和版本DOI同时写入：

- 论文Code availability；
- README；
- 参考文献；
- CITATION.cff；
- 给审稿人的说明。

### 第18步：把DOI徽章加入README

Zenodo记录页面通常提供DOI badge代码。复制Markdown形式并放在README标题下方。例如：

    [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

替换XXXXXXX后提交：

    git add README.md CITATION.cff
    git commit -m "Add Zenodo DOI"
    git push

注意：v1.0.0的Zenodo归档不会因为后续普通Commit而改变。如果希望DOI也存在于归档文件内部，需要创建后续Release，例如v1.0.1。也可以在首次Release前使用手动Zenodo上传并预留DOI，但GitHub–Zenodo自动集成本身不能预留DOI。

## 第七阶段：写入Nature Communications稿件

### Code availability

放置位置：

    Data availability
    Code availability
    References

推荐正文：

    Code availability

    The custom Python code used for the bootstrap Monte Carlo calculations, modelling of slab-derived fluid compositions, construction of mélange endmembers, mantle-mixing calculations and visualization of halogen-ratio distributions is archived on Zenodo at https://doi.org/10.5281/zenodo.XXXXXXX. The archived release (v1.0.0) contains the source code, software documentation and information required to reproduce Figs. 4–6. The corresponding development repository is available at https://github.com/USERNAME/subduction-halogen-recycling-models.

如果输入数据在同一个Zenodo记录中：

    The input datasets required to reproduce the analyses are included in the archived release.

如果数据另有DOI：

    The input and source datasets are available at https://doi.org/XXXX.

### 参考文献

建议把软件作为独立参考文献引用：

    Author Surname, Initials. et al. Halogen recycling models for subduction-zone metamorphism. Zenodo, version 1.0.0, https://doi.org/10.5281/zenodo.XXXXXXX (2026).

作者顺序必须与Zenodo和CITATION.cff一致。

### Figure legends

图例至少需要说明：

- Monte Carlo迭代次数；
- 固定随机种子；
- 数据点、箱线图、KDE和直方图分别代表什么；
- KDE等级是最大密度的比例还是累积概率；
- 中央95%如何筛选；
- 误差棒是±1σ、置信区间还是其他统计量；
- n代表原始样品数量还是模拟次数。

## 第八阶段：返修和新版本

如果审稿期间修改代码：

1. 在GitHub修改代码并提交；
2. 更新CHANGELOG或Release notes；
3. 根据修改程度创建v1.0.1或v1.1.0；
4. 发布新的GitHub Release；
5. Zenodo自动生成新的版本记录和版本DOI；
6. 在返修稿Code availability中引用真正生成最终结果的新版本DOI；
7. 不要继续引用已经不对应最终图件的v1.0.0 DOI。

建议版本规则：

- v1.0.1：注释、README或不改变结果的小修正；
- v1.1.0：增加分析或改变部分非核心输出；
- v2.0.0：改变核心算法、模型假设或输出结果。

## 双盲审稿注意事项

公开GitHub仓库、作者姓名、ORCID和Zenodo记录会直接暴露作者身份。如果选择Nature Communications的double-anonymized peer review：

1. 不要在匿名稿件中直接加入带作者姓名的公开仓库链接；
2. 按期刊匿名审稿指南准备代码访问方式；
3. 可向编辑提供私密审稿链接或使用能够生成匿名审稿链接的代码存档服务；
4. 在公开Zenodo记录前与编辑确认；
5. 接收后再公开带完整作者元数据的正式版本。

如果采用常规单盲审稿，则可在投稿前公开GitHub和Zenodo并直接引用。

## 最终检查清单

- [ ] 三个脚本与论文最终图件一一对应
- [ ] 输入数据文件齐全或具有永久访问地址
- [ ] 数据单位和列名已经核对
- [ ] Fig. 4随机种子问题已经处理或披露
- [ ] Fig. 5水含量单位已经确认
- [ ] Fig. 5缺失分配系数设为0的规则已经披露
- [ ] Fig. 6固定50:5:45比例具有依据
- [ ] DMM卤素组成具有文献来源
- [ ] requirements-lock.txt来自实际运行环境
- [ ] README没有占位符
- [ ] CITATION.cff没有占位符
- [ ] LICENSE已获所有作者同意
- [ ] 没有密码、Token、私人路径或未经授权数据
- [ ] 已在干净环境中复现全部结果
- [ ] Zenodo仓库开关已在Release前启用
- [ ] GitHub v1.0.0 Release已发布
- [ ] Zenodo版本DOI已经核对
- [ ] Code availability引用正确版本DOI
- [ ] 软件参考文献已经加入参考文献表
- [ ] 图例定义KDE等级、中央95%和±1σ
