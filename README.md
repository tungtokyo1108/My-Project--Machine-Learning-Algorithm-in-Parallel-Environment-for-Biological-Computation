Journal
------------------------------------------------------------------------------------------------------------------------------------------
Title: Stochastic Variational Inference for Bayesian Phylogenetics: A Case of CAT Model

Authors: Tung Dang and Hirohisa Kishino

Link: https://academic.oup.com/mbe/article/36/4/825/5305894 

Video: https://www.youtube.com/watch?v=_3u06pL-zds

Our presentation (included guides of software): https://drive.google.com/file/d/16bKeL7sQQtv82Tw17euQ-y6lZP9MlHzx/view 

概要
------------------------------------------------------------------------------------------------------------------------------------------
がんは、細胞が無制限に増殖する病気であり、DNAへの変異の蓄積がおもな要因です。多くのがん種では、細胞分裂のたびに、ゲノム中にさまざまな塩基変異を蓄積します。この、がんゲノムの変化の過程は、生物種分化の際のゲノム進化と似ており、機械学習など最先端の人工知能技術を用いて、がんの進行過程を調べることができます。私は、機械学習で使われる隠れマルコフモデル (Hidden Markov Model, HMM) や混合モデル (Mixture Model)を用いて、がんゲノムのシーケンシングデータから、がんの増殖過程を推定する手法の開発を進めています。HMM の状態数や混合モデルの混合数は本来、観測対象および観測データに適切に選択されなければならない。カテゴリが少なすぎればDNAの変異の粗い記述しかできず，多すぎれば極端な場合，一変数一カテゴリとなって，カテゴリ化が無意味となってしまう。私の方法は、全ての未知量を確率変数と見なし、それらを観測データを得た下での事後分布としてベイズ定理を用いて推定するノンパラメトリックベイズ法（ディリクレ過程）です。実際、ディリクレ過程では、データが観測される毎に、要素分布数が必要に応じて増える柔軟なデータ生成過程となっています。ディリクレ過程に基づいて、混合数が予め固定されているわけではなく観測データに応じて定めます。

しかしながら、ノンパラメトリックベイズモデルを学習する際、実用上、(1)局所最適性 (local optimality)、および、(2) モデルの複雑さ (model complexity) の決定、の問題に悩まされる。(1) はマルコフ連鎖モンテカルロ(MCMC)アルゴリズムが所望の大域的最適解に収束せず、初期解の近傍の局所最適解に収束するという問題である。また (2) は、混合モデルで は混合要素数の決定など、非線形モデルの構造決定問題を指す。モデルの複雑さを学習タスクの複雑さに応じて適切に定めないと汎化能力(未学習データに対する予測能力)の低下を招く。私は、変分ベイズ学習(variational Bayesian learning)に導入し、非線形モデルの上記 (1),(2)の問題を同時解決する新たな学習法(最良モデル探索のための変分ベイズ学習法)を提案し、混合モデルへの適用実験により手法の有効性を示す。変分ベイズでは複雑すぎて計算することできないような事後分布を、もっと簡単な分布（変分分布）に近似するとします。しかし、変分分布がどんな分布でも良いわけではなく、可能な限り事後分布に近づけることを選択します。確率分布を、変分分布と事後分布がどれだけ近い確率分布であるかは、KLダイバージェンス式で計算できます。つまり、上記のKLダイバージェンスの最小化問題を解くことで、事後分布に近所した変分分布を計算することが可能となります。


Abstract 
-----------------------------------------------------------------------------------------------------------------------------------------
The pattern of molecular evolution varies among gene sites and genes in a genome. By taking into account the complex heterogeneity of evolutionary processes among sites in a genome, Bayesian infinite mixturemodels of genomic evolution enable robust phylogenetic inference. With large modern data sets, however,the computational burden of Markov chain Monte Carlo sampling techniques becomes prohibitive. Here, we have developed a variational Bayesian procedure to speed up the widely used PhyloBayes MPI program, which deals with the heterogeneity of amino acid propensity. Rather than sampling from the posterior distribution, the procedure approximates the (unknown) posterior distribution using a manageable distribution called the variational distribution. The parameters in the variational distribution are estimated by minimizing Kullback-Leibler divergence. To examine performance, we analyzed three large data sets consisting of mitochondrial, plastid-encoded, and nuclear proteins. Our variational method accurately approximated the Bayesian phylogenetic tree, mixture proportions, and the amino acid propensity of each component of the mixture while using orders of magnitude less computational time.

-----------------------------------------------------------------------------------------------------------------------------------------

This package implements a scalable, parallel implementation of the fast algorithm for fitting 
a model of biological computation on large data set. 

Note: All files have name like "....VI.h/cpp" which we have updated new optimization algorithm for computating the parameters of model
