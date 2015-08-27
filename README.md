# STOD

This is a mirror from the repository located on Illimine
<a href="http://illimine.cs.uiuc.edu/software/">here</a>. This only
contains the matlab source files; it does not contain the data. The link
on Illimine contains the data.

**** Illimine Copyright ****

University of Illinois at Urbana-Champaign, 2015

illimine.cs.illinois.edu


**********************************************
Additional Copyright
**********************************************

This package contains the source code and the dataset used in the following paper:

Chi Wang, Xueqing Liu, Yanglei Song, Jiawei Han. Scalable Moment-based
Inference for Latent Dirichlet Allocation. ECMLPKDD, 2014.

If you use any contents in this package, please cite:

@inproceedings{wang14,
  title={Scalable Moment-based Inference for Latent Dirichlet Allocation},
  author={Wang, Chi and Liu, Xueqing and Song, Yanglei and Han, Jiawei},
  booktitle={ECMLPKDD},
  pages={290-305},
  year={2014},
}

**********************************************
Code explanation
**********************************************

(1). Input is in folder /Data. It contains two files: Data/test.corpus, which
 is document-word file (each line is in the format 'docID wordID') and
 Data/test.dict, which is vocabulary file (each line is in the format 'word')

(2). Output is in folder /Data. It contains two file:

     a) Data/test.stod.mat. It contains 3 parts:

     -- 1. alpha0, which is the summation of Dirichlet priors

     -- 2. alpha, which is the inferred Dirichlet prior of each topic

     -- 3. twmat, which is a k (number of topics) x V (vocabulary size) matrix,
	 the i-th row being the i-th inferred topic distribution p(w|topic=i)

     b) Data/test.topic.mat. It contains 1 part, a k x 1 cell, the i-th cell
	 being the topic representation of the i-th topic (which is the top 10
	 words ordered by p(w|topic=i))


**********************************************
Data explanation
**********************************************

Data/csabstract contains CS abstracts used in the paper.

For AP news dataset, please visit http://trec.nist.gov/data/docs_eng.html.


**** For More Questions ****

Please contact illimine.cs.illinois.edu or Chi Wang (chiw@microsoft.com)

