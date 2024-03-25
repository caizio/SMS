#include "dataDir.h"

// 0-S. pombe, 1-C. elegans, 2-S. cerevisiae, 3-A. thaliana, 4-D. Melanogaster, 5-HuRI
std::vector<std::string> LNBL3DataSet = {
    // 裂殖酵母
    "../dataset/LNBL3DataSet/S. pombe/S. pombe.txt",
    // 秀丽线虫
    "../dataset/LNBL3DataSet/C. elegans/C. elegans.txt",
    // 酿酒酵母
    "../dataset/LNBL3DataSet/S. cerevisiae/S. cerevisiae.txt",
    // 拟南芥
    "../dataset/LNBL3DataSet/A. thaliana/A. thaliana.txt",
    // 黑腹果蝇
    "../dataset/LNBL3DataSet/D. Melanogaster/D. Melanogaster.txt",
    // HuRI
    "../dataset/LNBL3DataSet/HuRI/HuRI.txt"
};

std::vector<std::string> seq = {
    "../dataset/seq/seq_S. pombeMap2.txt",
    "../dataset/seq/seq_C. elegansMap2.txt",
    "../dataset/seq/seq_S. cerevisiaeMap2.txt",
    "../dataset/seq/seq_A. thalianaMap2.txt",
    "../dataset/seq/seq_D. MelanogasterMap2.txt",
    "../dataset/seq/seq_HuRIMap2.txt"
};

std::vector<std::string>  datasetNames = {
    "S. pombe","C. elegans","S. cerevisiae", "A. thaliana", "D. Melanogaster","HuRI"
};
