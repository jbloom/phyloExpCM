INTEGRATION_PRECISION_FACTOR = 5.0e-6;
END_OF_FILE = 0;
LIKELIHOOD_FUNCTION_OUTPUT = 5;
ACCEPT_BRANCH_LENGTHS = 1;
#include "/home/jbloom/.local/lib/python2.7/site-packages/phyloExpCM/data//NTsCodonsAAs.ibf";
fprintf(stdout, "Running HYPHY script hyphy_cmds.bf...\n");
DataSet data = ReadDataFile("_codenames_aligned_KOSI07_TEM.fasta");
assert(data.sites % 3 == 0, "Sequence lengths not multiples of 3");
totalcodons = data.sites $ 3;
fprintf(stdout, "Read from _codenames_aligned_KOSI07_TEM.fasta a set of ", data.species, " sequences consisting of ", data.sites, " nucleotides corresponding to ", totalcodons, " codons each.\n");
fprintf(stdout, "The analysis will include the following 286 codon positions (sequential numbering starting with 1):\n1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286\n");
assert(totalcodons >= 286, "Largest included site exceeds sequence length");
DataSetFilter codonfilter = CreateFilter(data, 3, "0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,398,399,400,401,402,403,404,405,406,407,408,409,410,411,412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428,429,430,431,432,433,434,435,436,437,438,439,440,441,442,443,444,445,446,447,448,449,450,451,452,453,454,455,456,457,458,459,460,461,462,463,464,465,466,467,468,469,470,471,472,473,474,475,476,477,478,479,480,481,482,483,484,485,486,487,488,489,490,491,492,493,494,495,496,497,498,499,500,501,502,503,504,505,506,507,508,509,510,511,512,513,514,515,516,517,518,519,520,521,522,523,524,525,526,527,528,529,530,531,532,533,534,535,536,537,538,539,540,541,542,543,544,545,546,547,548,549,550,551,552,553,554,555,556,557,558,559,560,561,562,563,564,565,566,567,568,569,570,571,572,573,574,575,576,577,578,579,580,581,582,583,584,585,586,587,588,589,590,591,592,593,594,595,596,597,598,599,600,601,602,603,604,605,606,607,608,609,610,611,612,613,614,615,616,617,618,619,620,621,622,623,624,625,626,627,628,629,630,631,632,633,634,635,636,637,638,639,640,641,642,643,644,645,646,647,648,649,650,651,652,653,654,655,656,657,658,659,660,661,662,663,664,665,666,667,668,669,670,671,672,673,674,675,676,677,678,679,680,681,682,683,684,685,686,687,688,689,690,691,692,693,694,695,696,697,698,699,700,701,702,703,704,705,706,707,708,709,710,711,712,713,714,715,716,717,718,719,720,721,722,723,724,725,726,727,728,729,730,731,732,733,734,735,736,737,738,739,740,741,742,743,744,745,746,747,748,749,750,751,752,753,754,755,756,757,758,759,760,761,762,763,764,765,766,767,768,769,770,771,772,773,774,775,776,777,778,779,780,781,782,783,784,785,786,787,788,789,790,791,792,793,794,795,796,797,798,799,800,801,802,803,804,805,806,807,808,809,810,811,812,813,814,815,816,817,818,819,820,821,822,823,824,825,826,827,828,829,830,831,832,833,834,835,836,837,838,839,840,841,842,843,844,845,846,847,848,849,850,851,852,853,854,855,856,857", "", "TAA,TAG,TGA");
assert(data.species == codonfilter.species, "species number mismatch");
assert(codonfilter.sites == 286, "Codon filtered data does not contain the right number of sites");
fprintf(stdout, "Created a codon filter of ", codonfilter.sites, " sites.\n");
assert(totalcodons - (totalcodons - 286) - 0 == codonfilter.sites, "Codon filtered data is not the expected length. Do sequences contain stop codons?");
CheckCodonFilter("codonfilter");
fprintf(stdout, "Reading tree string from _codenames_codonphyml_tree_TEM.newick.\n");
fscanf("_codenames_codonphyml_tree_TEM.newick", String, treestring);
fprintf(stdout, "Using the Kosiol et al 2007 (KOSI07) codon model...\n");
#include "/home/jbloom/.local/lib/python2.7/site-packages/phyloExpCM/data//KOSI07.ibf";
CreateKOSI07Model("F", "ktv", "global", 4, 4, 1, "/home/jbloom/.local/lib/python2.7/site-packages/phyloExpCM/data//KOSI07_exchangeabilities.ibf");
UseModel(model);
ExecuteCommands("Tree tree = treestring;")
assert(codonfilter.species == TipCount(tree), "Number of species and number of tips differ");
LikelihoodFunction likelihood = (codonfilter, tree);
fprintf(stdout, "\nNow optimizing the likelihood function...\n");
Optimize(mlestimates, likelihood)
fprintf(stdout, "Completed likelihood optimization. Optimized ", mlestimates[1][1], " indpendent parameters and ", mlestimates[1][2], " shared parameters to obtain a log likelihood of ", mlestimates[1][0], ".\n");
fprintf(stdout, "Writing the results to hyphy_output.txt.\n");
fprintf("hyphy_output.txt", "Log likelihood: ", mlestimates[1][0], "\nindependent parameters (includes branch lengths): ", mlestimates[1][1], "\nshared parameters: ", mlestimates[1][2], "\nnumber of branch lengths: ", TipCount(tree) + BranchCount(tree), "\nnumber of tip nodes: ", TipCount(tree), "\nnumber of internal branches: ", BranchCount(tree), "\n",likelihood);
fprintf(stdout, "\nNow computing per-site likelihoods.\n");
fprintf(stdout, "\nFirst fixing all global variables to the maximum-likelihood values estimated on the entire tree.\n");
GetString(associativearray, likelihood, -1);
globalindependentvariables = associativearray["Global Independent"];
for (ivariable=0; ivariable<Columns(globalindependentvariables); ivariable=ivariable+1) {
  variable = globalindependentvariables[ivariable];
  cmdstring = variable + " := " + Format(variable, 0, 30) + ";";
  fprintf(stdout, "\nFixing variable as follows: ", cmdstring, "\n");
  ExecuteCommands(cmdstring);
}
persitelikelihoods = "sitelikelihoods.txt";
fprintf(stdout, "\nNow computing per-site likelihoods and writing to ", persitelikelihoods, "...\n");
fprintf(persitelikelihoods, "#SITE\tSITE_LOG_LIKELIHOODS\n");
fprintf(stdout, "\nComputing likelihood for site 1...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "0,1,2", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 1");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "1\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 2...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "3,4,5", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 2");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "2\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 3...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "6,7,8", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 3");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "3\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 4...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "9,10,11", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 4");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "4\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 5...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "12,13,14", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 5");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "5\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 6...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "15,16,17", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 6");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "6\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 7...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "18,19,20", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 7");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "7\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 8...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "21,22,23", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 8");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "8\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 9...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "24,25,26", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 9");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "9\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 10...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "27,28,29", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 10");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "10\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 11...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "30,31,32", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 11");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "11\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 12...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "33,34,35", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 12");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "12\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 13...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "36,37,38", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 13");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "13\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 14...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "39,40,41", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 14");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "14\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 15...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "42,43,44", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 15");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "15\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 16...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "45,46,47", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 16");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "16\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 17...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "48,49,50", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 17");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "17\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 18...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "51,52,53", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 18");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "18\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 19...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "54,55,56", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 19");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "19\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 20...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "57,58,59", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 20");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "20\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 21...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "60,61,62", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 21");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "21\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 22...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "63,64,65", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 22");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "22\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 23...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "66,67,68", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 23");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "23\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 24...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "69,70,71", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 24");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "24\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 25...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "72,73,74", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 25");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "25\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 26...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "75,76,77", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 26");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "26\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 27...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "78,79,80", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 27");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "27\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 28...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "81,82,83", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 28");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "28\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 29...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "84,85,86", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 29");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "29\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 30...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "87,88,89", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 30");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "30\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 31...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "90,91,92", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 31");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "31\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 32...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "93,94,95", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 32");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "32\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 33...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "96,97,98", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 33");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "33\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 34...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "99,100,101", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 34");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "34\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 35...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "102,103,104", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 35");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "35\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 36...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "105,106,107", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 36");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "36\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 37...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "108,109,110", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 37");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "37\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 38...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "111,112,113", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 38");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "38\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 39...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "114,115,116", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 39");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "39\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 40...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "117,118,119", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 40");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "40\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 41...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "120,121,122", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 41");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "41\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 42...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "123,124,125", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 42");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "42\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 43...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "126,127,128", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 43");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "43\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 44...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "129,130,131", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 44");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "44\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 45...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "132,133,134", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 45");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "45\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 46...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "135,136,137", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 46");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "46\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 47...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "138,139,140", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 47");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "47\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 48...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "141,142,143", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 48");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "48\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 49...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "144,145,146", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 49");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "49\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 50...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "147,148,149", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 50");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "50\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 51...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "150,151,152", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 51");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "51\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 52...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "153,154,155", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 52");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "52\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 53...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "156,157,158", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 53");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "53\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 54...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "159,160,161", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 54");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "54\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 55...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "162,163,164", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 55");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "55\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 56...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "165,166,167", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 56");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "56\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 57...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "168,169,170", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 57");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "57\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 58...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "171,172,173", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 58");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "58\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 59...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "174,175,176", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 59");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "59\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 60...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "177,178,179", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 60");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "60\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 61...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "180,181,182", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 61");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "61\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 62...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "183,184,185", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 62");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "62\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 63...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "186,187,188", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 63");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "63\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 64...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "189,190,191", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 64");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "64\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 65...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "192,193,194", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 65");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "65\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 66...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "195,196,197", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 66");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "66\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 67...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "198,199,200", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 67");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "67\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 68...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "201,202,203", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 68");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "68\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 69...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "204,205,206", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 69");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "69\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 70...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "207,208,209", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 70");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "70\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 71...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "210,211,212", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 71");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "71\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 72...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "213,214,215", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 72");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "72\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 73...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "216,217,218", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 73");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "73\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 74...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "219,220,221", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 74");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "74\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 75...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "222,223,224", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 75");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "75\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 76...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "225,226,227", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 76");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "76\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 77...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "228,229,230", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 77");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "77\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 78...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "231,232,233", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 78");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "78\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 79...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "234,235,236", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 79");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "79\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 80...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "237,238,239", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 80");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "80\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 81...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "240,241,242", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 81");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "81\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 82...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "243,244,245", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 82");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "82\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 83...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "246,247,248", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 83");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "83\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 84...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "249,250,251", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 84");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "84\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 85...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "252,253,254", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 85");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "85\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 86...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "255,256,257", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 86");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "86\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 87...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "258,259,260", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 87");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "87\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 88...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "261,262,263", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 88");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "88\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 89...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "264,265,266", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 89");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "89\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 90...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "267,268,269", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 90");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "90\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 91...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "270,271,272", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 91");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "91\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 92...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "273,274,275", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 92");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "92\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 93...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "276,277,278", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 93");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "93\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 94...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "279,280,281", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 94");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "94\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 95...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "282,283,284", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 95");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "95\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 96...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "285,286,287", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 96");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "96\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 97...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "288,289,290", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 97");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "97\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 98...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "291,292,293", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 98");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "98\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 99...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "294,295,296", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 99");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "99\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 100...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "297,298,299", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 100");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "100\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 101...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "300,301,302", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 101");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "101\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 102...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "303,304,305", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 102");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "102\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 103...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "306,307,308", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 103");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "103\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 104...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "309,310,311", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 104");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "104\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 105...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "312,313,314", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 105");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "105\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 106...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "315,316,317", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 106");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "106\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 107...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "318,319,320", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 107");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "107\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 108...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "321,322,323", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 108");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "108\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 109...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "324,325,326", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 109");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "109\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 110...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "327,328,329", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 110");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "110\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 111...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "330,331,332", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 111");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "111\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 112...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "333,334,335", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 112");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "112\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 113...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "336,337,338", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 113");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "113\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 114...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "339,340,341", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 114");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "114\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 115...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "342,343,344", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 115");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "115\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 116...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "345,346,347", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 116");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "116\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 117...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "348,349,350", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 117");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "117\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 118...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "351,352,353", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 118");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "118\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 119...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "354,355,356", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 119");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "119\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 120...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "357,358,359", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 120");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "120\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 121...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "360,361,362", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 121");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "121\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 122...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "363,364,365", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 122");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "122\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 123...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "366,367,368", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 123");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "123\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 124...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "369,370,371", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 124");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "124\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 125...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "372,373,374", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 125");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "125\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 126...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "375,376,377", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 126");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "126\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 127...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "378,379,380", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 127");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "127\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 128...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "381,382,383", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 128");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "128\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 129...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "384,385,386", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 129");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "129\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 130...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "387,388,389", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 130");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "130\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 131...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "390,391,392", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 131");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "131\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 132...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "393,394,395", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 132");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "132\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 133...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "396,397,398", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 133");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "133\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 134...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "399,400,401", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 134");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "134\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 135...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "402,403,404", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 135");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "135\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 136...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "405,406,407", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 136");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "136\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 137...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "408,409,410", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 137");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "137\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 138...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "411,412,413", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 138");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "138\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 139...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "414,415,416", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 139");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "139\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 140...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "417,418,419", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 140");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "140\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 141...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "420,421,422", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 141");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "141\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 142...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "423,424,425", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 142");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "142\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 143...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "426,427,428", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 143");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "143\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 144...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "429,430,431", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 144");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "144\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 145...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "432,433,434", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 145");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "145\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 146...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "435,436,437", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 146");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "146\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 147...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "438,439,440", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 147");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "147\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 148...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "441,442,443", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 148");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "148\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 149...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "444,445,446", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 149");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "149\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 150...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "447,448,449", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 150");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "150\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 151...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "450,451,452", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 151");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "151\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 152...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "453,454,455", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 152");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "152\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 153...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "456,457,458", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 153");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "153\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 154...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "459,460,461", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 154");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "154\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 155...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "462,463,464", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 155");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "155\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 156...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "465,466,467", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 156");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "156\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 157...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "468,469,470", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 157");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "157\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 158...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "471,472,473", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 158");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "158\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 159...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "474,475,476", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 159");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "159\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 160...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "477,478,479", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 160");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "160\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 161...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "480,481,482", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 161");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "161\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 162...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "483,484,485", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 162");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "162\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 163...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "486,487,488", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 163");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "163\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 164...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "489,490,491", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 164");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "164\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 165...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "492,493,494", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 165");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "165\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 166...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "495,496,497", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 166");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "166\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 167...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "498,499,500", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 167");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "167\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 168...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "501,502,503", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 168");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "168\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 169...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "504,505,506", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 169");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "169\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 170...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "507,508,509", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 170");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "170\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 171...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "510,511,512", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 171");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "171\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 172...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "513,514,515", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 172");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "172\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 173...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "516,517,518", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 173");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "173\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 174...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "519,520,521", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 174");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "174\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 175...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "522,523,524", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 175");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "175\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 176...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "525,526,527", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 176");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "176\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 177...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "528,529,530", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 177");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "177\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 178...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "531,532,533", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 178");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "178\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 179...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "534,535,536", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 179");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "179\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 180...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "537,538,539", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 180");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "180\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 181...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "540,541,542", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 181");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "181\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 182...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "543,544,545", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 182");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "182\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 183...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "546,547,548", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 183");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "183\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 184...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "549,550,551", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 184");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "184\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 185...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "552,553,554", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 185");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "185\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 186...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "555,556,557", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 186");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "186\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 187...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "558,559,560", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 187");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "187\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 188...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "561,562,563", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 188");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "188\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 189...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "564,565,566", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 189");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "189\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 190...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "567,568,569", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 190");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "190\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 191...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "570,571,572", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 191");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "191\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 192...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "573,574,575", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 192");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "192\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 193...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "576,577,578", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 193");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "193\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 194...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "579,580,581", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 194");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "194\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 195...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "582,583,584", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 195");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "195\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 196...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "585,586,587", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 196");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "196\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 197...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "588,589,590", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 197");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "197\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 198...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "591,592,593", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 198");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "198\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 199...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "594,595,596", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 199");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "199\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 200...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "597,598,599", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 200");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "200\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 201...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "600,601,602", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 201");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "201\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 202...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "603,604,605", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 202");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "202\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 203...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "606,607,608", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 203");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "203\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 204...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "609,610,611", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 204");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "204\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 205...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "612,613,614", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 205");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "205\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 206...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "615,616,617", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 206");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "206\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 207...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "618,619,620", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 207");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "207\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 208...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "621,622,623", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 208");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "208\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 209...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "624,625,626", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 209");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "209\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 210...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "627,628,629", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 210");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "210\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 211...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "630,631,632", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 211");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "211\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 212...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "633,634,635", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 212");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "212\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 213...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "636,637,638", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 213");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "213\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 214...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "639,640,641", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 214");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "214\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 215...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "642,643,644", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 215");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "215\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 216...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "645,646,647", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 216");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "216\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 217...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "648,649,650", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 217");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "217\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 218...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "651,652,653", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 218");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "218\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 219...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "654,655,656", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 219");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "219\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 220...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "657,658,659", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 220");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "220\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 221...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "660,661,662", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 221");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "221\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 222...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "663,664,665", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 222");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "222\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 223...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "666,667,668", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 223");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "223\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 224...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "669,670,671", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 224");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "224\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 225...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "672,673,674", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 225");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "225\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 226...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "675,676,677", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 226");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "226\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 227...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "678,679,680", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 227");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "227\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 228...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "681,682,683", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 228");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "228\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 229...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "684,685,686", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 229");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "229\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 230...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "687,688,689", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 230");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "230\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 231...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "690,691,692", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 231");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "231\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 232...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "693,694,695", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 232");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "232\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 233...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "696,697,698", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 233");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "233\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 234...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "699,700,701", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 234");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "234\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 235...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "702,703,704", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 235");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "235\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 236...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "705,706,707", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 236");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "236\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 237...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "708,709,710", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 237");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "237\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 238...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "711,712,713", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 238");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "238\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 239...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "714,715,716", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 239");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "239\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 240...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "717,718,719", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 240");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "240\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 241...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "720,721,722", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 241");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "241\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 242...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "723,724,725", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 242");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "242\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 243...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "726,727,728", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 243");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "243\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 244...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "729,730,731", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 244");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "244\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 245...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "732,733,734", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 245");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "245\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 246...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "735,736,737", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 246");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "246\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 247...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "738,739,740", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 247");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "247\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 248...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "741,742,743", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 248");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "248\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 249...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "744,745,746", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 249");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "249\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 250...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "747,748,749", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 250");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "250\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 251...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "750,751,752", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 251");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "251\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 252...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "753,754,755", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 252");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "252\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 253...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "756,757,758", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 253");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "253\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 254...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "759,760,761", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 254");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "254\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 255...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "762,763,764", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 255");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "255\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 256...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "765,766,767", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 256");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "256\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 257...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "768,769,770", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 257");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "257\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 258...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "771,772,773", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 258");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "258\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 259...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "774,775,776", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 259");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "259\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 260...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "777,778,779", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 260");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "260\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 261...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "780,781,782", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 261");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "261\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 262...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "783,784,785", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 262");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "262\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 263...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "786,787,788", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 263");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "263\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 264...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "789,790,791", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 264");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "264\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 265...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "792,793,794", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 265");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "265\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 266...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "795,796,797", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 266");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "266\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 267...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "798,799,800", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 267");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "267\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 268...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "801,802,803", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 268");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "268\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 269...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "804,805,806", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 269");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "269\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 270...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "807,808,809", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 270");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "270\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 271...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "810,811,812", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 271");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "271\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 272...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "813,814,815", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 272");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "272\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 273...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "816,817,818", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 273");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "273\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 274...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "819,820,821", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 274");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "274\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 275...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "822,823,824", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 275");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "275\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 276...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "825,826,827", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 276");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "276\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 277...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "828,829,830", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 277");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "277\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 278...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "831,832,833", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 278");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "278\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 279...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "834,835,836", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 279");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "279\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 280...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "837,838,839", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 280");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "280\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 281...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "840,841,842", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 281");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "281\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 282...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "843,844,845", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 282");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "282\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 283...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "846,847,848", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 283");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "283\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 284...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "849,850,851", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 284");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "284\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 285...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "852,853,854", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 285");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "285\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "\nComputing likelihood for site 286...\n");
DataSetFilter sitecodonfilter = CreateFilter(data, 3, "855,856,857", "", "TAA,TAG,TGA");
assert(sitecodonfilter.sites == 1, "Codon filtered data does not have one site for 286");
CheckCodonFilter("sitecodonfilter");
UseModel(model);
ExecuteCommands("Tree sitetree = treestring;");
assert(sitecodonfilter.species == TipCount(sitetree), "Number of species and number of tips differ");
assert(TipCount(tree) == TipCount(sitetree), "Number of tips differ");
for (ibranch=0; ibranch<BranchCount(tree); ibranch=ibranch+1) {
  branchname = BranchName(tree, ibranch);
  ExecuteCommands("branchlength = tree." + branchname + ".t;");
  ExecuteCommands("sitetree." + branchname + ".t := " + Format(branchlength, 0, 30) + ";");
}
for (itip=0; itip<TipCount(tree); itip=itip+1) {
  tipname = TipName(tree, itip);
  ExecuteCommands("tiplength = tree." + tipname + ".t;");
  ExecuteCommands("sitetree." + tipname + ".t := " + Format(tiplength, 0, 30) + ";");
}
LikelihoodFunction sitelikelihood = (sitecodonfilter, sitetree);
Optimize(sitemlestimates, sitelikelihood);
assert(sitemlestimates[1][1] == 0, "Found a variable optimized. Either a model or branch parameter must have not been fixed");
fprintf(persitelikelihoods, "286\t", sitemlestimates[1][0], "\n");
fprintf(stdout, "Completed HYPHY script hyphy_cmds.bf.\n");