
The 'bioclassifier_R' folder contains all the functions and fixed parameters.  You should not need to change anything in this folder, but there are some annotation files, the centroids file, and related information that you may want to be familiar with. 

The directory named 'bioclassifier_example' contains two files - an example input data matrix, and an example of the parameter file to run the algorithm.  All changes should be made to this R file, and I have labeled the variables accordingly.  As a first test, simply fill in the correct directory paths and source the file.  The output of this run should look like the set of files in the 'sampleOutput' sub-folder.

As you are probably aware, measurement bias and population bias are often confounded when we attempt cross-platform classification.  When the test sample draws from the same population as the training set, then the population bias is less of a problem and simple approaches to adjusting for measurement bias (gene centering) tend to work quite well.  However, as the test sample deviates from the training sample, the confounded effects become more difficult to overcome and allow for accurate classification.   An example of an extreme case would be a test sample that was all ER+.

One valuable resource in checking the pre-processing is if you have any true 'normal' samples.  When true normal samples are not called 'normal-like', then it is an indicator of bias.

You are likely more informed of such issues than myself, but I feel I must give this disclaimer whenever I provide the code.  I would be happy to discuss with you further regarding approaches for minimizing bias.  Please let me know if you are interested, and please let me know if you have trouble running or interpreting the code.

Joel Parker
parkerjs@email.unc.edu