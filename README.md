# Analysis Report 1 for CS640
# Using amplicon data to understand the microbiome
## Due as a pull request on Monday, November 6th, 2017 before 4:30 pm

The **goal of this assignment** is to have you dig even deeper into the Fierer *et al.* dataset, and build your skills in R, markdown, scientific writing, and exploratory data analysis.

The data we will be using is from the NCBI Sequence Read Archive study number ERP022657. A summary of the information is available [here](https://www.ncbi.nlm.nih.gov/Traces/study/?WebEnv=NCID_1_128047291_130.14.22.33_5555_1505945515_1626731749_0MetA0_S_HStore&query_key=5). The metadata from this study is included in the git repository for this assignment in a `data/metadata` directory.

Here is the abstract from the [original study](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?study=ERP022657) by Fierer et al.:

> Recent work has demonstrated that the diversity of skin-associated bacterial communities is far higher than previously recognized, with a high degree of interindividual variability in the composition of bacterial communities. Given that skin bacterial communities are personalized, we hypothesized that we could use the residual skin bacteria left on objects for forensic identification, matching the bacteria on the object to the skin-associated bacteria of the individual who touched the object. Here we describe a series of studies demonstrating the validity of this approach. We show that skin-associated bacteria can be readily recovered from surfaces (including single computer keys and computer mice) and that the structure of these communities can be used to differentiate objects handled by different individuals, even if those objects have been left untouched for up to 2 weeks at room temperature. Furthermore, we demonstrate that we can use a high-throughput pyrosequencing-based approach to quantitatively compare the bacterial communities on objects and skin to match the object to the individual with a high degree of certainty. Although additional work is needed to further establish the utility of this approach, this series of studies introduces a forensics approach that could eventually be used to independently evaluate results obtained using more traditional forensic practices.

For this assignment, you will only need to work in one file, the `Rmd` file entitled `Analysis_Report_01_amplicons.Rmd`. There is a good deal of starter code in that file, with comments on how the various functions work. You should modify this document as needed, and add code chunks and markdown text to complete the following set of tasks.

This assignment very much builds on the work you did several weeks ago with the BLAST script and your analysis of the results. In this assignment we are digging deeper into the results using R. Just like last time, you should use your results as a launching off point to do some additional research on the taxa you find. You can use Google Scholar or some other tool to search the peer-reviewed literature. You are now responsible for finding and citing at least **five** additional peer-reviewed studies in your analysis report. 

This report could focus in on, for example: different bacterial taxa on male vs female hands, or the similarity of mouse surfaces to student hands, or just what is found on computer mice, etc. It matters less to me what specific aspect of the data you focus on, and more important that you spend some time digging in to the data and developing and testing some hypotheses of your own about what is going on in this dataset.

Please follow the instructions carefully and read them all before getting started.

This first analysis report will be worth 50 points and will count for 12.5% of your final semester grade. So take this assignment seriously and start early!! 

The grading breakdown will be as follows:

* 25 points - The report describes an analysis of the focal microbial dataset -- the analysis should not just be a set of unrelated figures, but should be a focused inquiry into some particular question or questions you have about the data. It should include appropriate background to motivate the investigation, one or more hypotheses, your quantitative results (including somewhere around 4-5 tables/figures beyond what is included already in the template), and a thorough discussion where you interpret your results. Report text should be: between 2500 and 3000 words (not counting code), well organized, free of spelling or grammar errors, and written in proper scientific writing style. You should use the built-in spell check in RStudio (there's a little button with ABC and a check mark near the knit button) to catch any spelling errors.
* 5 points - At least 5 sources properly cited in the text of your report beyond the two examples I've already included. You should add proper BibTeX entries to the `references.bib` file. You can get these BibTeX reference entries from any bibliography management software or from Google Scholar. Then, once you've added each of them and given them a citekey on the first line of each reference, you can use the RStudio add-in from the `citr` package to cite these within your text as you're writing. Then when you knit the document, the sources will be properly formatted for you, and the Sources Cited will be properly populated for you at the end of the report. Doing this will be/was covered in class.
* 5 points - R code chunks are appropriately commented, individually named, and well organized.
* 5 points - Appropriate use of git to version control the steps, including adding files and making commits frequently as you work on this assignment, and writing informative and appropriately formatted commit messages. Can't have too many commits, but you can have too few! 
* 10 points - Pull Request passes automated checks for file being able to be sucessfully knitted on another machine, as well as having all code style errors fixed. This is an all or nothing set of points, so please make sure your report passes these checks! You can submit your PR early to catch errors. Contact me well before the deadline if you are having trouble with this part of the assignment.

You must submit your work as a Pull Request to the class organization ('2017-usfca-cs-640') on GitHub by 4:30 pm on Monday, November 6th for full credit. We will also be peer reviewing the reports after they are submitted, as usual.

Steps:

1. Fork this repository to your own GitHub account.
1. In a new RStudio session on your laptop (no open files, projects, etc), select "File -> New Project..." from the menu. Then click on "Version Control", then "Git", and then paste the URL from your **forked** repository (the one under your GitHub account, that should include something like YourName/yourid-analysis-report-01-amplicons) into the "Repository URL" field. Then choose where you would like the folder to be on your computer and click "Create Project". This will use `git clone` to download the repository from GitHub and get RStudio set up for you to work on it. This should also properly set up the git remote so that you can push and pull from within RStudio. 
1. Run through the existing code in the document manually to make sure each step works and to build the combined dataset that you will use for analysis. What I mean by this is to execute the code chunks by hand instead of knitting the document. You may have to install libraries to your machine if you do not have them already.
1. You should include in this analysis report all of the parts of a standard scientific paper, including at a minimum: Introduction, Methods, Results, Discussion. I have added top-level headers by using a single `#`; if you want to make additional subsection headers you can use `##` or `###`. Each of these sections has some rough guidelines for how long I estimate they should be when completed and other details. As you are working on creating your analyses/tables/figures in particular, you will probably want to make good use of the [DADA2](https://benjjneb.github.io/dada2/tutorial.html) and [phyloseq](https://joey711.github.io/phyloseq/) websites, which include very comprehensive tutorials available.
1. Commit the script as you work on it, whenever you make a good chunk of progress. Make sure you write an [appropriate commit message](https://chris.beams.io/posts/git-commit/).
1. After you have finished the analysis, your `Rmd` knits successfully, and you have fixed all errors or warnings (except the 'there is no variable with this name in scope' warnings you see when using `dplyr`, which is a bug that hasn't been fixed by RStudio), be sure to add a commit marking this milestone (and push to GitHub!).
1. Once that's all done, add, commit, and push the rendered `md` file and the folder of images that are generated when you knit. These should go back to your fork of the original repository on GitHub. Remember that you can only push what you have committed, so be sure all of your work is committed. Be sure to save your files often, and check the git tab in RStudio as you work.
2. Submit a Pull Request back to the organization repository to submit your assignment. Make sure the Pull Request (PR) has a useful description of the changes you made. 20% of your points this week will be based on your code passing the automated checks on Travis-CI. If you want to check for any style errors while you are working on the code, you can install the 'lintr' package by typing the following at the console:

```r
install.packages("devtools")
devtools::install_github("jimhester/lintr")
```

and then whenever you want to check your code you can run at the R console in RStudio:

```r
lint(filename = "Analysis_Report_01_amplicons.Rmd")
```

**Pro Tip** Save often, commit often, push often! If you have any questions or need clarification of what it is I'd like you to do, please ask me sooner rather than later so you stay on the right track for completing this assignment on time.

**NOTE:** There are some files in this repository you may not have seen before. They are used for automating testing. 

##### Infrastructure for Automated Software Testing

- `.travis.yml`: a configuration file for automatically running [continuous integration](https://travis-ci.com) checks when you submit your pull request, to verify reproducibility of all `.Rmd` notebooks in the repo.  If all `.Rmd` notebooks can render successfully and pass linting (or code style and syntax checks), then the "Build Status" badge above will be green (`build success`), otherwise it will be red (`build failure`).  
- `.Rbuildignore`: a configuration file telling the build script which files to ignore.
- `DESCRIPTION`: a metadata file for the repository, based on the R package standard. Its main purpose here is as a place to list any additional R packages/libraries needed for any of the `.Rmd` files to run.
- `tests/render_rmds.R`: an R script that is run to execute the above described tests, rendering and linting all `.Rmd` notebooks. 

