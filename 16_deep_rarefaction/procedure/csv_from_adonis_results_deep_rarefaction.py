from __future__ import division
from sys import argv

from os import sep,walk
from os.path import abspath

def parse_results_dir(root_results_dir,sep="\t"):
    """"Starting at filepath, find Adonis results"""
    header_line = sep.join(["Compartment","category","bdiv_metric","n_permuations",\
      "function_call","df","SS","F","R2","R2_adj","R2_adj_rounded","p","sig","residual df","residual ss",\
      "residual r2","total df","file"])+"\n"

    result_lines = [header_line]
    for path in iterate_subdirs(root_results_dir):

        compartment,bdiv_metric,category = get_analysis_type_from_filepath(path)
        if compartment is None:
            continue
        category_results,residual_results,total_results,n_permutations,function_call = parse_adonis_results(path)

        data_line = format_data_line(compartment,bdiv_metric,category,category_results,\
          residual_results,total_results,n_permutations,function_call,path,sep=sep)

        result_lines.append(data_line)

    return result_lines

def format_data_line(compartment,bdiv_metric,category,category_results,\
  residual_results,total_results,n_permutations,function_call,path,sep="\t"):
    """Format a csv using results"""
    #Headings are:
    if category_results["p"] < 0.001:
        category_results["sig"] = '***'
    elif category_results["p"] < 0.01:
        category_results["sig"] = '**'
    elif category_results["p"] < 0.05:
        category_results["sig"] = "*"
    elif category_results["p"] < 0.1:
        category_results["sig"] = "."
    else:
        category_results["sig"] = " "

    #adjust R2
    df = category_results["degrees_of_freedom"]
    total_df = total_results["degrees_of_freedom"]
    r2 = category_results["R2"]
    adjusted_r2 = adjust_r2(r2,df,total_df)
    rounded_adjusted_r2 = round(adjusted_r2,2)

    data_list = [compartment,category,bdiv_metric,n_permutations,function_call,\
      category_results["degrees_of_freedom"],category_results["sum_of_squares"],category_results["F"],category_results["R2"],\
      adjusted_r2,rounded_adjusted_r2,category_results["p"],category_results["sig"],residual_results["degrees_of_freedom"],\
      residual_results["sum_of_squares"],residual_results["R2"],\
      total_results["degrees_of_freedom"],path]
    return sep.join(map(str,data_list))+"\n"

def iterate_subdirs(rootdir,extension="adonis_results.txt"):
    for subdir, dirs, files in walk(rootdir,topdown=True):
        for f in files:
            filepath = subdir + sep + f
            if filepath.endswith(extension):
                yield filepath

def adjust_r2(r2,df,total_df):
    """Calculate the adjusted R2

    1 - (((1-r2)*(n-1))/n-k-1)
    """
    k = df #number of parameters we want to fit
    n = total_df + 1.0 #number of samples
    numerator = (1.0-r2)*(n-1.0)
    denominator = n-k-1.0
    return 1 - (numerator/denominator)

def parse_adonis_results(filepath):
    """Return F,P,R2 from Adonis results file
    """
    function_call = None
    n_permutations = None
    category_results = {"factor":None,"degrees_of_freedom":None,"sum_of_squares":None,"F":None,"R2":None,"p":None}
    residual_results = {"factor":None,"degrees_of_freedom":None,"sum_of_squares":None,"F":None,"R2":None,"p":None}
    total_results = {"factor":None,"degrees_of_freedom":None,"sum_of_squares":None,"F":None,"R2":None,"p":None}

    for line in open(filepath,"U").readlines():
        if "adonis" in line:
            #have to strip the line
            #otherwise newline chars will mess up output format
            function_call = line.strip()
            continue
        if line.startswith("Number of permutations:"):
            n_permutations = int(line.split(":")[1])
            continue
        #Don't actually need to parse the header line as it's consistent
        if line.rstrip().startswith("qiime.data$map"):
            fields = line.split()
            if len(fields) == 7: #no sig marker
                factor,df,ss,ms,F,R2,p = line.split()
                sig_marker = ''
            elif len(fields) == 8:
                factor,df,ss,ms,F,R2,p,sig_marker = line.split()

            category_results["degrees_of_freedom"]=int(df)
            category_results["sum_of_squares"]=float(ss)
            category_results["mean_square"] = float(ms)
            category_results["F"]=float(F)
            category_results["sig_marker"]=sig_marker
            category_results["R2"]=float(R2)
            category_results["p"]=float(p)


        elif "Residuals" in line:
            factor,df,ss,ms,R2 = line.split()
            residual_results["degrees_of_freedom"]=int(df)
            residual_results["sum_of_squares"]=float(ss)
            residual_results["mean_square"] = float(ms)
            residual_results["R2"]=float(R2)

        elif "Total" in line:
            factor,df,ss,R2 = line.split()
            total_results["degrees_of_freedom"]=int(df)
            total_results["sum_of_squares"]=float(ss)
            total_results["R2"]=float(R2)
            #After the Total lines repeat themselves, only
            #if there were sig. results.
            #To avoid having to parse that, break here
            break
    return category_results,residual_results,\
      total_results,n_permutations,function_call

def get_analysis_type_from_filepath(filepath):
    """Return the analysis details from a consistently structured filepath
    MUST be an adonis_results.txt file resulting from compare_categories.py

    NOTE: Given a filepath ending with:
    bdiv_otu_table_subset_mucus.biom/
    compare_categories_visibility_weighted_unifrac/
    adonis_results.txt

    the script will figure out this is mucus, testing visibility.
    """
    subfolders = filepath.split(sep)
    filename = subfolders[-1] #the last part of the path
    analysis_type_dir = subfolders[-2] #next to last part
    compartment_dir = subfolders[-3]

    compartment = None
    for curr_compartment in ['mucus','tissue','skeleton']:
        if curr_compartment in compartment_dir.lower():
            compartment = curr_compartment
            break

    #extract bdiv metric, category and script from analysis_type_dir
    if not analysis_type_dir.startswith('compare_categories'):
        raise ValueError("Can only parse compare_categories results")

    analysis_type_dir = analysis_type_dir[len('compare_categories_'):] #chop off the beginning, effectively
    bdiv_metric = None
    #Note: a bit hackish: the order matters
    for curr_bdiv_metric in ['unweighted_unifrac','weighted_unifrac','bray_curtis']:
        if analysis_type_dir.endswith(curr_bdiv_metric):
            bdiv_metric = curr_bdiv_metric
            analysis_type_dir = analysis_type_dir[:-len(curr_bdiv_metric)-1] #chop off the end and underscore
            break
    #what remains should be the category!
    category = analysis_type_dir

    return compartment,bdiv_metric,category



def main():
    """Loop over the file structure, finding Adonis results"""
    if len(argv) < 3:

        print("Usage: root directory for results you wish to parse as an argument")
        print(\
            "Example: python csv_from_adonis_results.py ../../output/compartment_tables/bdiv_otu_table_subset_mucus.biom/ ./adonis_results.txt")
        raise ValueError("must specify input and output file and 1st and 2nd arguments")
    #set the separator for output files to "," for csv, "\t" for tsv
    sep = "\t"
    raw_results = parse_results_dir(abspath(argv[1]),sep=sep)

    outfile = argv[2]
    f=open(outfile,"w+")
    f.writelines(raw_results)
    f.close()
if __name__ == "__main__":
    main()
