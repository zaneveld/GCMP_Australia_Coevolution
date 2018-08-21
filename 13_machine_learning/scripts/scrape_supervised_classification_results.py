from sys import argv
from os import walk
from os.path import join,abspath

"""
A quick and dirty script to compile the results of many QIIME supervised_classification.py
results into a single table

"""
def parse_supervised_classifier_summary_output(lines):
    """Return a dict from QIIME supervised_classification.py output"""
    result = {}
    for line in lines:
        field,value = line.strip().split("\t")
        result[field] = value
    return result

def parse_supervised_classifier_folder(dirpath):
    """Get results from a supervised_classifier.py output directory"""
    summary_filepath = join(str(dirpath),"summary.txt")
    summary_file = open(summary_filepath,"U")
    summary_result =\
       parse_supervised_classifier_summary_output(summary_file.readlines()) 
    return summary_result

def main(input_folder,summary_table_fields=['Model','Error type','Estimated error','Baseline error (for random guessing)',\
    'Ratio baseline error to observed error','Number of trees']):
    """Output a summary table of machine learning results"""
    required_files = ['summary.txt','mislabeling.txt','feature_importance_scores.txt','cv_probabilities.txt','confusion_matrix.txt'] 
    
    result_lines = ["\t".join(summary_table_fields)+"\n"]
    for root,subdirectories,files in walk(input_folder):
        print root,subdirectories,files
        is_machine_learning_dir = True
        for f in required_files:
            if f not in files:
                print "Couldn't find required file %s in files %s...skipping" %(f,str(files))
                is_machine_learning_dir = False
                break

        if not is_machine_learning_dir:
            continue
            
        print "Found a machine learning directory: %s" % root
        results = parse_supervised_classifier_folder(root)
        curr_result = "\t".join([root]+[results[curr_field] for curr_field in summary_table_fields])
        result_lines.append(curr_result+"\n")
    
    return result_lines
        
    
 
if __name__ == "__main__":
    if len(argv) < 3: 
        raise ValueError("Specify an input folder as the first parameter and output file as the second")
    input_folder = argv[1]
    output = argv[2]
    result_lines = main(input_folder)
    f = open(output,"w+")
    f.writelines(result_lines)
    f.close()
