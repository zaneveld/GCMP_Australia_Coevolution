
#Code from Antonio's GitHub issue response to detrend.py issue.
#I added a quick commandline interface

from sys import argv
from skbio.stats.ordination import OrdinationResults


# This function was directly taken from QIIME 1.8.0, tests and other
# information can be found there.
# https://github.com/biocore/qiime/blob/d4333e2ea06af942f1f61148c4ccb02ffc438d6b/qiime/format.py


def format_coords(coord_header, coords, eigvals, pct_var, headers = True):
    """formats coords given specified coords matrix etc."""
    result = []
    if (headers):
        result.append('pc vector number\t' +
           '\t'.join(map(str, range(1,len(coords[0])+1))))
        for name, row in zip(coord_header, coords):
            result.append('\t'.join([name] + map(str, row)))
        result.append('')
        result.append('')
        result.append('eigvals\t' + '\t'.join(map(str,eigvals)))
        result.append('% variation explained\t' +
           '\t'.join(map(str, pct_var)))
    else:
        result = ['\t'.join(map(str, row)) for row in coords]
        result.append('')
    return '\n'.join(result)


if __name__ == "__main__":
    old_file = argv[1]
    new_file = argv[2]
    with open(old_file, 'U') as infile:
        with open(new_file,'w') as outfile:
            res = OrdinationResults.from_file(infile)
            lines = format_coords(res.site_ids, res.site, res.eigvals, res.proportion_explained)
            outfile.write(lines)


